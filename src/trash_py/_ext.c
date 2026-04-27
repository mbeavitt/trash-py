#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// FNV-1a hash function
static uint64_t fnv1a(const char *p, Py_ssize_t n) {
    uint64_t h = 14695981039346656037ULL;
    for (Py_ssize_t i = 0; i < n; i++) {
        h ^= (uint8_t)p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

typedef struct {
    uint32_t *slots;
    Py_ssize_t *starts;
    size_t cap_mask;
    Py_ssize_t kmer_len;
    uint32_t next_id;
    const char *sequence;
} KmerIdMap;

static int id_map_init(KmerIdMap *m, uint32_t max_unique, Py_ssize_t kmer_len, const char *sequence) {
    size_t cap = 16;
    while (cap < (size_t)max_unique * 2u) cap <<= 1;
    m->slots = PyMem_Calloc(cap, sizeof(*m->slots));
    m->starts = PyMem_Malloc((size_t)max_unique * sizeof(*m->starts));
    if (!m->slots || !m->starts) {
        PyMem_Free(m->slots);
        PyMem_Free(m->starts);
        return 0;
    }
    m->cap_mask = cap - 1;
    m->kmer_len = kmer_len;
    m->next_id = 0;
    m->sequence = sequence;
    return 1;
}

static void id_map_free(KmerIdMap *m) {
    PyMem_Free(m->slots);
    PyMem_Free(m->starts);
}

static inline uint32_t id_map_intern(KmerIdMap *m, Py_ssize_t start) {
    size_t idx = (size_t)fnv1a(m->sequence + start, m->kmer_len) & m->cap_mask;
    while (m->slots[idx]) {
        uint32_t id = m->slots[idx] - 1u;
        if (memcmp(m->sequence + m->starts[id], m->sequence + start, (size_t)m->kmer_len) == 0) {
            return id;
        }
        idx = (idx + 1u) & m->cap_mask;
    }
    m->starts[m->next_id] = start;
    m->slots[idx] = m->next_id + 1u;
    return m->next_id++;
}

static int build_kmer_ids(
    const char *sequence,
    Py_ssize_t kmer_len,
    Py_ssize_t n_kmers,
    uint32_t **out_ids,
    uint32_t *out_unique_ids
) {
    KmerIdMap map;
    if (!id_map_init(&map, (uint32_t)n_kmers, kmer_len, sequence)) return 0;

    uint32_t *ids = PyMem_Malloc((size_t)n_kmers * sizeof(*ids));
    if (!ids) {
        id_map_free(&map);
        return 0;
    }

    for (Py_ssize_t i = 0; i < n_kmers; i++) {
        ids[i] = id_map_intern(&map, i);
    }

    *out_ids = ids;
    *out_unique_ids = map.next_id;
    id_map_free(&map);
    return 1;
}

static int parse_int_sequence(PyObject *obj, int32_t **out, Py_ssize_t *out_len, const char *name) {
    PyObject *seq = PySequence_Fast(obj, name);
    if (!seq) return 0;

    Py_ssize_t len = PySequence_Fast_GET_SIZE(seq);
    int32_t *values = PyMem_Malloc((size_t)len * sizeof(*values));
    if (!values) {
        Py_DECREF(seq);
        PyErr_NoMemory();
        return 0;
    }

    PyObject **items = PySequence_Fast_ITEMS(seq);
    for (Py_ssize_t i = 0; i < len; i++) {
        long value = PyLong_AsLong(items[i]);
        if (value == -1 && PyErr_Occurred()) {
            Py_DECREF(seq);
            PyMem_Free(values);
            return 0;
        }
        if (value < INT32_MIN || value > INT32_MAX) {
            Py_DECREF(seq);
            PyMem_Free(values);
            PyErr_Format(PyExc_OverflowError, "%s[%zd] is out of int32 range", name, i);
            return 0;
        }
        values[i] = (int32_t)value;
    }

    Py_DECREF(seq);
    *out = values;
    *out_len = len;
    return 1;
}

static int can_use_sliding_suffix(
    Py_ssize_t n_windows,
    int32_t n_kmers,
    const int32_t *ws,
    const int32_t *we,
    const int32_t *wsc,
    const int32_t *wec,
    uint32_t *out_step
) {
    if (n_windows == 0) {
        *out_step = 0;
        return 1;
    }

    int32_t a_len = we[0] - ws[0] + 1;
    if (a_len <= 0 || ws[0] < 1 || we[0] > n_kmers) return 0;
    if (wsc[0] != we[0] + 1 || wec[0] <= n_kmers) return 0;

    uint32_t step = 0;
    if (n_windows > 1) {
        int32_t delta = ws[1] - ws[0];
        if (delta <= 0 || delta > a_len) return 0;
        step = (uint32_t)delta;
    }

    for (Py_ssize_t i = 0; i < n_windows; i++) {
        if (ws[i] < 1 || we[i] > n_kmers) return 0;
        if (we[i] - ws[i] + 1 != a_len) return 0;
        if (wsc[i] != we[i] + 1) return 0;
        if (wec[i] <= n_kmers) return 0;
        if (i > 0) {
            if ((uint32_t)(ws[i] - ws[i - 1]) != step) return 0;
            if ((uint32_t)(we[i] - we[i - 1]) != step) return 0;
            if ((uint32_t)(wsc[i] - wsc[i - 1]) != step) return 0;
        }
    }

    *out_step = step;
    return 1;
}

static PyObject *build_zero_scores(Py_ssize_t n_windows, const int32_t *wsc, const int32_t *wec) {
    PyObject *scores = PyList_New(n_windows);
    if (!scores) return NULL;

    for (Py_ssize_t i = 0; i < n_windows; i++) {
        int32_t b_len = wec[i] - wsc[i] + 1;
        if (b_len <= 0) {
            Py_DECREF(scores);
            PyErr_SetString(PyExc_ZeroDivisionError, "comparison window has zero length");
            return NULL;
        }
        PyObject *value = PyFloat_FromDouble(0.0);
        if (!value) {
            Py_DECREF(scores);
            return NULL;
        }
        PyList_SET_ITEM(scores, i, value);
    }

    return scores;
}

static PyObject *run_generic_ids(
    const uint32_t *ids,
    uint32_t unique_ids,
    int32_t n_kmers,
    Py_ssize_t n_windows,
    const int32_t *ws,
    const int32_t *we,
    const int32_t *wsc,
    const int32_t *wec
) {
    if (unique_ids == 0) return build_zero_scores(n_windows, wsc, wec);

    uint32_t *marks = PyMem_Calloc((size_t)unique_ids, sizeof(*marks));
    if (!marks) {
        PyErr_NoMemory();
        return NULL;
    }

    PyObject *scores = PyList_New(n_windows);
    if (!scores) {
        PyMem_Free(marks);
        return NULL;
    }

    uint32_t epoch = 1;
    for (Py_ssize_t i = 0; i < n_windows; i++) {
        if (epoch == 0) {
            memset(marks, 0, (size_t)unique_ids * sizeof(*marks));
            epoch = 1;
        }

        int32_t b_len = wec[i] - wsc[i] + 1;
        if (b_len <= 0) {
            Py_DECREF(scores);
            PyMem_Free(marks);
            PyErr_SetString(PyExc_ZeroDivisionError, "comparison window has zero length");
            return NULL;
        }

        int32_t a0 = ws[i] < 1 ? 1 : ws[i];
        int32_t a1 = we[i] > n_kmers ? n_kmers : we[i];
        for (int32_t p = a0; p <= a1; p++) {
            marks[ids[p - 1]] = epoch;
        }

        int32_t b0 = wsc[i] < 1 ? 1 : wsc[i];
        int32_t b1 = wec[i] > n_kmers ? n_kmers : wec[i];
        int32_t count = 0;
        for (int32_t p = b0; p <= b1; p++) {
            count += (marks[ids[p - 1]] == epoch);
        }

        PyObject *value = PyFloat_FromDouble((double)count / (double)b_len);
        if (!value) {
            Py_DECREF(scores);
            PyMem_Free(marks);
            return NULL;
        }
        PyList_SET_ITEM(scores, i, value);
        epoch++;
    }

    PyMem_Free(marks);
    return scores;
}

static PyObject *run_sliding_suffix_ids(
    const uint32_t *ids,
    uint32_t unique_ids,
    int32_t n_kmers,
    Py_ssize_t n_windows,
    uint32_t step,
    const int32_t *ws,
    const int32_t *we,
    const int32_t *wsc,
    const int32_t *wec
) {
    if (unique_ids == 0) return build_zero_scores(n_windows, wsc, wec);

    uint32_t *window_counts = PyMem_Calloc((size_t)unique_ids, sizeof(*window_counts));
    uint32_t *suffix_counts = PyMem_Calloc((size_t)unique_ids, sizeof(*suffix_counts));
    if (!window_counts || !suffix_counts) {
        PyMem_Free(window_counts);
        PyMem_Free(suffix_counts);
        PyErr_NoMemory();
        return NULL;
    }

    PyObject *scores = PyList_New(n_windows);
    if (!scores) {
        PyMem_Free(window_counts);
        PyMem_Free(suffix_counts);
        return NULL;
    }

    for (int32_t p = wsc[0]; p <= n_kmers; p++) {
        suffix_counts[ids[p - 1]]++;
    }

    uint64_t total_match = 0;
    for (int32_t p = ws[0]; p <= we[0]; p++) {
        uint32_t id = ids[p - 1];
        if (window_counts[id]++ == 0u) total_match += suffix_counts[id];
    }

    for (Py_ssize_t i = 0; i < n_windows; i++) {
        int32_t b_len = wec[i] - wsc[i] + 1;
        if (b_len <= 0) {
            Py_DECREF(scores);
            PyMem_Free(window_counts);
            PyMem_Free(suffix_counts);
            PyErr_SetString(PyExc_ZeroDivisionError, "comparison window has zero length");
            return NULL;
        }

        PyObject *value = PyFloat_FromDouble((double)total_match / (double)b_len);
        if (!value) {
            Py_DECREF(scores);
            PyMem_Free(window_counts);
            PyMem_Free(suffix_counts);
            return NULL;
        }
        PyList_SET_ITEM(scores, i, value);

        if (i + 1 == n_windows) break;

        int32_t remove_a_stop = ws[i] + (int32_t)step;
        for (int32_t p = ws[i]; p < remove_a_stop; p++) {
            uint32_t id = ids[p - 1];
            if (--window_counts[id] == 0u) total_match -= suffix_counts[id];
        }

        int32_t add_a_stop = we[i] + (int32_t)step;
        for (int32_t p = we[i] + 1; p <= add_a_stop; p++) {
            uint32_t id = ids[p - 1];
            if (window_counts[id]++ == 0u) total_match += suffix_counts[id];
        }

        int32_t remove_b_stop = wsc[i] + (int32_t)step;
        for (int32_t p = wsc[i]; p < remove_b_stop; p++) {
            if (p > n_kmers) break;
            uint32_t id = ids[p - 1];
            suffix_counts[id]--;
            total_match -= (window_counts[id] != 0u);
        }
    }

    PyMem_Free(window_counts);
    PyMem_Free(suffix_counts);
    return scores;
}

static PyObject *seq_win_score_int_impl(PyObject *self, PyObject *args) {
    Py_ssize_t start, end;
    int kmer;
    PyObject *seq_obj;

    if (!PyArg_ParseTuple(args, "nniO", &start, &end, &kmer, &seq_obj)) {
        return NULL;
    }

    if (kmer <= 0) {
        PyErr_SetString(PyExc_ValueError, "kmer must be positive");
        return NULL;
    }

    if ((end - start) <= (Py_ssize_t)kmer) {
        return PyFloat_FromDouble(100.0);
    }

    PyObject *seq_bytes = PyUnicode_AsASCIIString(seq_obj);
    if (!seq_bytes) return NULL;

    Py_ssize_t seq_len = PyBytes_GET_SIZE(seq_bytes);
    const char *sequence = PyBytes_AS_STRING(seq_bytes);

    Py_ssize_t i_start = start - 1;
    Py_ssize_t i_end_excl = end - (Py_ssize_t)kmer;
    Py_ssize_t n_kmers = i_end_excl - i_start;

    if (i_start < 0 || i_end_excl + (Py_ssize_t)kmer - 1 > seq_len) {
        Py_DECREF(seq_bytes);
        PyErr_SetString(PyExc_ValueError, "window extends outside sequence bounds");
        return NULL;
    }
    if (n_kmers > (Py_ssize_t)UINT32_MAX) {
        Py_DECREF(seq_bytes);
        PyErr_SetString(PyExc_OverflowError, "too many kmers for the seq_win_score_int accelerator");
        return NULL;
    }

    KmerIdMap map;
    if (!id_map_init(&map, (uint32_t)n_kmers, (Py_ssize_t)kmer, sequence)) {
        Py_DECREF(seq_bytes);
        PyErr_NoMemory();
        return NULL;
    }

    uint32_t *counts = PyMem_Calloc((size_t)n_kmers, sizeof(*counts));
    if (!counts) {
        id_map_free(&map);
        Py_DECREF(seq_bytes);
        PyErr_NoMemory();
        return NULL;
    }

    for (Py_ssize_t i = i_start; i < i_end_excl; i++) {
        int has_n = 0;
        for (int j = 0; j < kmer; j++) {
            char c = sequence[i + j];
            if (c == 'n' || c == 'N') { has_n = 1; break; }
        }
        if (has_n) continue;
        uint32_t id = id_map_intern(&map, i);
        counts[id]++;
    }

    uint64_t total = 0;
    uint64_t singletons = 0;
    for (uint32_t id = 0; id < map.next_id; id++) {
        total += counts[id];
        if (counts[id] == 1u) singletons++;
    }

    id_map_free(&map);
    PyMem_Free(counts);
    Py_DECREF(seq_bytes);

    if (total < (uint64_t)kmer * 2u) {
        return PyFloat_FromDouble(100.0);
    }

    return PyFloat_FromDouble(100.0 * (double)singletons / (double)total);
}

typedef struct {
    const char *name;
    Py_ssize_t name_len;
    long count;
    Py_ssize_t orig_idx;
} CollapseEntry;

static int collapse_entry_cmp(const void *a, const void *b) {
    const CollapseEntry *ea = (const CollapseEntry *)a;
    const CollapseEntry *eb = (const CollapseEntry *)b;
    if (ea->count != eb->count) return ea->count < eb->count ? -1 : 1;
    Py_ssize_t min_len = ea->name_len < eb->name_len ? ea->name_len : eb->name_len;
    int c = memcmp(ea->name, eb->name, (size_t)min_len);
    if (c != 0) return c;
    if (ea->name_len != eb->name_len) return ea->name_len < eb->name_len ? -1 : 1;
    return 0;
}

static PyObject *collapse_kmers_impl(PyObject *self, PyObject *args) {
    PyObject *names_obj;
    PyObject *counts_obj;
    int max_edit;

    if (!PyArg_ParseTuple(args, "OOi", &names_obj, &counts_obj, &max_edit)) {
        return NULL;
    }
    if (max_edit < 0) {
        PyErr_SetString(PyExc_ValueError, "max_edit must be non-negative");
        return NULL;
    }

    PyObject *names_seq = PySequence_Fast(names_obj, "names must be a sequence");
    if (!names_seq) return NULL;
    PyObject *counts_seq = PySequence_Fast(counts_obj, "counts must be a sequence");
    if (!counts_seq) {
        Py_DECREF(names_seq);
        return NULL;
    }

    Py_ssize_t n = PySequence_Fast_GET_SIZE(names_seq);
    if (n != PySequence_Fast_GET_SIZE(counts_seq)) {
        Py_DECREF(names_seq);
        Py_DECREF(counts_seq);
        PyErr_SetString(PyExc_ValueError, "names and counts must have the same length");
        return NULL;
    }

    if (n == 0) {
        Py_DECREF(names_seq);
        Py_DECREF(counts_seq);
        return PyList_New(0);
    }

    PyObject **name_items = PySequence_Fast_ITEMS(names_seq);
    PyObject **count_items = PySequence_Fast_ITEMS(counts_seq);

    CollapseEntry *entries = PyMem_Malloc((size_t)n * sizeof(*entries));
    if (!entries) {
        Py_DECREF(names_seq);
        Py_DECREF(counts_seq);
        PyErr_NoMemory();
        return NULL;
    }

    for (Py_ssize_t i = 0; i < n; i++) {
        Py_ssize_t name_len;
        const char *s = PyUnicode_AsUTF8AndSize(name_items[i], &name_len);
        if (!s) {
            PyMem_Free(entries);
            Py_DECREF(names_seq);
            Py_DECREF(counts_seq);
            return NULL;
        }
        long c = PyLong_AsLong(count_items[i]);
        if (c == -1 && PyErr_Occurred()) {
            PyMem_Free(entries);
            Py_DECREF(names_seq);
            Py_DECREF(counts_seq);
            return NULL;
        }
        entries[i].name = s;
        entries[i].name_len = name_len;
        entries[i].count = c;
        entries[i].orig_idx = i;
    }

    qsort(entries, (size_t)n, sizeof(*entries), collapse_entry_cmp);

    char *consumed = PyMem_Calloc((size_t)n, sizeof(*consumed));
    PyObject *result = PyList_New(0);
    if (!consumed || !result) {
        PyMem_Free(entries);
        PyMem_Free(consumed);
        Py_XDECREF(result);
        Py_DECREF(names_seq);
        Py_DECREF(counts_seq);
        if (!PyErr_Occurred()) PyErr_NoMemory();
        return NULL;
    }

    for (Py_ssize_t i = 0; i < n; i++) {
        if (consumed[i]) continue;
        consumed[i] = 1;
        long cluster_count = entries[i].count;
        const char *anchor = entries[i].name;
        Py_ssize_t anchor_len = entries[i].name_len;

        PyObject *cluster_names = PyList_New(0);
        if (!cluster_names) goto error;
        if (PyList_Append(cluster_names, name_items[entries[i].orig_idx]) < 0) {
            Py_DECREF(cluster_names);
            goto error;
        }

        for (Py_ssize_t j = i + 1; j < n; j++) {
            if (consumed[j]) continue;
            if (entries[j].name_len != anchor_len) continue;
            const char *b = entries[j].name;
            int d = 0;
            int ok = 1;
            for (Py_ssize_t k = 0; k < anchor_len; k++) {
                if (anchor[k] != b[k]) {
                    if (++d > max_edit) { ok = 0; break; }
                }
            }
            if (!ok) continue;
            consumed[j] = 1;
            cluster_count += entries[j].count;
            if (PyList_Append(cluster_names, name_items[entries[j].orig_idx]) < 0) {
                Py_DECREF(cluster_names);
                goto error;
            }
        }

        PyObject *count_obj = PyLong_FromLong(cluster_count);
        if (!count_obj) {
            Py_DECREF(cluster_names);
            goto error;
        }
        PyObject *tuple = PyTuple_Pack(2, cluster_names, count_obj);
        Py_DECREF(cluster_names);
        Py_DECREF(count_obj);
        if (!tuple) goto error;
        if (PyList_Append(result, tuple) < 0) {
            Py_DECREF(tuple);
            goto error;
        }
        Py_DECREF(tuple);
    }

    PyMem_Free(entries);
    PyMem_Free(consumed);
    Py_DECREF(names_seq);
    Py_DECREF(counts_seq);
    return result;

error:
    PyMem_Free(entries);
    PyMem_Free(consumed);
    Py_DECREF(result);
    Py_DECREF(names_seq);
    Py_DECREF(counts_seq);
    return NULL;
}

static PyObject *shift_scores_impl(PyObject *self, PyObject *args) {
    PyObject *seq_obj;
    int k;

    if (!PyArg_ParseTuple(args, "Oi", &seq_obj, &k)) return NULL;
    if (k <= 0) {
        PyErr_SetString(PyExc_ValueError, "k must be positive");
        return NULL;
    }

    PyObject *seq_bytes = PyUnicode_AsASCIIString(seq_obj);
    if (!seq_bytes) return NULL;

    Py_ssize_t n = PyBytes_GET_SIZE(seq_bytes);
    if (n == 0) {
        Py_DECREF(seq_bytes);
        return PyList_New(0);
    }

    const char *seq = PyBytes_AS_STRING(seq_bytes);

    long *scores = PyMem_Malloc((size_t)n * sizeof(*scores));
    if (!scores) {
        Py_DECREF(seq_bytes);
        PyErr_NoMemory();
        return NULL;
    }

    // scores[i] = sum_{j=0..k-1} 4^(j+1) * value(seq[(i+j) mod n])
    // bases: a=0, c=1, t=2, g=3, anything else=0 (matches Python's dict.get(...,0))
    for (Py_ssize_t i = 0; i < n; i++) {
        long h = 0;
        long w = 4;
        Py_ssize_t p = i;
        for (int j = 0; j < k; j++) {
            int v;
            switch (seq[p]) {
                case 'a': case 'A': v = 0; break;
                case 'c': case 'C': v = 1; break;
                case 't': case 'T': v = 2; break;
                case 'g': case 'G': v = 3; break;
                default: v = 0; break;
            }
            h += w * v;
            w *= 4;
            if (++p == n) p = 0;
        }
        scores[i] = h;
    }

    // y for shift j is `scores` rotated by j → my, sum(y^2), syy are constant across j.
    // Compute via the same accumulation order as the Python reference so fp errors match.
    double mx = (double)(n + 1) / 2.0;
    double sxx = 0.0;
    for (Py_ssize_t i = 1; i <= n; i++) {
        double dx = (double)i - mx;
        sxx += dx * dx;
    }

    double sum_scores = 0.0;
    for (Py_ssize_t i = 0; i < n; i++) sum_scores += (double)scores[i];
    double my = sum_scores / (double)n;

    double syy = 0.0;
    for (Py_ssize_t i = 0; i < n; i++) {
        double dy = (double)scores[i] - my;
        syy += dy * dy;
    }

    PyObject *result = PyList_New(n);
    if (!result) {
        PyMem_Free(scores);
        Py_DECREF(seq_bytes);
        return NULL;
    }

    if (sxx == 0.0 || syy == 0.0) {
        for (Py_ssize_t j = 0; j < n; j++) {
            PyObject *value = PyFloat_FromDouble(NAN);
            if (!value) {
                Py_DECREF(result);
                PyMem_Free(scores);
                Py_DECREF(seq_bytes);
                return NULL;
            }
            PyList_SET_ITEM(result, j, value);
        }
        PyMem_Free(scores);
        Py_DECREF(seq_bytes);
        return result;
    }

    double denom = sqrt(sxx * syy);

    for (Py_ssize_t j = 0; j < n; j++) {
        double sxy = 0.0;
        Py_ssize_t p = j;
        for (Py_ssize_t i = 1; i <= n; i++) {
            double dx = (double)i - mx;
            double dy = (double)scores[p] - my;
            sxy += dx * dy;
            if (++p == n) p = 0;
        }
        PyObject *value = PyFloat_FromDouble(sxy / denom);
        if (!value) {
            Py_DECREF(result);
            PyMem_Free(scores);
            Py_DECREF(seq_bytes);
            return NULL;
        }
        PyList_SET_ITEM(result, j, value);
    }

    PyMem_Free(scores);
    Py_DECREF(seq_bytes);
    return result;
}

// primary function to be passed up to Python
static PyObject *window_compare_scores(PyObject *self, PyObject *args) {
    PyObject *sequence_obj;
    int kmer;
    PyObject *ws_obj;
    PyObject *we_obj;
    PyObject *wsc_obj;
    PyObject *wec_obj;

    if (!PyArg_ParseTuple(
            args,
            "OiOOOO",
            &sequence_obj,
            &kmer,
            &ws_obj,
            &we_obj,
            &wsc_obj,
            &wec_obj
        )) {
        return NULL;
    }

    if (kmer <= 0) {
        PyErr_SetString(PyExc_ValueError, "kmer must be positive");
        return NULL;
    }

    PyObject *sequence_bytes = PyUnicode_AsASCIIString(sequence_obj);
    if (!sequence_bytes) return NULL;

    int32_t *ws = NULL;
    int32_t *we = NULL;
    int32_t *wsc = NULL;
    int32_t *wec = NULL;
    Py_ssize_t n_windows = 0;
    Py_ssize_t n_window_ends = 0;
    Py_ssize_t n_window_starts_compare = 0;
    Py_ssize_t n_window_ends_compare = 0;

    if (!parse_int_sequence(ws_obj, &ws, &n_windows, "window_starts") ||
        !parse_int_sequence(we_obj, &we, &n_window_ends, "window_ends") ||
        !parse_int_sequence(wsc_obj, &wsc, &n_window_starts_compare, "window_starts_compare") ||
        !parse_int_sequence(wec_obj, &wec, &n_window_ends_compare, "window_ends_compare")) {
        Py_DECREF(sequence_bytes);
        PyMem_Free(ws);
        PyMem_Free(we);
        PyMem_Free(wsc);
        PyMem_Free(wec);
        return NULL;
    }

    if (n_windows != n_window_ends ||
        n_windows != n_window_starts_compare ||
        n_windows != n_window_ends_compare) {
        Py_DECREF(sequence_bytes);
        PyMem_Free(ws);
        PyMem_Free(we);
        PyMem_Free(wsc);
        PyMem_Free(wec);
        PyErr_SetString(PyExc_ValueError, "window arrays must have the same length");
        return NULL;
    }

    if (n_windows == 0) {
        Py_DECREF(sequence_bytes);
        PyMem_Free(ws);
        PyMem_Free(we);
        PyMem_Free(wsc);
        PyMem_Free(wec);
        return PyList_New(0);
    }

    Py_ssize_t sequence_len = PyBytes_GET_SIZE(sequence_bytes);
    Py_ssize_t n_kmers = sequence_len - (Py_ssize_t)kmer;
    if (n_kmers < 0) n_kmers = 0;
    if (n_kmers > INT32_MAX) {
        Py_DECREF(sequence_bytes);
        PyMem_Free(ws);
        PyMem_Free(we);
        PyMem_Free(wsc);
        PyMem_Free(wec);
        PyErr_SetString(PyExc_OverflowError, "sequence is too large for the stage06 accelerator");
        return NULL;
    }

    const char *sequence = PyBytes_AS_STRING(sequence_bytes);
    uint32_t *ids = NULL;
    uint32_t unique_ids = 0;
    if (n_kmers > 0 && !build_kmer_ids(sequence, (Py_ssize_t)kmer, n_kmers, &ids, &unique_ids)) {
        Py_DECREF(sequence_bytes);
        PyMem_Free(ws);
        PyMem_Free(we);
        PyMem_Free(wsc);
        PyMem_Free(wec);
        PyErr_NoMemory();
        return NULL;
    }

    uint32_t step = 0;
    PyObject *scores;
    if (can_use_sliding_suffix(n_windows, (int32_t)n_kmers, ws, we, wsc, wec, &step)) {
        scores = run_sliding_suffix_ids(
            ids, unique_ids, (int32_t)n_kmers, n_windows, step, ws, we, wsc, wec
        );
    } else {
        scores = run_generic_ids(
            ids, unique_ids, (int32_t)n_kmers, n_windows, ws, we, wsc, wec
        );
    }

    Py_DECREF(sequence_bytes);
    PyMem_Free(ids);
    PyMem_Free(ws);
    PyMem_Free(we);
    PyMem_Free(wsc);
    PyMem_Free(wec);
    return scores;
}

static PyMethodDef module_methods[] = {
    {
        "window_compare_scores",
        window_compare_scores,
        METH_VARARGS,
        PyDoc_STR("Compute stage06 window-comparison scores for one sequence.")
    },
    {
        "seq_win_score_int",
        seq_win_score_int_impl,
        METH_VARARGS,
        PyDoc_STR("Score a window by the proportion of singleton kmers.")
    },
    {
        "collapse_kmers",
        collapse_kmers_impl,
        METH_VARARGS,
        PyDoc_STR("Greedy Hamming-distance clustering of kmer names.")
    },
    {
        "shift_scores",
        shift_scores_impl,
        METH_VARARGS,
        PyDoc_STR("Pearson correlations of position-weighted kmer hashes against rotations.")
    },
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "_ext",
    "Native accelerator for trash2.stage06.",
    -1,
    module_methods,
};

PyMODINIT_FUNC PyInit__ext(void) {
    return PyModule_Create(&module_def);
}
