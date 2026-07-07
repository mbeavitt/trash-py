# HOR module — bugs & quirks in the upstream source

While porting the TRASH HOR module (`HORT.R` + the `HOR.V3.3` binary) to
`trash-py`, I reverse-engineered the binary (its shipped source `main.c` is only
a partial draft — see #1) and diffed the port against the reference `hor_output`
golden. This is the catalogue of upstream bugs and quirks found along the way.

**The Python port reproduces every one of these deliberately**, so its HOR
tables are byte-for-byte identical to the original tool (the one exception being
#6, a formatting choice we intentionally change — full integers instead of R's
scientific contractions, per project preference).

Legend: 🐞 = genuine bug (wrong result), ⚠️ = off-by-one, 📎 = quirk / not a bug.

---

## 1. 📎 `main.c` is not the real source and does not compile

The `main.c` the author shared is a `V3.0` draft; the shipped binary is `V3.3`.
Two concrete divergences:

* **`isSimilar` is declared inside the `while` loop but used after it** (the
  "edge-close" `if (isSimilar && openHOR >= cutoff)` at line 175). This does not
  compile under any modern gcc (`'isSimilar' undeclared`). The real binary
  declares it at function scope, so the edge-close sees the *last* value from the
  loop. The port matches the binary.
* **`main.c` keeps newlines in the sequence buffer; the binary strips them**
  (see #5). Another silent behavioural difference.

Consequence: `main.c` cannot be trusted as a spec. The port was derived from the
disassembly of `HOR.V3.3` and validated against the golden output.

## 2. ⚠️🐞 Edge-close under-counts the last pair's SNVs (`total_variant`)

The SNV accumulator `snvBack` is captured at the **top** of each iteration,
*before* the current pair is compared. When a HOR run ends because a pair fails
(mid-loop close), that's correct — `snvBack` holds the sum over the whole run.
But when a run runs off the end of a diagonal (**edge close**), the code emits
`snvBack` from the final iteration, which was captured *before* the last matched
pair was added. So **edge-terminated HORs report `total_variant` short by the
last pair's SNV count.**

Demonstration (three matching pairs, 1 SNV each):

| case | how the run ends | reported `total_variant` | correct |
|------|------------------|--------------------------|---------|
| run reaches array boundary | edge close | **2** | 3 |
| a 7th repeat breaks the run first | mid close | 3 | 3 |

Covered by `tests/test_hor_core.py::test_edge_close_undercounts_last_pair` and
`::test_mid_close_counts_all_pairs`.

## 3. 🐞 SNVs carry over between HORs (inflated `total_variant`)

`compareAB` adds a pair's SNV count to the running total **whenever
`snv <= threshold`**, independent of whether the pair is actually "similar"
(similarity also requires the strand/direction test). The accumulator only
resets when an *open* HOR closes. Therefore:

* SNVs from pairs that pass the threshold but fail the direction test, and
* SNVs from any threshold-passing pairs seen *before* a HOR opens on that
  diagonal,

both **leak into the next HOR's `total_variant`**. This is why `total_variant`
routinely exceeds `threshold × block_size` — e.g. values up to ~732 appear with
`threshold_SNV = 44`. It also makes `total_variant` (and hence `SNV_per_kbp`)
depend on where the diagonal walk started, not just on the HOR itself.

## 4. 🐞 `SNV_per_kbp` divides by block A twice (`HORT.R`)

```r
hors$SNV_per_kbp <- 1000 * hors$total_variant / ((hors$block.A.size.bp + hors$block.A.size.bp) / 2)
```

The second `block.A.size.bp` should be `block.B.size.bp`. As written the average
is just `block.A.size.bp`, so block B's length never enters the normalisation:
effectively `SNV_per_kbp = 1000 * total_variant / block.A.size.bp`.

## 5. 📎 Wrapped-line handling (clarification, not a bug in the binary)

Aligned FASTAs are written wrapped at 60 columns (seqinr default). `main.c`
would store the embedded newlines and compare over `alilength` *bytes*, silently
dropping the last few real columns. The real binary **skips newlines** while
reading, comparing over real alignment columns. The port strips newlines to
match the binary; wrapping width is therefore irrelevant to results.
(`tests/test_hor_core.py::test_wrapped_alignment_strips_newlines`.)

## 6. 📎 Numeric columns are R doubles → scientific contraction in the CSV

`block.size.in.units`, `block.A.size.bp`, `block.B.size.bp` are computed with
`... + 1`; the literal `1` is a double in R, so the whole column is `double`.
On write, R (`scipen = 0`) contracts round numbers to the shorter form, e.g.
`200000` → `2e+05`, but leaves `199828` alone — so a single column mixes
`"535"` and `"2e+05"`. **The port intentionally writes the full integer
(`200000`)** rather than reproduce this contraction (project preference; these
are still numerically identical).

## 7. ⚠️🐞 `horb.R` region-B overlap count is off by one

```r
B_repeats_found_in_A = sum(repeats$hors_formed_count[split_after : nrow(repeats)] > 0)
```

Region B is repeats `(split_after + 1) : nrow`, but the index starts at
`split_after`, so the **last region-A repeat is counted as part of region B**.
(The port reproduces the upstream index so the `summary_of_hors` numbers match.)

## 8. 🐞 A lone HOR is silently discarded (`nrow(hors) > 1`)

Both `hort.R` and `horb.R` only process results when `nrow(hors) > 1`. If exactly
**one** HOR is found it is dropped and the run reports "no HORs" (and in
`hort.R`, via the typo'd message `"Noh HORs identified"`). The port keeps this
`<= 1 ⇒ skip` behaviour for parity.

## 9. 📎 Method-1 repetitiveness counts only the A block

In method 1, `hors_formed_count` is incremented for the A-block repeats only
(`start_A : end_A`), never the B block, whereas method 2 counts both. So a repeat
that appears only as the *B* side of a HOR is not credited in the method-1
`repeats_with_hors` table. Faithfully reproduced.

## 10. 📎 Dead `p <- p - 1` and single-record infinite loop

* `filter_different_blocks` (default `FALSE`) contains `p <- p - 1` inside a
  `for (p in rev(...))` loop — a no-op, since `p` is overwritten each iteration.
* `main.c`'s length-counting loop `for (c = fgetc(...); c != '>'; ...)` has no
  `EOF` guard, so a single-sequence alignment loops forever. Not reachable in
  practice (HOR inputs always have ≥ 2 repeats).

---

### Validation

* `_ext.find_hors` matches the `HOR.V3.3` binary on **3000 randomised fuzz
  cases** (methods 1 and 2, varied lengths, gaps, thresholds, cutoffs, wrapping).
* The full pipeline reproduces the upstream `hor_output` golden **byte-for-byte**
  for all five *A. thaliana* chromosomes (up to 1.45M HOR rows each), modulo the
  #6 scientific-notation formatting we deliberately keep as full integers, plus
  ~5 rows in ~3.9M where R's and Python's float formatters round a last-digit tie
  in `SNV_per_kbp` in opposite directions (e.g. `…049` vs `…0491`). This is the
  same `%.15g` handling the rest of trash-py uses for R-parity floats.
