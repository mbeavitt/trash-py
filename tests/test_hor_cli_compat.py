"""The `trash-py hor` subcommand accepts the upstream HORT.R getopt surface, so
replacing `Rscript HORT.R` with `trash-py hor` in an existing script Just Works."""
from __future__ import annotations

from pathlib import Path

from trash_py.hor_cli import build_hor_parser


def test_hortr_style_flags_parse():
    ns = build_hor_parser().parse_args([
        "-r", "reps.csv", "-o", "out", "-t", "25", "-l", "3",
        "-c", "178_1", "-m", "1", "-A", "Chr1", "-g", "Ath",
    ])
    assert ns.repeats == Path("reps.csv")
    assert ns.output == Path("out")
    assert ns.hor_threshold == 25 and ns.hor_min_len == 3
    assert ns.hor_class == "178_1" and ns.method == 1
    assert ns.chrA == "Chr1" and ns.genomeA == "Ath"


def test_underscore_long_names_parse():
    ns = build_hor_parser().parse_args([
        "--repeats", "r.csv", "--output_folder", "o", "--hor_threshold", "10",
        "--hor_min_len", "5", "--chrA", "Chr1", "--chrB", "Chr2",
        "--saveR", "n", "--plot_simple", "y",
    ])
    assert ns.hor_threshold == 10 and ns.hor_min_len == 5
    assert ns.chrA == "Chr1" and ns.chrB == "Chr2"
    assert ns.saveR == "n" and ns.plot_simple == "y"


def test_positional_repeats_and_chr_list():
    ns = build_hor_parser().parse_args(["reps.csv", "--chr-list", "Chr1,Chr2"])
    assert ns.repeats_pos == Path("reps.csv")
    assert ns.repeats is None
    assert ns.chr_list == "Chr1,Chr2"
