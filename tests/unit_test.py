#
# Copyright (c) 2021
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of pygas.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
import pytest

from pygas.aligner import Aligner
from pygas.alignercpu import AlignerCpu

# from pygas.alignergpu import AlignerGpu
from pygas.classes import Backtrack

READ_A = "ACGTAAAAAAAAAAAACGT"
READ_C = "ACGTCCCCCCCCCCCCCGT"
READ_G = "ACGTGGGGGGGGGGGACGT"
READ_T = "ACGTTTTTTTTTTTTACGT"
READ_A_MM = "ACGTAAAAATAAAAAACGT"  # mismatch
READ_A_DMM = "ACGTAAAATTAAAAAACGT"  # double mismatch
READ_BAD = "AAAAAAAAAAAAAAAAAAA"
READ_D = "ACGTAAAAAAAAAAACGTT"
READ_I = "ACGTAAAAAAAAAAAAACG"
READ_A_REV = READ_T  # crafted like this for a reason
TARGETS = [READ_A, READ_C, READ_G, READ_T]
MIN_SCORE = 15  # likely need a test to validate input reads aren't fully excluded by penalty vs read length
RULES_NONE = []  # exact match required to get a hit
RULES_MM = ["M"]
RULES_I = ["I"]
RULES_D = ["D"]

MATRIX = " 01  01  00  00 \n 01  02  01  01 \n 01  02  02  02 \n 01  02  02  03 \n"
BT_STR = "Score: 3, Cigar: 4M, MD: 2A1, TargetId: 0, TargetPos: 1\nEvents (D/I/M): 0/0/1\nT: AAAA\nM: || |\nQ: AATA\n"


def test_01_test_instantiate_aligner():
    with pytest.raises(TypeError) as e_info:
        a = Aligner(targets=TARGETS, rules=[], score_min=MIN_SCORE)


def test_02_instantiate_alignerCpu():
    a = AlignerCpu(targets=TARGETS, rules=[], score_min=MIN_SCORE, rev_comp=True)
    assert isinstance(a, Aligner), "AlignerCpu is instance of Aligner"
    assert len(a.targets) == 4, "Input target set == 4 seqs"


@pytest.mark.parametrize(
    "test_input, min_expected, max_expected",
    [
        ([], 0, 0),
        (["MM"], 2, 2),
        (["IM"], 3, 3),
        (["DM"], 3, 3),
        (["IDM"], 5, 5),
        (["MM", "IDM"], 2, 5),
        (["MM", "M"], 1, 2),
    ],
)
def test_03_aligner_penalties(test_input, min_expected, max_expected):
    a = AlignerCpu(
        targets=TARGETS,
        rules=test_input,
        score_min=MIN_SCORE,
        rev_comp=True,
    )
    assert a.min_penalty == min_expected
    assert a.max_penalty == max_expected


@pytest.mark.parametrize(
    "test_reads, expected_tid, expected_count",
    [
        (["AAAAAAAC"], 0, 1),
        (["CCCCCCCA"], 1, 1),
        (["GGGGGGGT"], 2, 1),
        (["TTTTTTTG"], 3, 1),
    ],
)
def test_04_aligner_tids(test_reads, expected_tid, expected_count):
    a = AlignerCpu(
        targets=["AAAAAAAC", "CCCCCCCA", "GGGGGGGT", "TTTTTTTG"],
        rules=[],
        score_min=8,
        rev_comp=True,
    )
    ab = a.align_queries(test_reads)
    assert len(ab.mapped[0]) == expected_count
    assert ab.mapped[0][0].sm.target_id == expected_tid


@pytest.mark.parametrize(
    "test_reads, expected_count",
    [
        (["AAAAAAAA"], 1),
        (["CCCCCCCC"], 1),
        (["GGGGGGGG"], 1),
        (["TTTTTTTT"], 1),
    ],
)
def test_05_aligner_revcomp(test_reads, expected_count):
    a = AlignerCpu(
        targets=["AAAAAAAA", "CCCCCCCC"],
        rules=[],
        score_min=8,
        rev_comp=True,
    )
    ab = a.align_queries(test_reads)
    assert len(ab.mapped[0]) == expected_count


@pytest.mark.parametrize(
    "test_reads, expected_count, pos, cigar, md",
    [
        (["TTTTTTTT"], 1, 1, "8M", "8"),
    ],
)
def test_06_aligner_nomatch(test_reads, expected_count, pos, cigar, md):
    a = AlignerCpu(targets=["AAAAAAAA"], rules=[], score_min=8, rev_comp=True)
    ab = a.align_queries(test_reads)
    assert len(ab.mapped) == expected_count
    bt = ab.mapped[0][0]
    assert bt.t_pos == pos
    assert bt.cigar == cigar
    assert bt.md == md


@pytest.mark.parametrize(
    "test_input, rules, expected",
    [
        ([READ_A], RULES_NONE, 1.0),
        ([READ_C], RULES_NONE, 1.0),
        ([READ_G], RULES_NONE, 1.0),
        ([READ_T], RULES_NONE, 1.0),
        ([READ_A, READ_C, READ_G, READ_T, READ_BAD], RULES_NONE, 0.8),
        ([READ_A, READ_A_MM], RULES_NONE, 0.5),
        # these do not distinguish between type of event in this test
        ([READ_A, READ_A_MM], RULES_MM, 1),
        ([READ_A, READ_A_MM], RULES_I, 1),  # insert will allow 2xMM
        ([READ_A, READ_A_MM], RULES_D, 1),  # del will allow 2xMM
        ([READ_D], RULES_NONE, 0.0),  # del not allowed if no rule
        ([READ_D], RULES_MM, 0.0),  # del not allowed when single MM
        ([READ_D], RULES_D, 1.0),  # del OK
        ([READ_A_DMM], RULES_D, 1.0),  # del allows double MM
        ([READ_I], RULES_NONE, 0.0),  # ins not allowed if no rule
        ([READ_I], RULES_MM, 0.0),  # ins not allowed when single MM
        ([READ_I], RULES_I, 1.0),  # ins OK
        ([READ_A_DMM], RULES_I, 1.0),  # ins allows double MM
    ],
)
def test_10_align_reads_cpu(test_input, rules, expected):
    a = AlignerCpu(
        targets=TARGETS,
        rules=rules,
        score_min=MIN_SCORE,
        rev_comp=True,
    )
    ab = a.align_queries(test_input)
    assert ab.mapped_fraction() == expected, "Aligned queries"


@pytest.mark.parametrize(
    "match_type, num_expected",
    [
        (3, 2),
        (0, 1),
        (1, 2),
        (2, 1),
    ],
)
def test_11_multialign_cpu(match_type, num_expected):
    a = AlignerCpu(
        targets=["AAAAAAAA", "TAAAAAAAAT"],
        rules=["MM"],
        score_min=8,
        rev_comp=True,
        match_type=match_type,
    )
    ab = a.align_queries(["AAAAAAAA"])
    # reads mapped is still 1
    assert len(ab.mapped) == 1
    # but number of mappings for the reads is 2
    assert len(ab.mapped[0]) == num_expected


@pytest.mark.parametrize(
    "targets, reads, exp_t, exp_m, exp_q",
    [
        (["TAAAAT"], ["TAAAAT"], "TAAAAT", "||||||", "TAAAAT"),  # full match
        (["TAAAAT"], ["TATAAT"], "TAAAAT", "|| |||", "TATAAT"),  # single mismatch
        (
            ["ACCATTACCATTACC"],
            ["ACCATTACCATACC"],
            "ACCATTACCATTACC",
            "|||||||||| ||||",
            "ACCATTACCA-TACC",
        ),  # deletion
        (
            ["ACCATTACCATTACC"],
            ["ACCATTACCCATTACC"],
            "ACCATTA-CCATTACC",
            "||||||| ||||||||",
            "ACCATTACCCATTACC",
        ),  # insert
        (["AAAACCCC"], ["CCCC"], "AAAACCCC", "    ||||", "    CCCC"),  # late start
        (
            ["AAAACCCC"],
            ["CCCCAA"],
            "AAAACCCC  ",
            "    ||||  ",
            "    CCCCAA",
        ),  # late start, hangover
        (["CCCC"], ["AACCCC"], "  CCCC", "  ||||", "AACCCC"),  # early start, hangover
        (["CCCCCC"], ["CCCC"], "CCCCCC", "||||", "CCCC"),  # early start, early stop
        (["AGAGGG"], ["GGAGGG"], "AGAGGG", " |||||", "GGAGGG"),  # first base mismatch
    ],
)
def test_15_backtrack(targets, reads, exp_t, exp_m, exp_q):
    # easiest to build the thing we want and call bt
    a = AlignerCpu(targets=targets, rules=["MDI"], score_min=4, rev_comp=True)
    bt = a.align_queries(reads).mapped[0][0]
    assert bt.align_target == exp_t
    assert bt.align_match == exp_m
    assert bt.align_query == exp_q


def test_16_backtrack_no_matrix():
    a = AlignerCpu(targets=["AAAA"], rules=[], score_min=4, rev_comp=True)
    bt = a.align_queries(["AAAA"]).mapped[0][0]
    assert bt.cigar == "4M"
    assert bt.md == "4"
    assert bt.t_pos == 1
    assert bt.align_target == "AAAA"
    assert bt.align_match == "||||"
    assert bt.align_query == "AAAA"


def test_17_backtrack_no_matrix_not_exact():
    a = AlignerCpu(targets=["ACCA"], rules=["MM"], score_min=2, rev_comp=True)
    bt = a.align_queries(["AAAA"], keep_matrix=False).mapped[0][0]
    # slightly circular but confirms all bt operations were performed and the matrix was cleared
    with pytest.raises(AssertionError):
        Backtrack(bt.sm)


@pytest.mark.parametrize(
    "target, query, e_tpos, e_cigar, e_md",
    [
        (
            "GAGCATTCGGATTTCCCGA",
            "GAGCATTCGGATTTCCCGT",
            1,
            "18M1S",
            "18",
        ),  # a trailing soft clip for equal length
        # T:  AGCATTCGGATTTCCCGAA
        # M:  ||||||||||||||||||
        # Q: TAGCATTCGGATTTCCCGA
        ("AGCATTCGGATTTCCCGAA", "TAGCATTCGGATTTCCCGA", 1, "1S18M", "18"),
        # T: GAGCATTCGGATTTCCCGA
        # M:  ||||||||||||||||||
        # Q: TAGCATTCGGATTTCCCGA
        ("GAGCATTCGGATTTCCCGA", "TAGCATTCGGATTTCCCGA", 2, "1S18M", "18"),
        ("AAACCCTTTGGG", "AAACCCTTTGGG", 1, "12M", "12"),  # exact
        ("AAACCCGGGTTT", "AAACCCCGGTTT", 1, "12M", "6G5"),  # 1 mismatch
        ("ACCCTTTGGG", "AAACCCTTTGGG", 1, "2S10M", "10"),  # soft-clip start
        ("AAACCCTTTGGG", "ACCCTTTGGG", 3, "10M", "10"),  # not start of target
        ("AAACCCTTTG", "AAACCCTTTGGG", 1, "10M2S", "10"),  # soft-clip end
        (
            "AAACCCTTTGGG",
            "AAACCCTTGGG",
            1,
            "11M",
            "8T2",
        ),  # deletion, but mismatch is better score
        (
            "AAAAAAAAAACCCTTTCGCGCGCGCG",
            "AAAAAAAAAACCCTTCGCGCGCGCG",
            1,
            "13M1D12M",
            "13^T12",
        ),  # deletion of base
        ("AAACCCTTTTCACACA", "AAACCCTTTTTCACACA", 1, "6M1I10M", "16"),  # Insert 1
        (
            "AATTTATATATATAACGTCGCGCGCGAAA",
            "AATTTATATATATGGTCGCGCGCGAAA",
            1,
            "13M2D14M",
            "13^AA0C13",
        ),  # deletion precedes to mismatch
        (
            "CGCGCGCGTCGCGCGCG",
            "CGCGCGCGCCCGCGCGCG",
            1,
            "8M1I9M",
            "8T8",
        ),  # insertion precedes mismatch
        # T: CTTACTG-CGTCAACGGCTA
        # M: ||||||| ||||||||||
        # Q: CTTACTGCCGTCAACGGCN
        ("CTTACTGCGTCAACGGCTA", "CTTACTGCCGTCAACGGCN", 1, "7M1I10M1S", "17"),
        # rules of scoring will never allow del/ins to precede a mismatch as when tracing back the mismatch has a lower penalty
    ],
)
def test_18_backtrack_cigar(target, query, e_tpos, e_cigar, e_md):
    a = AlignerCpu(
        targets=[target],
        rules=["MDDDDDI"],
        score_min=1,
        rev_comp=True,
    )
    bt = a.align_queries([query], keep_matrix=True).mapped[0][0]
    assert bt.t_pos == e_tpos
    assert bt.cigar == e_cigar
    assert bt.md == e_md


def test_20_backtrack_out(capsys):
    a = AlignerCpu(targets=["AAAA"], rules=["M"], score_min=3, rev_comp=True)
    bt = a.align_queries(["AATA"], keep_matrix=True).mapped[0][0]
    bt.print_matrix()  # don't delete captured by test
    out, err = capsys.readouterr()
    assert out == MATRIX
    print(bt)  # don't delete captured by test
    out, err = capsys.readouterr()
    assert out == BT_STR


@pytest.mark.parametrize(
    "target, query, mode, exp_map, exp_unmap",
    [
        # exact match, equal length
        # T: AAAAAA
        # M: ||||||
        # Q: AAAAAA
        ("AAAAAA", "AAAAAA", 0, 1, 0),  # pass all modes
        ("AAAAAA", "AAAAAA", 1, 1, 0),  # pass all modes
        ("AAAAAA", "AAAAAA", 2, 1, 0),  # pass all modes
        ("AAAAAA", "AAAAAA", 3, 1, 0),  # pass all modes
        # target-3' overhang
        # T: AAAAAA
        # M: |||||
        # Q: AAAAA
        ("AAAAAA", "AAAAA", 0, 0, 1),  # fail as exact required
        ("AAAAAA", "AAAAA", 1, 1, 0),  # pass QinT
        ("AAAAAA", "AAAAA", 2, 0, 1),  # fail T>Q
        ("AAAAAA", "AAAAA", 3, 1, 0),  # any, pass same as QinT
        # query-3' overhang
        # T: AAAAA
        # M: |||||
        # Q: AAAAAA
        ("AAAAA", "AAAAAA", 0, 0, 1),  # fail as exact required
        ("AAAAA", "AAAAAA", 1, 0, 1),  # fail as Q>T
        ("AAAAA", "AAAAAA", 2, 1, 0),  # pass TinQ
        ("AAAAA", "AAAAAA", 3, 1, 0),  # any, pass same as TinQ
        # Offset due to 5' mismatch
        # T:  AAAAAA
        # M:  |||||
        # Q: GAAAAA
        ("AAAAAA", "GAAAAA", 0, 0, 1),  # fail as exact required
        ("AAAAAA", "GAAAAA", 1, 0, 1),  # fail as QinT required
        ("AAAAAA", "GAAAAA", 2, 0, 1),  # fail as TinQ required
        ("AAAAAA", "GAAAAA", 3, 1, 0),  # pass as "anything"
        # Trailing query soft-clip
        # T: AAAAAA
        # M: |||||
        # Q: AAAAAG
        ("AAAAAA", "AAAAAG", 0, 0, 1),  # fail as exact required
        ("AAAAAA", "AAAAAG", 1, 1, 0),  # pass as QinT required (within the bounds)
        ("AAAAAA", "AAAAAG", 2, 1, 0),  # pass as TinQ required (within the bounds)
        ("AAAAAA", "AAAAAG", 3, 1, 0),  # pass as "anything"
        # mismatch mapping, exact but forced start!=1
        # T: AAAACA
        # M:  |||||
        # Q:  AAACA
        ("AAAACA", "AAACA", 0, 0, 1),  # fail as requires start=1
        ("AAAACA", "AAACA", 1, 1, 0),  # TinQ should pass
        ("AAAACA", "AAACA", 2, 0, 1),  # fail as Q<T
        ("AAAACA", "AAACA", 3, 1, 0),  # pass if 0-2 pass
        # mapping with central mismatch
        # T: AAAAAAA
        # M: ||| |||
        # Q: AAACAAA
        ("AAAAAAA", "AAACAAA", 0, 1, 0),  # pass as full overlap
        ("AAAAAAA", "AAACAAA", 1, 1, 0),  # pass as full overlap
        ("AAAAAAA", "AAACAAA", 2, 1, 0),  # pass as full overlap
        ("AAAAAAA", "AAACAAA", 3, 1, 0),  # pass if 0-2 pass
    ],
)
def test_21_match_mode_fuzzy(target, query, mode, exp_map, exp_unmap):
    a = AlignerCpu(targets=[target], rules=["M"], rev_comp=False, match_type=mode, score_min=3)
    results = a.align_queries([query], keep_matrix=False)
    # print(results)
    assert len(results.mapped) == exp_map
    assert len(results.unmapped) == exp_unmap
