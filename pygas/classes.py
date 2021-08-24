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
from dataclasses import dataclass
from typing import List
from pygas.constants import GAP

BACKTRACK_STR_FMT = "Score: {score}, Cigar: {cigar}, MD: {md}, TargetId: {tid}, TargetPos: {tpos}\nEvents (D/I/M): {d}/{i}/{m}\nT: {tseq}\nM: {match}\nQ: {qseq}"

SPACE = " "
DASH = "-"
MD_S = "S"
MD_I = "I"
MD_D = "D"
MD_M = "M"


@dataclass
class ScoreMatrix:
    query: str
    target: str
    target_id: int
    score: int
    reversed: bool
    original_seq: str
    matrix: List[List[int]] = None
    exact: bool = False


@dataclass
class Backtrack:
    sm: ScoreMatrix
    match_mode: int = 3

    def __post_init__(self):
        self.events = {
            "D": 0,
            "I": 0,
            "M": 0,
        }
        if self.sm.matrix:
            self.backtrack()
            self.cigar_md()
        elif self.sm.exact:
            seq_len = len(self.sm.query)
            self.md = str(seq_len)
            self.cigar = f"{seq_len}M"
            self.nm = 0
            self.align_target = self.sm.target
            self.align_match = "|" * seq_len
            self.align_query = self.sm.query
            self.t_pos = 1
            self.pass_mode = True
        else:
            assert False, "Backtrack object must be created with matrix or be an exact mapping."

    def __str__(self):
        return BACKTRACK_STR_FMT.format(
            score=self.sm.score,
            cigar=self.cigar,
            md=self.md,
            tid=self.sm.target_id,
            tpos=self.t_pos,
            tseq=self.align_target,
            match=self.align_match,
            qseq=self.align_query,
            d=self.events["D"],
            i=self.events["I"],
            m=self.events["M"],
        )

    def print_matrix(self):
        for r in self.sm.matrix:
            for c in r:
                print(" {:02d} ".format(c), end="")
            print()

    def matrix_limits(self):
        i = None
        j = None
        score = self.sm.score
        f = self.sm.matrix
        for i_t, r in enumerate(f):
            for j_t, c in enumerate(r):
                if c == score:
                    i = i_t
                    j = j_t
                    break
            if i is not None:
                break
        return (i, j)

    def pass_rules(self, rules: List[str]):
        # if no rules, i.e. exact within mapping
        if not rules:
            if self.events["D"] == 0 and self.events["I"] == 0 and self.events["M"] == 0:
                return True
        else:
            for r in rules:
                if (
                    self.events["D"] <= r.count("D")
                    and self.events["I"] <= r.count("I")
                    and self.events["M"] <= r.count("M")
                ):
                    return True
        return False

    def backtrack(self):
        # find the max score location, only finding first instance
        target = self.sm.target
        query = self.sm.query
        (i, j) = self.matrix_limits()
        f = self.sm.matrix

        t_align = target[i + 1 :]  # i/X
        m_align = ""
        q_align = query[j + 1 :]  # j/Y
        while i >= 0 or j >= 0:
            m_chr = " "
            if i < 0:
                m_chr = " " * (j + 1)
                t_align = m_chr + t_align
                q_align = query[0 : j + 1] + q_align
                j = -1
            elif j < 0:
                m_chr = " " * (i + 1)
                t_align = target[0 : i + 1] + t_align
                q_align = m_chr + q_align
                i = -1
            elif ((i == 0 and j != 0) or (i != 0 and j == 0)) and target[i] == query[j]:
                m_chr = "|"
                t_align = target[i] + t_align
                q_align = query[j] + q_align
                i = i - 1
                j = j - 1
            elif i == 0 and j == 0:
                if target[i] == query[j]:
                    m_chr = "|"
                t_align = target[i] + t_align
                q_align = query[j] + q_align
                i = i - 1
                j = j - 1
            elif i >= 0 and j >= 0 and f[i][j] == f[i - 1][j - 1] + int(target[i] == query[j]):
                t_align = target[i] + t_align
                q_align = query[j] + q_align
                if target[i] == query[j]:
                    m_chr = "|"
                i = i - 1
                j = j - 1
            elif i > 0 and f[i][j] == f[i - 1][j] + GAP:
                t_align = target[i] + t_align
                q_align = "-" + q_align
                i = i - 1
            else:
                t_align = "-" + t_align
                q_align = query[j] + q_align
                j = j - 1
            m_align = m_chr + m_align

        # Padding is to simplify cigar/md, it stops as soon as the query is exhausted
        pad = len(q_align)
        self.align_target = t_align.ljust(pad, SPACE)
        self.align_match = m_align.ljust(pad, SPACE)
        self.align_query = q_align.ljust(pad, SPACE)

    def cigar_md(self):
        """
        CIGAR: See README for link to doc
        MD: See README for link to doc
        """
        # use align_* to create the cigar string
        ops = []
        md_ops = []
        nm = 0

        t_seq = self.align_target
        m_chr = self.align_match
        q_seq = self.align_query
        q_len = len(q_seq)
        # calculate start in target for pos, 1-based
        t_pos = (len(q_seq) - len(q_seq.lstrip())) + 1

        for i in range(q_len):
            (t, m, q) = (t_seq[i], m_chr[i], q_seq[i])
            if t == SPACE and q not in [SPACE, DASH]:
                ops.append(MD_S)
                continue
            if t == DASH:
                self.events["I"] += 1
                ops.append(MD_I)
                nm += 1
                continue
            if q == DASH:
                self.events["D"] += 1
                ops.append(MD_D)
                md_ops.append(t.lower())
                nm += 1
                continue
            if t == q:
                # match
                ops.append(MD_M)
                md_ops.append(MD_M)
                continue
            if q == SPACE:
                continue
            if m == SPACE:
                # if all mismatch to end are SPACE we can soft clip and break
                if m_chr[i:].isspace():
                    ops.extend([MD_S] * (q_len - i))
                    break
                if t not in [SPACE, DASH] and q not in [SPACE, DASH] and (len(ops) == 0 or ops[-1] == MD_S):
                    ops.append(MD_S)
                    t_pos += 1
                    continue
                self.events["M"] += 1
                ops.append(MD_M)
                md_ops.append(t)
                nm += 1
                continue

        (cigar, cigar_ops, cigar_len) = self._ops_to_cigar(ops)
        md = self._opt_to_md(md_ops)
        self.t_pos = t_pos
        self.cigar = cigar
        self.md = md
        self.nm = nm
        self.pass_mode = self._match_mode_validation(cigar_ops, cigar_len)

    def _match_mode_validation(self, cigar_ops, cigar_len):
        if self.match_mode == 3:
            # anything
            return True
        if self.match_mode == 0:
            # these need merging
            if (
                self.t_pos != 1
                or self.align_target.startswith(" ")
                or "S" in cigar_ops
                or len(self.align_query) != len(self.align_target.strip())
            ):
                return False
            return True
        if self.match_mode == 1:
            # QinT
            if self.align_target.startswith(" ") or len(self.align_query) > len(self.align_target.rstrip()):
                return False
            return True
        # only match_mode=2 left
        if (
            self.align_query.startswith(" ")
            or self.align_query.endswith(" ")
            or len(self.align_query) != len(self.align_target)
        ):
            return False
        return True

    def _ops_to_cigar(self, ops: List[str]) -> str:
        cigar_ops = [ops[0]]
        cigar_len = [1]
        for i in range(1, len(ops)):
            if ops[i] == cigar_ops[-1]:
                cigar_len[-1] += 1
            else:
                cigar_ops.append(ops[i])
                cigar_len.append(1)
        cigar = ""
        for i, type in enumerate(cigar_ops):
            cigar += f"{cigar_len[i]}{type}"
        return (cigar, cigar_ops, cigar_len)

    def _opt_to_md(self, md_ops: List[str]) -> str:
        md_clean = [md_ops[0]]
        for i in range(1, len(md_ops)):
            op = md_ops[i]
            # match, numeric
            if op == "M":
                if md_clean[-1].startswith("M"):
                    md_clean[-1] += op
                else:
                    md_clean.append(op)
                continue
            # deletions
            if op.islower():
                if md_clean[-1].islower():
                    md_clean.append(op)
                else:
                    md_clean.extend(["^", op])
                continue
            # mismatch
            md_clean.append(op)

        md_str = ""
        for i, mde in enumerate(md_clean):
            if mde.startswith("M"):
                md_str += str(len(mde))
                continue
            if mde == "^":
                md_str += mde
                continue
            if mde.isupper():
                if i > 0 and md_clean[i - 1].islower():
                    md_str += "0"
                md_str += mde
            else:
                # inverse of preceding block is not relevant as rules of scoring will never allow del/ins to precede a
                # mismatch as when tracing back the mismatch has a lower penalty, always get "--m"
                md_str += mde.upper()
        return md_str


@dataclass
class AlignmentBatch:
    unmapped: List[str]
    mapped: List[List[Backtrack]]

    def __post_init__(self):
        self.total_reads = len(self.mapped) + len(self.unmapped)

    def mapped_fraction(self):
        return len(self.mapped) / self.total_reads

    def unmapped_fraction(self):
        return len(self.unmapped) / self.total_reads
