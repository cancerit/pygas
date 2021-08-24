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
from typing import List
import gzip
import sys

from pygas.alignercpu import AlignerCpu
from pygas.classes import Backtrack


def _simple_seq_load(ifile: str) -> List[str]:  # pragma: no cover
    try:
        with gzip.open(ifile, "rt") as i_fh:
            lines = i_fh.readlines()
    except gzip.BadGzipFile:
        with open(ifile, "rt") as i_fh:
            lines = i_fh.readlines()
    lines = [x.strip() for x in lines]
    return lines


def run(targets, queries, output, minscore, rules, allow_rev_comp):  # pragma: no cover
    target_seqs = _simple_seq_load(targets)
    query_seqs = _simple_seq_load(queries)
    a = AlignerCpu(
        targets=target_seqs,
        rules=rules,
        score_min=minscore,
        rev_comp=allow_rev_comp,
    )
    ab = a.align_queries(query_seqs)
    ofh = open(output, "w") if output else sys.stdout
    print(
        "\t".join(
            [
                "#query",
                "reversed",
                "t_id",
                "t_pos",
                "seq",
                "cigar",
                "md",
                "repeat_2-7...",
            ]
        ),
        file=ofh,
    )
    for unmapped_seq in ab.unmapped:
        print(f"{unmapped_seq}\t.\t.\t.\t.", file=ofh)
    for results in ab.mapped:
        # get max score
        max_score = max([bt.sm.score for bt in results])

        res_ele = [results[0].sm.original_seq]
        for bt in results:
            sm = bt.sm
            if sm.score == max_score:
                res_ele.extend(
                    [
                        str(sm.reversed),
                        str(sm.target_id),
                        str(bt.t_pos),
                        sm.query,
                        bt.cigar,
                        bt.md,
                    ]
                )
        print("\t".join(res_ele), file=ofh)

    if ofh is not sys.stdout:
        ofh.close()
