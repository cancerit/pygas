#
# Copyright (C) 2021
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
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List
from pygas.classes import AlignmentBatch

DNA_TRANS = ("ACGT", "TGCA")


@dataclass
class Aligner(ABC):
    """
    match_type:
        0 = true exact
        1 = query in target
        2 = target in query
        3 = all (1+2)
    """

    targets: List[str]
    rules: List[str]
    score_min: int
    rev_comp: bool = True
    match_type: int = 3

    def __post_init__(self):
        self._process_rules()

    @property
    def exact_only(self) -> bool:
        return self._exact_only

    @property
    def max_penalty(self) -> int:
        return self._max_penalty

    @property
    def min_penalty(self) -> int:
        return self._min_penalty

    def _process_rules(self):
        """
        Rule are passed in as a list of strings, they are exclusive rules to indicate the level of fuzzy matching
        permitted.  Each type of error is denoted as follows:

        * I = Insertion of 1 base
        * D = Deletion of 1 base
        * M = Mismatch of 1 base

        e.g. IMM = Allow an insertion and 2 mismatches, equivalent to a minimum mapping score of:

            query_len - (2+1+1)

        The higest value from the set of restrictions is used to determine the minimum score of a mapping to progress
        beyond the alignment matrix (backtrack to produce alignment used to confirm above).

        The rules are also used when determining which alignments of sufficient score are valid.
        """
        max_p = 0
        min_p = 1000000
        for rule in self.rules:
            uc_rule = rule.upper()
            d = uc_rule.count("D")
            i = uc_rule.count("I")
            m = uc_rule.count("M")
            # D/I are eqivalent to -2, mismatch -1 as a match scores 1
            # but as penalty reverse sign
            val = (d * 2) + (i * 2) + m
            if val > max_p:
                max_p = val
            if val < min_p:
                min_p = val
        if len(self.rules) == 0:
            min_p = 0
            self._exact_only = True
        else:
            self._exact_only = False
        self._max_penalty = max_p
        self._min_penalty = min_p

    @abstractmethod
    def align_queries(self, queries: List[str], keep_matrix=True) -> AlignmentBatch:  # pragma: no cover
        pass
