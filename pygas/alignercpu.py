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
from pygas.aligner import Aligner
from typing import List
from dataclasses import dataclass
from pygas.matrix import map_queries
from pygas.classes import AlignmentBatch


@dataclass
class AlignerCpu(Aligner):
    def align_queries(self, queries: List[str], keep_matrix=True) -> AlignmentBatch:
        self.alignment_batch = map_queries(
            targets=self.targets,
            penalty_max=self.max_penalty,
            hard_min=self.score_min,
            queries=queries,
            keep_matrix=keep_matrix,
            do_revcomp=self.rev_comp,
            exact_only=self.exact_only,
            match_type=self.match_type,
        )
        return self.alignment_batch
