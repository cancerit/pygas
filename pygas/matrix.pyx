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
from typing import List, Tuple
from pygas.classes import ScoreMatrix, Backtrack, AlignmentBatch

cdef int GAP = -1
cdef int MATCH = 1
cdef int MISMATCH = 0

REVCOMP_TABLE = ''.maketrans({"A":"T", "C":"G", "G":"C", "T":"A"})

def revcomp(str seq) -> str:
    cdef str rev_seq = seq.upper()[::-1]
    return rev_seq.translate(REVCOMP_TABLE)

def map_queries(
    targets: List[str],
    queries: List[str],
    int penalty_max,
    int hard_min,
    keep_matrix=True,
    do_revcomp=True,
    exact_only=False,
    match_type: int=0,
) -> AlignmentBatch:
    """
    max_penalty is applied on a per-read basis
    if it is below hard_min then a warning is emitted and the item is not sent for alignment.
    """
    cdef str query, target, rev_query
    cdef int t_idx, q_len, t_len, min_score, max_score, max_t_len
    cdef dict targets_by_seq = {}
    cdef list target_lengths = []
    cdef list unmapped = []
    cdef list mapped = []

    max_t_len = -1

    do_substr = False
    for t_idx, target in enumerate(targets):
        target_lengths.append(len(target))
        if target in targets_by_seq:
            targets_by_seq[target].append(t_idx)
        else:
            targets_by_seq[target] = [t_idx]
            t_len = len(target)
            if t_len > max_t_len:
                if max_t_len != -1:
                    do_substr = True
                max_t_len = t_len

    for query in queries:
        q_len = len(query)
        if q_len < hard_min:
            unmapped.append(query)
            continue
        q_penalty_score = q_len - penalty_max
        if q_len != max_t_len:
            do_substr = True

        if do_revcomp:
            rev_query = revcomp(query)

        result = []  # store the best scoring results we see
        best_score = 0

        if do_substr and match_type != 0:  # try substr
            for t_idx, target in enumerate(targets):
                if ((match_type & 1) == 1 and query in target) or ((match_type & 2) == 2 and target in query):
                    t_len = target_lengths[t_idx]
                    # use matrix to build the bits we need quickly
                    min_score = max(hard_min, best_score, min(t_len - penalty_max, q_penalty_score))
                    (matrix, max_score) = fill(target, query, min_score)
                    if max_score > best_score:
                        result = []
                        best_score = max_score
                    if max_score < min_score:
                        continue
                    result.append(
                        ScoreMatrix(
                            query=query,
                            target=target,
                            target_id=t_idx,
                            score=max_score,
                            matrix=matrix,
                            reversed=False,
                            original_seq=query,
                        )
                    )
        elif query in targets_by_seq:
            best_score = q_len
            for t_idx in targets_by_seq[query]:
                result.append(
                    ScoreMatrix(
                        query=query,
                        target=query,  # they are equal
                        target_id=t_idx,
                        score=q_len,
                        reversed=False,
                        original_seq=query,
                        exact=True,
                    )
                )

        if do_revcomp:
            if do_substr:  # try substr
                for t_idx, target in enumerate(targets):
                    if target in rev_query or rev_query in target:
                        t_len = target_lengths[t_idx]
                        # use matrix to build the bits we need quickly
                        min_score = max(hard_min, best_score, min(t_len - penalty_max, q_penalty_score))
                        (matrix, max_score) = fill(target, rev_query, min_score)
                        if max_score > best_score:
                            result = []
                            best_score = max_score
                        if max_score < min_score:
                            continue
                        result.append(
                            ScoreMatrix(
                                query=rev_query,
                                target=target,
                                target_id=t_idx,
                                score=max_score,
                                matrix=matrix,
                                reversed=True,
                                original_seq=query,
                            )
                        )
            elif rev_query in targets_by_seq:
                best_score = q_len
                for t_idx in targets_by_seq[rev_query]:
                    result.append(
                        ScoreMatrix(
                            query=rev_query,
                            target=rev_query,  # they are equal
                            target_id=t_idx,
                            score=q_len,
                            reversed=True,
                            original_seq=query,
                            exact=True,
                        )
                    )

        if exact_only:
            if len(result) == 0:
                unmapped.append(query)
                continue
            mapped.append(result)
            continue

        # if we get to here and the query has been "exact" or "substr" mapped to a target of the maximum length no point in processing
        # the matrix, regardless of mismatch options
        if len(result) > 0:
            mapped.append(result)
            continue

        for t_idx, target in enumerate(targets):
            t_len = target_lengths[t_idx]
            if best_score >= t_len:
                # we only need to try for a better score if a target is longer, it'll be a multi-map though
                continue
            if t_len < hard_min:
                continue

            # user can define the hard minimum
            # slightly painful that this has to be calculated for each target
            min_score = max(hard_min, best_score, min(t_len - penalty_max, q_penalty_score))

            (matrix, max_score) = fill(target, query, min_score)

            # this will be used in min_score calc above on each loop
            if max_score > best_score:
                best_score = max_score
                result = []  # clear any old entries if the score increases

            if max_score >= min_score:
                result.append(
                    ScoreMatrix(
                        query=query,
                        target=target,
                        target_id=t_idx,
                        score=max_score,
                        matrix=matrix,
                        reversed=False,
                        original_seq=query,
                    )
                )

            if do_revcomp:
                (matrix, max_score) = fill(target, rev_query, min_score)
                if max_score >= min_score:
                    result.append(
                        ScoreMatrix(
                            query=rev_query,
                            target=target,
                            target_id=t_idx,
                            score=max_score,
                            matrix=matrix,
                            reversed=True,
                            original_seq=query,
                        )
                    )

        if len(result) == 0:
            unmapped.append(query)
            continue
        mapped.append(result)

    bt_mapped = []
    for r_set in mapped:
        clean_set = []
        for sm in r_set:
            bt = Backtrack(sm, match_type)
            if bt.pass_mode:
                if not keep_matrix:
                    bt.sm.matrix = None
                clean_set.append(bt)

        if len(clean_set) > 0:
            bt_mapped.append(clean_set)
        else:
            unmapped.append(r_set[0].original_seq)

    return AlignmentBatch(unmapped=unmapped, mapped=bt_mapped)


def fill(str target, str query, int min_score) -> Tuple[List[List[int]], int]:
    ## loosely based on https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
    cdef int t_len = len(target)
    cdef int q_len = len(query)
    cdef int max_score_seen = -1
    cdef int i_max = t_len - 1
    cdef int j_max = q_len - 1
    cdef int i, j, m, insert, delete, m_val

    f = [[]] * t_len  # list of lists Y-axis
    last_i = [0] * q_len
    i = 0
    while i < t_len:
        this_i = [0] * q_len
        first_j = True
        j = 0
        while j < q_len:
            # last scores
            if first_j:
                m = 0
                delete = GAP
                insert = GAP
                first_j = False
            else:
                m = last_i[j - 1]
                delete = last_i[j]
                insert = this_i[j - 1]
            # this score
            if target[i] == query[j]:
                m += MATCH
                # no else as mismatch += 0
            delete += GAP  # deletion
            insert += GAP  # insertion
            m_val = max([m, delete, insert])  # faster than sorted()[-1]
            this_i[j] = m_val
            if m_val > max_score_seen:
                max_score_seen = m_val
            if (j_max - j) + max_score_seen < min_score:
                break
            j += 1

        f[i] = this_i
        last_i = this_i
        if (i_max - i) + max_score_seen < min_score:
            break

        i += 1

    return (f, max_score_seen)
