#!/usr/bin/env bash
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

temp=$(mktemp -d)

# deletes the temp directory
function cleanup {
  rm -rf $temp
}
# register the cleanup function to be called on the EXIT signal
trap cleanup EXIT

pytest tests/ -x --cov-branch --cov-report term --cov-report html --cov=pygas --junitxml=junit.xml --cov-fail-under=88

set -e

pygas run -t examples/targets.txt.gz -q examples/queries.txt.gz -o ${temp}/large_check.tsv
zdiff -sq examples/result.tsv.gz ${temp}/large_check.tsv

pygas run -m 8 -t examples/multi_target.txt -q examples/multi_query.txt -o ${temp}/multi_check.tsv
zdiff -sq examples/multi_result.tsv.gz ${temp}/multi_check.tsv
