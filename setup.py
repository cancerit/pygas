#!/usr/bin/env python3
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

from setuptools import setup
from Cython.Build import cythonize

config = {
    "name": "pygas",
    "description": "python Guide aligned Sequences",
    "long_description": open("README.md").read(),
    "author": "Keiran M Raine",
    "url": "https://github.com/cancerit/pygas",
    "author_email": "cgphelp@sanger.ac.uk",
    "version": "0.4.3",
    "license": "AGPL-3.0",
    "python_requires": ">= 3.9",
    "install_requires": ["click", "click-option-group"],
    "packages": ["pygas"],
    "setup_requires": ["click"],
    "test_suite": "tests",
    "tests_require": ["pytest"],
    "entry_points": {
        "console_scripts": ["pygas=pygas.cli:cli"],
    },
    "ext_modules": cythonize(
        ["pygas/matrix.pyx"],
        compiler_directives={"language_level": "3", "embedsignature": True},
    ),
}

setup(**config)
