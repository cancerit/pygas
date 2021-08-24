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
import click
from click_option_group import OptionGroup
import logging
import pkg_resources  # part of setuptools
from pygas.main import run as pygas_run

LOG_LEVELS = ("WARNING", "INFO", "DEBUG")


def _log_setup(loglevel):  # pragma: no cover
    logging.basicConfig(level=getattr(logging, loglevel.upper()), format="%(levelname)s: %(message)s")


def io_opts(function):  # pragma: no cover
    function = click.option(
        "-q",
        "--queries",
        required=True,
        type=click.Path(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
        ),
        help="Query/read sequences, one per line",
    )(function)
    function = click.option(
        "-t",
        "--targets",
        required=True,
        type=click.Path(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
        ),
        help="target/guide sequences, one per line",
    )(function)
    function = click.option(
        "-o",
        "--output",
        required=False,
        type=click.Path(
            exists=False,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
        ),
        help="Output to file, omit for stdout",
    )(function)
    function = click.option(
        "-m",
        "--minscore",
        required=False,
        type=int,
        default=15,
        help="Minimum score to retain, regardless of rule penalties.  Perfect match has score equal to query length.",
        show_default=True,
    )(function)
    function = click.option(
        "-r",
        "--rules",
        required=False,
        type=str,
        default=["M"],
        multiple=True,
        help="Rules for decreasing score based on allowed differences, e.g. MDI allows 1 base of mismatch, deletion and insertion.",
        show_default=True,
    )(function)
    return function


optgroup_debug = OptionGroup("\nDebug options", help="Options specific to troubleshooting, testing and debugging")


@click.group()
@click.version_option(pkg_resources.require(__name__.split(".")[0])[0].version)
def cli():  # pragma: no cover
    pass


@cli.command()
@io_opts
@click.option(
    "--rc",
    required=False,
    default=True,
    type=bool,
    help="Try both orientations of reads",
    show_default=True,
)
@optgroup_debug.option(
    "-l",
    "--loglevel",
    required=False,
    default="INFO",
    show_default=True,
    type=click.Choice(LOG_LEVELS, case_sensitive=False),
    help="Set logging verbosity",
)
def run(loglevel, targets, queries, output, minscore, rules, rc):  # pragma: no cover
    """
    Very basic command line for limited use cases, packages is intended to be used as an API
    """
    _log_setup(loglevel)
    pygas_run(targets, queries, output, minscore, rules, rc)
