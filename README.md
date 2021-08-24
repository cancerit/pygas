# pyGaS

python Guide aligned Sequences

- [Command example](#command-example)
  - [Inputs](#inputs)
  - [Output format](#output-format)
- [Development](#development)
  - [Install](#install)
- [Testing](#testing)
  - [Local `venv` testing](#local-venv-testing)
  - [Local `pre-commit` hooks](#local-pre-commit-hooks)
  - [Docker testing](#docker-testing)
  - [CI tests](#ci-tests)

## Command example

The code is intended to be used as an API, not through this command line, however limited use is possible.

```bash
pygas run -t examples/targets.txt.gz -q examples/queries.txt.gz -o your_result.tsv
```

### Inputs

- `queries.txt`
  - A unique list of sequences (for performance reasons), one per line
    - This could be reworked to handle internally, however memory is a consideration
  - Matching sequences back to real input data and related information would be the responsibility of wrapping code
- `targets.txt`
  - One target sequence per line
  - Reverse compliment is handled automatically, see output format.
  - Targets need to be unique during mapping, expand out for things like dual guide permutations in your application

### Output format

Very simple text output of values that are available in API:

```text
#query	reversed	t_id	t_pos	cigar	seq	md	repeat_2-7...
AAAAATCGCTGCTACAGGT	False	48566	1	AAAAATCGCTGCTACAGGT	M19	19
CTGGTCTCGCACCCCAGGC	False	65601	1	CTGGTCTCGCACCCCAGGC	M19	18T
GGCGCGGTACTTGCCCAGA	False	34773	1	GGCGCGGTACTTGCCCAGA	S1M18	18
AAAAAAAAAAAAAAAAAAA	False	0	1	AAAAAAAAAAAAAAAAAAA	M19	19	True	1	1	TTTTTTTTTTTTTTTTTTT	M19	19
...
```

Where:

| Column     | Description                              | Interpretation                                                   |
|------------|------------------------------------------|------------------------------------------------------------------|
| `query`    | Original query sequence                  |                                                                  |
| `reversed` | Read was reversed to match the target    | following fields are based on this orientation                   |
| `t_id`     | ID of target mapped to                   | 0-based numbering in order targets passed                        |
| `t_pos`    | Start position within target sequence    | 1-based                                                          |
| `seq`      | Query in mapped orientation              | Corresponds to `cigar` and `md` orientation                      |
| `cigar`    | `cigar` string for use in SAM like files | For details see the [SAM specification][sam-spec]                |
| `md`       | `MD` string for use in SAM like files    | For details see the [SAM optional field specification][sam-opts] |

## Development

### Install

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python3 setup.py develop

# see later
pre-commit install

# remember to update requirements
pip freeze | grep -v virtualenv > requirements.txt
```

## Testing

There are 4 layers to testing and standards:

1. Local `venv` testing
1. Local `pre-commit` hooks
1. Tests embedded in `docker build`
1. `CI` tests

### Local `venv` testing

```bash
/tests/scripts/run_unit_tests.sh
```

### Local `pre-commit` hooks

This project additionally uses git pre-commit hooks via the [pre-commit tool](https://pre-commit.com/).  These are concerned
with file formats and standards, not the actual execution of code.  See `./.pre-commit-config.yaml`.

### Docker testing

The Docker build includes the unit tests, but removes many of the libraries before the final build stage.  Mainly for CI tests.

### CI tests

CI includes 2 additional tests, each based on the 2 datasets in the `./examples` directory.

```
Copyright (C) 2021

Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>

This file is part of pygas.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’.
```

<!-- refs -->

[sam-opts]: https://samtools.github.io/hts-specs/SAMtags.pdf
[sam-spec]: https://samtools.github.io/hts-specs/SAMv1.pdf
