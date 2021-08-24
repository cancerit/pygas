# CHANGES

## 1.0.1 - 1.0.3

Teething issues with CI/CDand PyPi

## 1.0.0

First public release.

## Releases below this point are not publicly available

These items pre-date migration to github.  History has not been preserved.

## 0.4.3

- Various corrections on exact matching
- Addition of options to deal with boundary logic

## 0.4.2

Licensing

## 0.4.1

- Exact flag is deduced by absence of rules

## 0.4.0

- Exact matching also attempts substring matching for compatibility with SGE reads.
- Substring matching checked before falling through to matrix alignments.
- Minor efficiency improvements.

## 0.3.0

- Make revcomp of reads optional.
- Exact matching only option
- Exact matching speedup

## 0.2.5

CI test fix: Correct multi-hit test so matches usage on guides being unique regardless of revcomp

## 0.2.4

Reorganized exact matching code so that matrix match as last resort.

## 0.2.3

Handle when no rules are specified (no mismatch, ins or del)

## 0.2.2

Correction to the rule checks

## 0.2.1

Introduce use of score to allow:

- only best-score alignments to be captured
- matrix generation abandoned as soon as best score cannot be equalled.

## 0.2.0

Adds exact matching as a shortcut (not substr)

## 0.1.1

Add counts of D/I/M to the backtrack object to aid rule validation.

## 0.1.0

Base release for some initial manual use case testing.
