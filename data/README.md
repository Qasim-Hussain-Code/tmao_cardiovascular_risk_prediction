# Data

No participant level data are stored in this repository, and none should be
added to it. The `.gitignore` blocks the tabular formats an export is likely
to arrive in, and it does so both by file extension and by directory, so
that a file saved in the wrong place is still caught.

## Directory layout

| Directory | Contents |
| --- | --- |
| `raw/` | The cohort export exactly as received, never edited by hand. |
| `interim/` | Intermediate files written during processing. Disposable. |
| `processed/` | Analysis ready tables produced by the pipeline. Disposable. |

## Supplying a cohort

Place the export at `data/raw/cohort.csv`, matching the schema in
`docs/data_dictionary.md`. The pipeline picks it up automatically on the next
run and reports its filename as the data provenance in every table, figure
and manifest it writes.

If that file is absent the pipeline simulates a cohort instead, and labels
every output as simulated. This is a convenience for developing and testing
the code. It is not a substitute for data, and a simulated run tells you
about the software rather than about TMAO.

## Provenance

When a real cohort is used, record here where it came from, the version or
extraction date of the export, the ethical approval under which it was
collected, and any restriction on redistribution. Anyone rerunning the
analysis later needs that information, and it is easier to write down at the
time than to reconstruct afterwards.
