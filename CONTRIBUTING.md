# Contributing

GraphSplit's maintained surface is the TOML interface, native `.gstt` format,
and the catalog/theta file contract documented in this repository. Changes
should preserve serial-ID joins and standard-library-only operation.

Before opening a pull request:

1. run `julia --project=. -e 'using Pkg; Pkg.test()'`;
2. add a focused test for numerical or format changes;
3. update the complete TOML and manual when adding an option;
4. report peak memory as well as runtime for large-catalog changes;
5. avoid committing generated `.gstt` tables or relocation output.

Travel-time format changes require a format-version bump and a documented
migration path. Changes to default inversion, filtering, or gauge settings
should include benchmark evidence rather than only a lower residual on one
catalog.
