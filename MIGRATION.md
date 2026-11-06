# Migrating to N-TRIAGE 0.3

The software formerly named ACREAGE is now **N-TRIAGE**: nested tandem repeat
inference & analysis of genomic evolution.

New integrations should use:

- executable: `ntriage`
- C++ namespace and headers: `ntriage`, `include/ntriage/`
- CMake project and test prefix: `NTRIAGE`, `ntriage_*`
- bundled component references: `references/dfam_2.fasta`

The old executable and C++ include path are not built as aliases. Reconfigure a
fresh CMake build tree so cached target names do not leave a stale `acreage`
binary in place.

`ACRE`, `ACRO1`, `ACRE_000001`, `acre_id`, `complete_acres`, and `.acres.*` are
biological or output-schema terms and intentionally have not been renamed.
Historical analysis folders, frozen command records, and legacy schema fields
whose names contain `acreage` are retained as provenance; active commands now
invoke `ntriage`.

The main new analysis is enabled with `--deep-subfamilies`. For large cohorts,
combine it with `--pangenome` to use the bounded-batch streaming path. That path
retains only fixed-size sketches around called switches and writes structural QC
plus switch outputs; it intentionally skips the global feature, annotation, and
k-means products. See [README.md](README.md) for switch thresholds, output files,
and the distinction between whole-ACRE switch subfamilies and ACRO1 monomer
structural labels.
