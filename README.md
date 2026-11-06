# N-TRIAGE

**N-TRIAGE** (nested tandem repeat inference & analysis of genomic evolution)
is a dependency-light C++17 program for annotating, clustering, and detecting
sequence-state switches in nested tandem repeats. Input BED records are treated as **ACRE
arrays**, not as single ACREs. A fast lagged canonical-k-mer comparison estimates
the repeat period, then each array is split into complete 6.2–7.1 kb ACRE units
before annotation, feature construction, clustering, and optional deep subfamily
analysis around intra-array switches.

## Build

Requirements are a C++17 compiler, CMake 3.16+, and zlib.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

There are no runtime libraries beyond the platform C++ library and zlib.

## Run

```bash
./build/ntriage \
  --bed arrays.bed \
  --fasta genome.fa.gz \
  --refs references/dfam_2.fasta \
  --acro1-monomers acro1_monomers.bed \
  --out sample
```

The genome must have a samtools-compatible `.fai`. A compressed genome must be
BGZF and have a `.gzi`. Coordinates use standard zero-based, half-open BED
convention. Input BED rows are assigned `ARRAY_000001`, ... identifiers, and
the complete periods inferred inside them receive global `ACRE_000001`, ...
identifiers. Edge and short fragments are reported separately rather than
silently discarded.

Useful controls:

- `--min-period 6200`, `--max-period 7100`: permitted length of a complete ACRE.
- `--period-target 6650`: fallback period when sequence evidence is weak.
- `--boundary-span 1024`: bases sampled in the fast k-mer boundary comparison.
- `-k 15`: canonical k-mer length (valid range 1–31).
- `--window 200`: component annotation resolution in bases.
- `--bins 256`: hashed frequency-sketch dimensions.
- `-c N`: k-means cluster count; otherwise a deterministic 2–8 cluster heuristic.
- `-t N`: BGZF interval-fetch workers; auto mode is capped at eight.
- `--min-hits` and `--min-density`: annotation sensitivity controls.
- `--dump-exclusive`: emit the potentially large exclusive-k-mer table.
- `--write-fasta`: emit one sequence per complete ACRE.
- `--acro1-monomers FILE`: classify adjacent ACRO1 monomer starts as subfamily
  1, 2, or 3 and summarize the composition of every complete ACRE. This is
  optional and does not run an aligner.
- `--deep-subfamilies`: find coherent left/right sequence-state transitions and
  infer whole-ACRE subfamilies in the units surrounding each switch.
- `--sketch-size 512`: retained canonical bottom-k hashes per orientation and
  ACRE. Runtime is linear in input k-mers and retained memory is fixed per ACRE.
- `--sketch-position-bins 16`: coarse within-ACRE phase bins mixed into each
  k-mer token. This distinguishes repeats that share vocabulary but differ in
  internal order; both orientations are retained.
- `--switch-flank 3` and `--switch-radius 3`: units used to call a switch and
  units emitted around each nonredundant event.
- `--subfamily-ani 0.94`: minimum Mash-style phase-sketch identity proxy used
  by deterministic representative/leader clustering of switch-region ACREs.
- `--pangenome`: run the switch analysis as a bounded-batch stream. Each worker
  releases an array sequence after segmentation/sketching, and the process
  retains only fixed-size sketches in called switch contexts. This mode requires
  `--deep-subfamilies`; `--refs` is not needed.

## Deep switch-region subfamilies

`--deep-subfamilies` builds a deterministic phase-aware bottom-k sketch for every
complete ACRE, then evaluates boundaries with the same three-unit flank design used by
the research ANI workflow. The default primary call requires minimum within-flank
phase-sketch identity >=0.90, cross-flank phase-sketch identity <=0.92, a drop
>=0.04, and adjacent-unit phase-sketch identity <=0.92. Qualifying boundaries
less than three units apart are collapsed by keeping the strongest drop.

The reported identity values are Mash transforms of phase-binned canonical
k-mer Jaccard. They are a fast **phase-sketch identity proxy**, not conventional
full edit-distance ANI; thresholds copied from an alignment workflow should be
calibrated on the target data. This mode is intended to screen pangenomes for
candidate ANI discontinuities. Nearby ACREs are classified against the left and
right flank states and assigned a deterministic, order-stable
`switch_subfamily` by representative/leader clustering through an inverted
sketch index. These whole-ACRE states are deliberately separate from the
144/147-bp ACRO1 structural labels below.

Use a full sequence-alignment workflow when global edit-distance ANI
confirmation is required; the streaming sketch screen is not a replacement for
that exact analysis.

For a pangenome-scale switch analysis, use:

```bash
./build/ntriage \
  --bed arrays.bed --fasta pangenome.fa.gz \
  --out cohort --deep-subfamilies --pangenome --threads 8
```

The streaming path writes the six switch outputs plus `PREFIX.acres.bed`,
`PREFIX.partials.bed`, and `PREFIX.arrays.tsv`. It intentionally does not build
the all-unit annotation/feature/k-means outputs. Run without `--pangenome` when
those matrices or ACRO1 monomer classification are required.

## ACRO1 subfamily labels

When `--acro1-monomers FILE` is supplied, N-TRIAGE reads a BED3-or-BED6 file of
ACRO1 monomer anchors. BED6 strand is strongly recommended. Monomer names ending
in `_monomer_NUMBER` are grouped by the preceding name, which prevents anchors
from different source arrays from being paired. Otherwise, anchors are grouped
by chromosome and strand.

The deliberately simple numeric labels are:

| Label | Biological description | Fast classification rule |
|---|---|---|
| `subfamily_1` | canonical ACRO1 | exact 147-bp start-to-start spacing |
| `subfamily_2` | CCC/internal deletion | 144-bp spacing, closer to the ACRO1 template with positions 6–8 removed |
| `subfamily_3` | terminal/end deletion | 144-bp spacing, closer to the ACRO1 template with positions 145–147 removed |

Only exact 144- and 147-bp adjacent-anchor spacings receive these three labels.
Other spacings are excluded from the numeric composition rather than forced into
a class. The 144-bp comparison is a direct mismatch count against two templates
derived from the supplied 147-bp reference named `ACRO1`; ambiguous bases are
ignored. No dynamic programming or whole-sequence alignment is performed.

This preserves the important caveat from the primate analysis: subfamilies 2 and
3 are strand/phase-linked sequence forms. The numeric names should not by
themselves be interpreted as three independently evolved phylogenetic lineages.

## RGB component colors

`PREFIX.annotations.bed` is BED9+ and always has an `itemRgb` value in column 9.
Core defaults include red (`255,0,0`) for `ACRO1` and blue (`0,0,255`) for
`COMP-subunit_ACRO_rnd-5_family-38`. Other references receive stable,
name-derived colors. The complete mapping is written to `PREFIX.colors.tsv`.

Override any color with `--colors FILE`. Accepted tab-separated forms are:

```text
component	itemRgb
ACRO1	255,0,0
COMP-subunit_ACRO_rnd-5_family-38	0	0	255
```

Component names must match the reference FASTA. In a UCSC custom-track header,
set `itemRgb=On` to display column 9.

## Outputs

For `--out PREFIX`, a standard N-TRIAGE run writes:

- `PREFIX.acres.bed`: one complete inferred ACRE per line. BED columns 1–6 are
  followed by parent `array_id`, one-based `unit_in_array`, estimated array
  period, periodicity score, and the original BED name.
- `PREFIX.partials.bed`: material that cannot form a complete ACRE under the
  configured period bounds. Extra columns identify the parent array and whether
  it is a short array or an edge residual.
- `PREFIX.arrays.tsv`: one row per input BED interval, with complete-unit and
  residual counts/bases, segmentation method, estimated period, and periodicity
  score. For every array, `complete_bp + residual_bp = input_length`.
- `PREFIX.annotations.bed`: component calls for the complete ACRE units. Adjacent
  windows with the same component and orientation are merged.
- `PREFIX.features.tsv`: the numeric per-ACRE matrix used for clustering.
- `PREFIX.clusters.tsv`: per-ACRE cluster assignments, parent array, coordinates,
  periodicity metadata, dominant component, and k-mer QC statistics.
- `PREFIX.colors.tsv`: component-to-RGB mapping used in the annotation BED.
- `PREFIX.exclusive_kmers.tsv`: only with `--dump-exclusive`.
- `PREFIX.acres.fasta`: only with `--write-fasta`.
- `PREFIX.acro1_subfamilies.bed`: only with `--acro1-monomers`; one exact
  144/147-bp call per row with numeric label, strand, ACRE/array identifiers,
  classification method, and both 144-bp mismatch scores.
- `PREFIX.acro1_composition.tsv`: only with `--acro1-monomers`; one row per
  complete ACRE with counts, fractions, dominant numeric subfamily, pure/mixed
  composition class, and an orientation-aware run-length signature such as
  `S1x2|S2x1`.
- `PREFIX.switch_boundaries.tsv`: only with `--deep-subfamilies`; every
  evaluable boundary and its sketch-based flank/coherence statistics.
- `PREFIX.switches.tsv`: nonredundant primary switch events with genomic
  coordinates and phase-sketch identity-proxy metrics.
- `PREFIX.switches.bed`: one single-base BED marker per nonredundant switch,
  with the parent array and phase-sketch identity drop in extra columns.
- `PREFIX.switch_arrays.bed`: one deduplicated BED row per input array containing
  at least one switch, including switch count and event IDs.
- `PREFIX.switch_subfamilies.bed`: BED9+ switch-context ACREs, their left/right
  local state, global subfamily, confidence margin, and sketch QC.
- `PREFIX.switch_subfamily_summary.tsv`: global switch-subfamily representatives,
  unique-unit counts, arrays, and events.

In streaming `--pangenome` mode, only the three structural/QC outputs
(`acres.bed`, `partials.bed`, `arrays.tsv`) and the six `switch*` outputs are
written.

### Annotation BED columns

`PREFIX.annotations.bed` has 14 columns:

1. chromosome
2. start (zero-based)
3. end (half-open)
4. component name, or `unclassified`
5. BED score: hit density scaled to 0–1000
6. inferred component orientation (`+`, `-`, or `.`)
7. thick start (same as hit start)
8. thick end (same as hit end)
9. `itemRgb` color
10. child `acre_id`
11. parent `array_id`
12. weighted reference-k-mer hits
13. weighted hit density
14. confidence margin between the best and second-best component calls

## Feature and annotation model

A component-diagnostic k-mer occurs in exactly one supplied component reference.
In a standard run, an ACRE-exclusive k-mer occurs in exactly one **segmented
complete ACRE** in the current run; repeated occurrences within that unit are
allowed. The feature
vector combines component-diagnostic frequencies and presence flags, a
fixed-width hashed canonical-k-mer frequency sketch, and exclusive-k-mer
fraction/presence.

Annotation uses all component-reference k-mers. A k-mer shared by `n` references
contributes `1/n` to each, which prevents conserved transposable-element sequence
from dominating component-specific evidence. Period finding and annotation are
exact k-mer operations; no aligner or dynamic programming is used.

Only requested BGZF blocks are decompressed. Feature extraction and component
annotation run in parallel, annotation windows use zero-copy sequence views,
and k-means standardizes rows without duplicating the full feature matrix.
The dedicated `--pangenome` path processes at most `2 * threads` arrays per
batch. It does not materialize the global distinct-k-mer owner map, all-unit
sequence collection, annotation collection, or dense feature matrix. Only
called switch contexts survive batch release, at no more than
`16 * sketch-size` bytes of hash payload per retained ACRE (two orientations,
excluding small vector/metadata overhead). Peak memory therefore scales with a
bounded input batch plus the number of switch-context rows that must be emitted,
not with all pangenome ACREs. Runtime is linear in selected bases
apart from bounded period searches, local flank comparisons, hash-table
constants, and representative-clustering candidates.
