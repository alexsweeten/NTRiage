#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace ntriage {

// A fixed-memory, phase-aware KMV (bottom-k) sketch. Canonical k-mers are paired
// with a coarse relative-position bin so tandem-repeat units with the same
// k-mer vocabulary but different internal order remain distinguishable. Both
// orientations are retained to keep comparisons reverse-complement invariant.
struct KmerSketch {
    uint64_t valid_kmers = 0;
    uint64_t distinct_kmers = 0;
    std::vector<uint64_t> hashes;
    uint64_t reverse_distinct_kmers = 0;
    std::vector<uint64_t> reverse_hashes;
};

// Build a deterministic sketch using a stable 64-bit hash. The temporary work
// is proportional to the sequence length; the retained two-orientation
// representation is at most 2 * sketch_size hashes (16 * sketch_size bytes,
// excluding vector overhead).
KmerSketch make_bottom_sketch(std::string_view sequence, unsigned k,
                              size_t sketch_size, size_t position_bins = 16);

// Estimate phase-aware set Jaccard using the smallest hashes in the coordinated
// union. Redundant forward/reverse samples are averaged before choosing the
// better relative orientation. If neither sketch was truncated, each
// directional estimate is exact.
double sketch_jaccard(const KmerSketch& left, const KmerSketch& right);

// Convert k-mer Jaccard to the Mash identity proxy 1 - D, where
// D = -log(2J / (1+J)) / k. The result is clamped to [0, 1].
double mash_ani(double jaccard, unsigned k);

struct SwitchThresholds {
    size_t flank_units = 3;
    size_t sketch_size = 512;
    double min_coherence = 0.90;
    double max_cross_ani = 0.92;
    double min_ani_drop = 0.04;
    double max_adjacent_ani = 0.92;
};

// boundary_index is the zero-based index of the first unit on the right side;
// equivalently, the boundary lies after unit boundary_index in one-based output.
struct SwitchBoundary {
    size_t boundary_index = 0;
    double adjacent_ani = 0.0;
    double left_flank_ani = 0.0;
    double right_flank_ani = 0.0;
    double minimum_flank_coherence = 0.0;
    double cross_flank_ani = 0.0;
    double ani_drop = 0.0;
    bool candidate = false;
};

// Evaluate every boundary having flank_units complete units on both sides.
// Pairwise flank summaries are medians: C(flank,2) comparisons within each
// side and flank^2 comparisons across sides. Boundaries involving an empty
// sketch are not evaluable and are omitted.
std::vector<SwitchBoundary> detect_switch_boundaries(
    const std::vector<KmerSketch>& sketches, unsigned k,
    const SwitchThresholds& thresholds);

// Greedy deterministic non-maximum suppression. Candidate boundaries are
// considered by decreasing ANI drop (then increasing boundary index); a weaker
// candidate is removed when it is fewer than min_separation units from a kept
// candidate. The returned boundaries are in genomic/unit order.
std::vector<SwitchBoundary> collapse_switch_boundaries(
    const std::vector<SwitchBoundary>& boundaries, size_t min_separation);

// Convenience overload using the flank width as the separation, matching the
// three-unit-flank switch definition.
std::vector<SwitchBoundary> collapse_switch_boundaries(
    const std::vector<SwitchBoundary>& boundaries,
    const SwitchThresholds& thresholds);

}  // namespace ntriage
