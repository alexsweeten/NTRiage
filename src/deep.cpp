#include "ntriage/deep.hpp"

#include "ntriage/kmer.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>

namespace ntriage {
namespace {

uint64_t mix64(uint64_t value) noexcept {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

bool sorted_unique(const std::vector<uint64_t>& values) noexcept {
    return std::adjacent_find(values.begin(), values.end(),
                              std::greater_equal<uint64_t>()) == values.end();
}

void validate_sketch(const KmerSketch& sketch) {
    if (sketch.hashes.size() > sketch.distinct_kmers) {
        throw std::invalid_argument("sketch retains more hashes than distinct k-mers");
    }
    if (!sorted_unique(sketch.hashes)) {
        throw std::invalid_argument("sketch hashes must be sorted and unique");
    }
    if (sketch.reverse_hashes.size() > sketch.reverse_distinct_kmers) {
        throw std::invalid_argument(
            "reverse sketch retains more hashes than distinct tokens");
    }
    if (!sorted_unique(sketch.reverse_hashes)) {
        throw std::invalid_argument("reverse sketch hashes must be sorted and unique");
    }
}

double directional_jaccard(const std::vector<uint64_t>& left, uint64_t left_distinct,
                           const std::vector<uint64_t>& right, uint64_t right_distinct) {
    const bool left_complete = left.size() == left_distinct;
    const bool right_complete = right.size() == right_distinct;
    if (left_complete && right_complete && left.empty() && right.empty()) return 1.0;

    size_t sample_size = 0;
    if (left_complete && right_complete) {
        if (left.size() > std::numeric_limits<size_t>::max() - right.size()) {
            throw std::overflow_error("sketch union is too large");
        }
        sample_size = left.size() + right.size();
    } else if (!left_complete && !right_complete) {
        if (left.size() != right.size()) {
            throw std::invalid_argument(
                "truncated sketches must use the same bottom-k size");
        }
        sample_size = left.size();
    } else {
        sample_size = left_complete ? right.size() : left.size();
    }

    size_t left_index = 0;
    size_t right_index = 0;
    size_t sampled_union = 0;
    size_t sampled_shared = 0;
    while (sampled_union < sample_size &&
           (left_index < left.size() || right_index < right.size())) {
        if (right_index >= right.size() ||
            (left_index < left.size() && left[left_index] < right[right_index])) {
            ++left_index;
        } else if (left_index >= left.size() || right[right_index] < left[left_index]) {
            ++right_index;
        } else {
            ++sampled_shared;
            ++left_index;
            ++right_index;
        }
        ++sampled_union;
    }
    if (sampled_union == 0U) return 0.0;
    return static_cast<double>(sampled_shared) / static_cast<double>(sampled_union);
}

double median(std::vector<double> values) {
    if (values.empty()) {
        throw std::invalid_argument("cannot take the median of an empty collection");
    }
    const size_t middle = values.size() / 2U;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if (values.size() % 2U != 0U) return upper;
    const double lower = *std::max_element(values.begin(), values.begin() + middle);
    return (lower + upper) / 2.0;
}

double sketch_ani(const KmerSketch& left, const KmerSketch& right, unsigned k) {
    return mash_ani(sketch_jaccard(left, right), k);
}

void validate_thresholds(const SwitchThresholds& thresholds) {
    if (thresholds.flank_units < 2U) {
        throw std::invalid_argument("switch flank must contain at least two units");
    }
    if (thresholds.flank_units > std::numeric_limits<size_t>::max() / 2U ||
        thresholds.flank_units >
            std::numeric_limits<size_t>::max() / thresholds.flank_units) {
        throw std::invalid_argument("switch flank is too large");
    }
    if (thresholds.sketch_size == 0U) {
        throw std::invalid_argument("switch sketch size must be positive");
    }
    const auto probability = [](double value) { return value >= 0.0 && value <= 1.0; };
    if (!probability(thresholds.min_coherence) ||
        !probability(thresholds.max_cross_ani) ||
        !probability(thresholds.min_ani_drop) ||
        !probability(thresholds.max_adjacent_ani)) {
        throw std::invalid_argument("switch thresholds must lie in [0, 1]");
    }
}

}  // namespace

KmerSketch make_bottom_sketch(std::string_view sequence, unsigned k,
                              size_t sketch_size, size_t position_bins) {
    if (k == 0U || k > 31U) {
        throw std::invalid_argument("k-mer size must be in [1, 31]");
    }
    if (sketch_size == 0U) {
        throw std::invalid_argument("sketch size must be positive");
    }
    if (position_bins == 0U) {
        throw std::invalid_argument("sketch position-bin count must be positive");
    }

    KmerSketch result;
    std::vector<uint64_t> hashes;
    std::vector<uint64_t> reverse_hashes;
    if (sequence.size() >= k) hashes.reserve(sequence.size() - k + 1U);
    if (sequence.size() >= k) reverse_hashes.reserve(sequence.size() - k + 1U);
    const size_t possible = sequence.size() >= k ? sequence.size() - k + 1U : 1U;
    for_each_canonical_kmer(sequence, k,
        [&](size_t position, uint64_t code, bool) {
            ++result.valid_kmers;
            const size_t forward_bin = std::min(
                position_bins - 1U, position * position_bins / possible);
            const size_t reverse_position = possible - 1U - position;
            const size_t reverse_bin = std::min(
                position_bins - 1U, reverse_position * position_bins / possible);
            hashes.push_back(mix64(code ^ mix64(forward_bin + 0x4e545249ULL)));
            reverse_hashes.push_back(
                mix64(code ^ mix64(reverse_bin + 0x4e545249ULL)));
        });

    std::sort(hashes.begin(), hashes.end());
    hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
    result.distinct_kmers = hashes.size();
    const size_t retained_forward = std::min(sketch_size, hashes.size());
    result.hashes.reserve(retained_forward);
    result.hashes.insert(result.hashes.end(), hashes.begin(),
                         hashes.begin() + static_cast<std::ptrdiff_t>(retained_forward));
    std::sort(reverse_hashes.begin(), reverse_hashes.end());
    reverse_hashes.erase(std::unique(reverse_hashes.begin(), reverse_hashes.end()),
                         reverse_hashes.end());
    result.reverse_distinct_kmers = reverse_hashes.size();
    const size_t retained_reverse = std::min(sketch_size, reverse_hashes.size());
    result.reverse_hashes.reserve(retained_reverse);
    result.reverse_hashes.insert(
        result.reverse_hashes.end(), reverse_hashes.begin(),
        reverse_hashes.begin() + static_cast<std::ptrdiff_t>(retained_reverse));
    return result;
}

double sketch_jaccard(const KmerSketch& left, const KmerSketch& right) {
    validate_sketch(left);
    validate_sketch(right);
    const double forward_forward = directional_jaccard(
        left.hashes, left.distinct_kmers, right.hashes, right.distinct_kmers);
    if (!left.reverse_hashes.empty() && !right.reverse_hashes.empty()) {
        const double reverse_reverse = directional_jaccard(
            left.reverse_hashes, left.reverse_distinct_kmers,
            right.reverse_hashes, right.reverse_distinct_kmers);
        const double forward_reverse = directional_jaccard(
            left.hashes, left.distinct_kmers,
            right.reverse_hashes, right.reverse_distinct_kmers);
        const double reverse_forward = directional_jaccard(
            left.reverse_hashes, left.reverse_distinct_kmers,
            right.hashes, right.distinct_kmers);
        // Mirrored phase bins create two coordinated samples for each relative
        // orientation. Average redundant samples before choosing orientation;
        // taking the maximum of all four would bias truncated sketches upward.
        const double same_orientation = (forward_forward + reverse_reverse) / 2.0;
        const double opposite_orientation = (forward_reverse + reverse_forward) / 2.0;
        return std::max(same_orientation, opposite_orientation);
    }
    return forward_forward;
}

double mash_ani(double jaccard, unsigned k) {
    if (k == 0U) throw std::invalid_argument("Mash k-mer size must be positive");
    if (!std::isfinite(jaccard) || jaccard < 0.0 || jaccard > 1.0) {
        throw std::invalid_argument("Jaccard similarity must lie in [0, 1]");
    }
    if (jaccard <= 0.0) return 0.0;
    if (jaccard >= 1.0) return 1.0;
    const double shared_fraction = (2.0 * jaccard) / (1.0 + jaccard);
    const double distance = -std::log(shared_fraction) / static_cast<double>(k);
    return std::max(0.0, std::min(1.0, 1.0 - distance));
}

std::vector<SwitchBoundary> detect_switch_boundaries(
    const std::vector<KmerSketch>& sketches, unsigned k,
    const SwitchThresholds& thresholds) {
    if (k == 0U || k > 31U) {
        throw std::invalid_argument("switch k-mer size must be in [1, 31]");
    }
    validate_thresholds(thresholds);
    for (const auto& sketch : sketches) {
        validate_sketch(sketch);
        if (sketch.hashes.size() != sketch.distinct_kmers &&
            sketch.hashes.size() != thresholds.sketch_size) {
            throw std::invalid_argument(
                "truncated sketch size does not match switch thresholds");
        }
    }

    const size_t flank = thresholds.flank_units;
    if (sketches.size() < 2U * flank) return {};

    // Adjacent boundary windows overlap heavily. Cache each pair out to the
    // largest possible distance in a 2*flank window so every sketch comparison
    // is merged at most once rather than once per boundary that contains it.
    const size_t maximum_delta = 2U * flank - 1U;
    if (sketches.size() > std::numeric_limits<size_t>::max() / maximum_delta) {
        throw std::overflow_error("switch similarity cache is too large");
    }
    const double missing = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> similarities(sketches.size() * maximum_delta, missing);
    auto pair_ani = [&](size_t first, size_t second) -> double {
        if (second < first) std::swap(first, second);
        if (first == second) return 1.0;
        const size_t delta = second - first;
        if (delta > maximum_delta) {
            throw std::logic_error("switch comparison exceeds cache width");
        }
        double& value = similarities[first * maximum_delta + delta - 1U];
        if (std::isnan(value)) value = sketch_ani(sketches[first], sketches[second], k);
        return value;
    };

    std::vector<SwitchBoundary> result;
    result.reserve(sketches.size() - 2U * flank + 1U);
    for (size_t boundary = flank; boundary + flank <= sketches.size(); ++boundary) {
        bool evaluable = true;
        for (size_t unit = boundary - flank; unit < boundary + flank; ++unit) {
            if (sketches[unit].hashes.empty()) {
                evaluable = false;
                break;
            }
        }
        if (!evaluable) continue;

        std::vector<double> left_values;
        std::vector<double> right_values;
        std::vector<double> cross_values;
        left_values.reserve(flank * (flank - 1U) / 2U);
        right_values.reserve(flank * (flank - 1U) / 2U);
        cross_values.reserve(flank * flank);

        const size_t left_begin = boundary - flank;
        const size_t right_end = boundary + flank;
        for (size_t first = left_begin; first < boundary; ++first) {
            for (size_t second = first + 1U; second < boundary; ++second) {
                left_values.push_back(pair_ani(first, second));
            }
        }
        for (size_t first = boundary; first < right_end; ++first) {
            for (size_t second = first + 1U; second < right_end; ++second) {
                right_values.push_back(pair_ani(first, second));
            }
        }
        for (size_t left_index = left_begin; left_index < boundary; ++left_index) {
            for (size_t right_index = boundary; right_index < right_end; ++right_index) {
                cross_values.push_back(pair_ani(left_index, right_index));
            }
        }

        SwitchBoundary call;
        call.boundary_index = boundary;
        call.adjacent_ani = pair_ani(boundary - 1U, boundary);
        call.left_flank_ani = median(std::move(left_values));
        call.right_flank_ani = median(std::move(right_values));
        call.minimum_flank_coherence =
            std::min(call.left_flank_ani, call.right_flank_ani);
        call.cross_flank_ani = median(std::move(cross_values));
        call.ani_drop = call.minimum_flank_coherence - call.cross_flank_ani;
        call.candidate =
            call.minimum_flank_coherence >= thresholds.min_coherence &&
            call.cross_flank_ani <= thresholds.max_cross_ani &&
            call.ani_drop >= thresholds.min_ani_drop &&
            call.adjacent_ani <= thresholds.max_adjacent_ani;
        result.push_back(call);
    }
    return result;
}

std::vector<SwitchBoundary> collapse_switch_boundaries(
    const std::vector<SwitchBoundary>& boundaries, size_t min_separation) {
    if (min_separation == 0U) {
        throw std::invalid_argument("switch separation must be positive");
    }

    std::vector<SwitchBoundary> ranked;
    ranked.reserve(boundaries.size());
    for (const auto& boundary : boundaries) {
        if (boundary.candidate) ranked.push_back(boundary);
    }
    std::sort(ranked.begin(), ranked.end(), [](const auto& left, const auto& right) {
        if (left.ani_drop != right.ani_drop) return left.ani_drop > right.ani_drop;
        return left.boundary_index < right.boundary_index;
    });

    std::vector<SwitchBoundary> kept;
    kept.reserve(ranked.size());
    for (const auto& candidate : ranked) {
        const bool separated = std::all_of(
            kept.begin(), kept.end(), [&](const SwitchBoundary& other) {
                const size_t distance = candidate.boundary_index > other.boundary_index
                    ? candidate.boundary_index - other.boundary_index
                    : other.boundary_index - candidate.boundary_index;
                return distance >= min_separation;
            });
        if (separated) kept.push_back(candidate);
    }
    std::sort(kept.begin(), kept.end(), [](const auto& left, const auto& right) {
        return left.boundary_index < right.boundary_index;
    });
    return kept;
}

std::vector<SwitchBoundary> collapse_switch_boundaries(
    const std::vector<SwitchBoundary>& boundaries,
    const SwitchThresholds& thresholds) {
    if (thresholds.flank_units == 0U) {
        throw std::invalid_argument("switch flank must be positive");
    }
    return collapse_switch_boundaries(boundaries, thresholds.flank_units);
}

}  // namespace ntriage
