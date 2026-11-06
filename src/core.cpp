#include "ntriage/core.hpp"

#include "ntriage/deep.hpp"
#include "ntriage/fasta.hpp"
#include "ntriage/kmer.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cctype>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace ntriage {
namespace {

using Clock = std::chrono::steady_clock;

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    size_t begin = 0;
    while (true) {
        const size_t tab = line.find('\t', begin);
        fields.emplace_back(line.substr(begin, tab - begin));
        if (tab == std::string::npos) break;
        begin = tab + 1;
    }
    return fields;
}

uint64_t parse_coordinate(const std::string& text, const std::string& context) {
    if (text.empty() || text.front() == '-') {
        throw std::runtime_error("invalid BED coordinate in " + context + ": " + text);
    }
    size_t used = 0;
    uint64_t value = 0;
    try {
        value = std::stoull(text, &used);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid BED coordinate in " + context + ": " + text);
    }
    if (used != text.size()) {
        throw std::runtime_error("invalid BED coordinate in " + context + ": " + text);
    }
    return value;
}

unsigned popcount64(uint64_t value) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    return static_cast<unsigned>(__builtin_popcountll(value));
#else
    unsigned count = 0;
    while (value) { value &= value - 1U; ++count; }
    return count;
#endif
}

unsigned trailing_zeroes64(uint64_t value) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    return static_cast<unsigned>(__builtin_ctzll(value));
#else
    unsigned count = 0;
    while (!(value & 1U)) { value >>= 1U; ++count; }
    return count;
#endif
}

uint64_t mix64(uint64_t value) noexcept {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

struct AnnotationSegment {
    size_t relative_start = 0;
    size_t relative_end = 0;
    int reference = -1;
    char strand = '.';
    double weighted_hits = 0;
    uint64_t valid_kmers = 0;
    double confidence_sum = 0;
};

struct RegionStats {
    uint64_t valid_kmers = 0;
    uint64_t distinct_kmers = 0;
    uint64_t exclusive_kmers = 0;
    std::string top_component = "unclassified";
    double component_coverage = 0;
};

struct Rgb {
    unsigned red = 0;
    unsigned green = 0;
    unsigned blue = 0;

    std::string text() const {
        return std::to_string(red) + "," + std::to_string(green) + "," +
               std::to_string(blue);
    }
};

struct AcreUnit {
    std::string chrom;
    uint64_t start = 0;
    uint64_t end = 0;
    uint64_t array_start = 0;
    uint64_t array_end = 0;
    std::string source_name;
    std::string id;
    std::string array_id;
    size_t unit_number = 0;
    size_t estimated_period = 0;
    double periodicity_score = 0;
    std::string sequence;
};

struct ResidualFragment {
    std::string chrom;
    uint64_t start = 0;
    uint64_t end = 0;
    std::string array_id;
    std::string reason;
};

struct ArraySegmentation {
    std::vector<std::pair<size_t, size_t>> units;
    std::vector<std::pair<size_t, size_t>> residuals;
    size_t estimated_period = 0;
    double periodicity_score = 0;
    std::string method = "no_complete_acre";
};

struct Acro1MonomerAnchor {
    std::string chrom;
    uint64_t start = 0;
    uint64_t end = 0;
    std::string group;
    char strand = '.';
};

struct Acro1SubfamilyCall {
    uint64_t start = 0;
    uint64_t end = 0;
    size_t unit_index = 0;
    unsigned subfamily = 0;
    char strand = '.';
    size_t spacing = 0;
    size_t internal_mismatches = 0;
    size_t terminal_mismatches = 0;
    size_t compared_bases = 0;
    std::string method;
};

struct DeepBoundaryRecord {
    std::string array_id;
    std::string source_name;
    std::string chrom;
    uint64_t array_start = 0;
    uint64_t array_end = 0;
    size_t boundary_after_unit = 0;
    size_t boundary_before_unit = 0;
    uint64_t boundary_position = 0;
    SwitchBoundary metrics;
};

struct DeepSwitchEvent {
    std::string id;
    size_t array_begin = 0;
    size_t array_end = 0;
    size_t boundary_global = 0;
    DeepBoundaryRecord record;
};

struct StreamingArrayAnalysis {
    uint64_t input_bases = 0;
    uint64_t complete_bases = 0;
    uint64_t residual_bases = 0;
    size_t complete_units = 0;
    size_t retained_hashes = 0;
    ArraySegmentation segmentation;
    std::vector<DeepBoundaryRecord> boundaries;
    std::vector<SwitchBoundary> switches;
    std::vector<size_t> retained_local_indices;
    std::vector<AcreUnit> retained_units;
    std::vector<KmerSketch> retained_sketches;
};

char complement_base(char base) {
    switch (base) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default: return 'N';
    }
}

std::string reverse_complement(const std::string& sequence) {
    std::string result(sequence.size(), 'N');
    for (size_t i = 0; i < sequence.size(); ++i) {
        result[sequence.size() - i - 1U] = complement_base(sequence[i]);
    }
    return result;
}

std::vector<Acro1MonomerAnchor> read_acro1_monomers(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open ACRO1 monomer BED: " + path);
    std::vector<Acro1MonomerAnchor> anchors;
    std::string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line.front() == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 3U) {
            throw std::runtime_error("ACRO1 monomer BED line " +
                                     std::to_string(line_number) + " has fewer than 3 columns");
        }
        Acro1MonomerAnchor anchor;
        anchor.chrom = fields[0];
        anchor.start = parse_coordinate(fields[1], path + " line " + std::to_string(line_number));
        anchor.end = parse_coordinate(fields[2], path + " line " + std::to_string(line_number));
        if (anchor.end < anchor.start) {
            throw std::runtime_error("ACRO1 monomer BED end precedes start at line " +
                                     std::to_string(line_number));
        }
        const std::string name = fields.size() >= 4U ? fields[3] : anchor.chrom;
        const size_t marker = name.rfind("_monomer_");
        anchor.group = marker == std::string::npos ? anchor.chrom : name.substr(0, marker);
        if (fields.size() >= 6U && !fields[5].empty() &&
            (fields[5][0] == '+' || fields[5][0] == '-')) {
            anchor.strand = fields[5][0];
        }
        anchors.push_back(std::move(anchor));
    }
    if (anchors.empty()) throw std::runtime_error("ACRO1 monomer BED is empty: " + path);
    return anchors;
}

std::pair<size_t, size_t> mismatch_count(const std::string& query,
                                         const std::string& target) {
    if (query.size() != target.size()) {
        throw std::runtime_error("internal ACRO1 template length mismatch");
    }
    size_t mismatches = 0;
    size_t compared = 0;
    for (size_t i = 0; i < query.size(); ++i) {
        const char q = static_cast<char>(std::toupper(static_cast<unsigned char>(query[i])));
        const char t = static_cast<char>(std::toupper(static_cast<unsigned char>(target[i])));
        if ((q != 'A' && q != 'C' && q != 'G' && q != 'T') ||
            (t != 'A' && t != 'C' && t != 'G' && t != 'T')) continue;
        ++compared;
        if (q != t) ++mismatches;
    }
    return {mismatches, compared};
}

unsigned classify_acro1_144(const std::string& query, const std::string& internal_template,
                            const std::string& terminal_template, size_t& internal_mismatches,
                            size_t& terminal_mismatches, size_t& compared_bases) {
    const auto internal = mismatch_count(query, internal_template);
    const auto terminal = mismatch_count(query, terminal_template);
    internal_mismatches = internal.first;
    terminal_mismatches = terminal.first;
    compared_bases = std::min(internal.second, terminal.second);
    return internal.first <= terminal.first ? 2U : 3U;
}

std::string acro1_composition_class(const std::array<size_t, 3>& counts) {
    std::vector<unsigned> present;
    for (unsigned i = 0; i < counts.size(); ++i) if (counts[i]) present.push_back(i + 1U);
    if (present.empty()) return "no_subfamily_call";
    if (present.size() == 1U) return "pure_subfamily_" + std::to_string(present.front());
    std::string result = "mixed";
    for (unsigned label : present) result += "_subfamily_" + std::to_string(label);
    return result;
}

std::string acro1_architecture(const std::vector<const Acro1SubfamilyCall*>& calls,
                               char strand) {
    if (calls.empty()) return "";
    std::vector<const Acro1SubfamilyCall*> ordered = calls;
    std::sort(ordered.begin(), ordered.end(), [strand](const auto* left, const auto* right) {
        return strand == '-' ? left->start > right->start : left->start < right->start;
    });
    std::ostringstream result;
    unsigned previous = 0;
    size_t run = 0;
    auto flush = [&] {
        if (!run) return;
        if (result.tellp() > 0) result << '|';
        result << 'S' << previous << 'x' << run;
    };
    for (const auto* call : ordered) {
        if (call->subfamily == previous) {
            ++run;
        } else {
            flush();
            previous = call->subfamily;
            run = 1;
        }
    }
    flush();
    return result.str();
}

void write_acro1_subfamilies(const Options& options, const ReferenceIndex& reference_index,
                             const std::vector<AcreUnit>& units) {
    size_t acro1_reference = reference_index.references.size();
    for (size_t i = 0; i < reference_index.references.size(); ++i) {
        if (reference_index.references[i].name == "ACRO1") {
            acro1_reference = i;
            break;
        }
    }
    if (acro1_reference == reference_index.references.size()) {
        throw std::runtime_error("--acro1-monomers requires a reference named ACRO1");
    }
    const std::string canonical = reference_index.references[acro1_reference].sequence;
    if (canonical.size() != 147U) {
        throw std::runtime_error("--acro1-monomers requires the ACRO1 reference to be 147 bp");
    }
    const std::string internal_template = canonical.substr(0, 5) + canonical.substr(8);
    const std::string terminal_template = canonical.substr(0, 144);

    std::unordered_map<std::string, std::vector<size_t>> units_by_chrom;
    for (size_t i = 0; i < units.size(); ++i) units_by_chrom[units[i].chrom].push_back(i);
    for (auto& item : units_by_chrom) {
        std::sort(item.second.begin(), item.second.end(), [&](size_t left, size_t right) {
            return units[left].start < units[right].start;
        });
    }
    auto locate_unit = [&](const std::string& chrom, uint64_t midpoint) -> size_t {
        const auto found = units_by_chrom.find(chrom);
        if (found == units_by_chrom.end()) return units.size();
        const auto& indices = found->second;
        auto position = std::upper_bound(indices.begin(), indices.end(), midpoint,
            [&](uint64_t value, size_t index) { return value < units[index].start; });
        if (position == indices.begin()) return units.size();
        const size_t index = *--position;
        return midpoint < units[index].end ? index : units.size();
    };

    auto anchors = read_acro1_monomers(options.acro1_monomers_path);
    std::unordered_map<std::string, std::vector<Acro1MonomerAnchor>> groups;
    for (auto& anchor : anchors) {
        const std::string key = anchor.chrom + "\t" + anchor.group + "\t" + anchor.strand;
        groups[key].push_back(std::move(anchor));
    }
    std::vector<Acro1SubfamilyCall> calls;
    size_t skipped_outside_acres = 0;
    size_t skipped_other_spacing = 0;
    for (auto& item : groups) {
        auto& group = item.second;
        std::sort(group.begin(), group.end(), [](const auto& left, const auto& right) {
            return left.start < right.start;
        });
        for (size_t i = 1; i < group.size(); ++i) {
            const auto& left = group[i - 1U];
            const auto& right = group[i];
            const uint64_t spacing = right.start - left.start;
            if (spacing != 144U && spacing != 147U) {
                ++skipped_other_spacing;
                continue;
            }
            const size_t unit_index = locate_unit(left.chrom, left.start + spacing / 2U);
            if (unit_index == units.size() || left.start < units[unit_index].start ||
                right.start > units[unit_index].end) {
                ++skipped_outside_acres;
                continue;
            }
            Acro1SubfamilyCall call;
            call.start = left.start;
            call.end = right.start;
            call.unit_index = unit_index;
            call.strand = left.strand;
            call.spacing = static_cast<size_t>(spacing);
            if (spacing == 147U) {
                call.subfamily = 1;
                call.method = "exact_147_spacing";
                call.compared_bases = 147;
            } else {
                std::string query = units[unit_index].sequence.substr(
                    static_cast<size_t>(left.start - units[unit_index].start), 144U);
                if (left.strand == '-') query = reverse_complement(query);
                call.subfamily = classify_acro1_144(
                    query, internal_template, terminal_template, call.internal_mismatches,
                    call.terminal_mismatches, call.compared_bases);
                if (call.subfamily == 2U) {
                    call.method = call.internal_mismatches == call.terminal_mismatches
                                      ? "exact_144_tie_assigned_internal"
                                      : "exact_144_internal_deletion_template";
                } else {
                    call.method = "exact_144_terminal_deletion_template";
                }
            }
            calls.push_back(std::move(call));
        }
    }

    const std::array<Rgb, 3> colors = {{{78, 121, 167}, {228, 26, 28}, {255, 157, 0}}};
    std::ofstream bed(options.output_prefix + ".acro1_subfamilies.bed");
    std::ofstream composition(options.output_prefix + ".acro1_composition.tsv");
    if (!bed || !composition) {
        throw std::runtime_error("cannot create ACRO1 subfamily outputs with prefix: " +
                                 options.output_prefix);
    }
    bed << "#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb"
           "\tacre_id\tarray_id\tunit_in_array\tspacing\tsubfamily_id\tsubfamily_label"
           "\tbiological_name"
           "\tclassification_method\tinternal_mismatches\tterminal_mismatches"
           "\tcompared_bases\n";
    const std::array<const char*, 3> biological_names = {{
        "canonical", "ccc_removed", "end_removed"
    }};
    std::vector<std::vector<const Acro1SubfamilyCall*>> by_unit(units.size());
    for (const auto& call : calls) {
        by_unit[call.unit_index].push_back(&call);
        const size_t difference = call.internal_mismatches > call.terminal_mismatches
                                      ? call.internal_mismatches - call.terminal_mismatches
                                      : call.terminal_mismatches - call.internal_mismatches;
        const int score = call.subfamily == 1 ? 1000 : static_cast<int>(std::min<size_t>(
            1000U, call.compared_bases ? 1000U * difference / call.compared_bases : 0U));
        const auto& unit = units[call.unit_index];
        bed << unit.chrom << '\t' << call.start << '\t' << call.end << "\tsubfamily_"
            << call.subfamily << '\t' << score << '\t' << call.strand << '\t'
            << call.start << '\t' << call.end << '\t' << colors[call.subfamily - 1U].text()
            << '\t' << unit.id << '\t' << unit.array_id << '\t' << unit.unit_number
            << '\t' << call.spacing << '\t' << call.subfamily << "\tsubfamily_"
            << call.subfamily << '\t' << biological_names[call.subfamily - 1U] << '\t'
            << call.method << '\t' << call.internal_mismatches
            << '\t' << call.terminal_mismatches << '\t' << call.compared_bases << '\n';
    }

    composition << "acre_id\tarray_id\tunit_in_array\tchrom\tstart\tend"
                   "\tsubfamily_1\tsubfamily_2\tsubfamily_3\tclassified_total"
                   "\tfraction_subfamily_1\tfraction_subfamily_2\tfraction_subfamily_3"
                   "\tdominant_subfamily\tcomposition_class\tarchitecture_rle\n"
                << std::setprecision(9);
    for (size_t i = 0; i < units.size(); ++i) {
        std::array<size_t, 3> counts{{0, 0, 0}};
        size_t plus = 0;
        size_t minus = 0;
        for (const auto* call : by_unit[i]) {
            ++counts[call->subfamily - 1U];
            if (call->strand == '-') ++minus; else ++plus;
        }
        const size_t total = counts[0] + counts[1] + counts[2];
        size_t maximum = 0;
        for (size_t value : counts) maximum = std::max(maximum, value);
        std::string dominant = "none";
        if (maximum) {
            std::vector<unsigned> tied;
            for (unsigned label = 0; label < counts.size(); ++label) {
                if (counts[label] == maximum) tied.push_back(label + 1U);
            }
            dominant = tied.size() == 1U ? "subfamily_" + std::to_string(tied.front()) : "tie";
        }
        const auto& unit = units[i];
        composition << unit.id << '\t' << unit.array_id << '\t' << unit.unit_number << '\t'
                    << unit.chrom << '\t' << unit.start << '\t' << unit.end;
        for (size_t value : counts) composition << '\t' << value;
        composition << '\t' << total;
        for (size_t value : counts) {
            if (total) composition << '\t' << static_cast<double>(value) / total;
            else composition << "\tNA";
        }
        const char orientation = minus > plus ? '-' : '+';
        composition << '\t' << dominant << '\t' << acro1_composition_class(counts) << '\t'
                    << acro1_architecture(by_unit[i], orientation) << '\n';
    }
    std::cerr << "Classified " << calls.size() << " ACRO1 units as subfamily 1/2/3; "
              << skipped_other_spacing << " adjacent pairs had spacing other than 144/147 and "
              << skipped_outside_acres << " fell outside complete ACRE boundaries\n";
}

double boundary_similarity(const std::vector<uint64_t>& codes, size_t first,
                           size_t second, size_t span) {
    constexpr uint64_t invalid = std::numeric_limits<uint64_t>::max();
    constexpr size_t stride = 4;
    size_t matches = 0;
    size_t valid = 0;
    for (size_t offset = 0; offset < span && first + offset < codes.size() &&
                            second + offset < codes.size(); offset += stride) {
        const uint64_t a = codes[first + offset];
        const uint64_t b = codes[second + offset];
        if (a == invalid || b == invalid) continue;
        ++valid;
        if (a == b) ++matches;
    }
    return valid ? static_cast<double>(matches) / static_cast<double>(valid) : 0.0;
}

ArraySegmentation segment_sequence(const std::string& sequence, const Options& options) {
    ArraySegmentation result;
    const size_t length = sequence.size();
    result.estimated_period = options.period_target;

    constexpr uint64_t invalid = std::numeric_limits<uint64_t>::max();
    std::vector<uint64_t> codes(length, invalid);
    for_each_canonical_kmer(sequence, options.k,
        [&](size_t position, uint64_t code, bool) { codes[position] = code; });

    // Estimate the repeat length before deciding how many ACREs are present.
    // Comparing canonical k-mers at a lag is much cheaper than alignment and
    // avoids forcing every array to the midpoint of the allowed period range.
    double period_evidence = -1.0;
    if (length >= options.min_period * 2U) {
        const size_t upper_period = std::min(options.max_period, length - options.min_period);
        for (size_t period = options.min_period; period <= upper_period; ++period) {
            const double score = boundary_similarity(codes, 0, period, options.boundary_span);
            const size_t candidate_distance = period > options.period_target
                                                  ? period - options.period_target
                                                  : options.period_target - period;
            const size_t best_distance = result.estimated_period > options.period_target
                                             ? result.estimated_period - options.period_target
                                             : options.period_target - result.estimated_period;
            if (score > period_evidence + 1e-12 ||
                (std::abs(score - period_evidence) <= 1e-12 &&
                 candidate_distance < best_distance)) {
                period_evidence = score;
                result.estimated_period = period;
            }
        }
        // A handful of matching canonical k-mers is already strong evidence at this
        // scale; weaker maxima are treated as no detectable periodic signal.
        if (period_evidence < 0.03) result.estimated_period = options.period_target;
    }

    const size_t minimum_units = (length + options.max_period - 1U) / options.max_period;
    const size_t maximum_units = length / options.min_period;
    size_t unit_count = 0;
    size_t covered_length = 0;

    if (maximum_units >= minimum_units && maximum_units > 0) {
        const size_t rounded = static_cast<size_t>(
            std::llround(static_cast<double>(length) /
                         static_cast<double>(result.estimated_period)));
        unit_count = std::max(minimum_units, std::min(maximum_units, rounded));
        covered_length = length;
        result.method = unit_count > 1 && period_evidence >= 0.03
                            ? "self_similarity"
                            : "target_length";
    } else {
        unit_count = length / result.estimated_period;
        covered_length = unit_count * result.estimated_period;
        if (unit_count > 0) {
            result.method = unit_count > 1 && period_evidence >= 0.03
                                ? "self_similarity_with_residual"
                                : "target_length_with_residual";
        }
    }

    if (unit_count == 0) {
        if (length > 0) result.residuals.emplace_back(0, length);
        return result;
    }

    std::vector<size_t> boundaries;
    boundaries.reserve(unit_count + 1U);
    boundaries.push_back(0);
    std::vector<double> boundary_scores;
    boundary_scores.reserve(unit_count > 0 ? unit_count - 1U : 0U);
    for (size_t unit = 1; unit < unit_count; ++unit) {
        const size_t previous = boundaries.back();
        const size_t remaining_units = unit_count - unit;
        size_t lower = previous + options.min_period;
        const size_t remaining_maximum = remaining_units * options.max_period;
        if (covered_length > remaining_maximum) {
            lower = std::max(lower, covered_length - remaining_maximum);
        }
        const size_t remaining_minimum = remaining_units * options.min_period;
        const size_t upper = std::min(previous + options.max_period,
                                      covered_length - remaining_minimum);
        if (lower > upper) throw std::runtime_error("internal ACRE segmentation error");

        const size_t ideal = static_cast<size_t>(std::llround(
            static_cast<double>(covered_length) * static_cast<double>(unit) /
            static_cast<double>(unit_count)));
        const size_t boundary_slop = std::max<size_t>(
            1U, (options.max_period - options.min_period) / 4U);
        const size_t ideal_lower = ideal > boundary_slop ? ideal - boundary_slop : 0U;
        lower = std::max(lower, ideal_lower);
        const size_t ideal_upper = ideal > std::numeric_limits<size_t>::max() - boundary_slop
                                       ? std::numeric_limits<size_t>::max()
                                       : ideal + boundary_slop;
        const size_t constrained_upper = std::min(upper, ideal_upper);
        if (lower > constrained_upper) {
            throw std::runtime_error("internal ACRE boundary constraint failure");
        }
        size_t best_boundary = lower;
        double best_score = -1.0;
        for (size_t candidate = lower; candidate <= constrained_upper; ++candidate) {
            const double score = boundary_similarity(codes, previous, candidate,
                                                       options.boundary_span);
            const size_t candidate_distance = candidate > ideal ? candidate - ideal
                                                                 : ideal - candidate;
            const size_t best_distance = best_boundary > ideal ? best_boundary - ideal
                                                                : ideal - best_boundary;
            if (score > best_score + 1e-12 ||
                (std::abs(score - best_score) <= 1e-12 &&
                 candidate_distance < best_distance)) {
                best_score = score;
                best_boundary = candidate;
            }
        }
        boundaries.push_back(best_boundary);
        boundary_scores.push_back(best_score);
    }
    boundaries.push_back(covered_length);
    for (size_t unit = 0; unit < unit_count; ++unit) {
        const size_t unit_length = boundaries[unit + 1U] - boundaries[unit];
        if (unit_length < options.min_period || unit_length > options.max_period) {
            throw std::runtime_error("internal ACRE period constraint failure");
        }
        result.units.emplace_back(boundaries[unit], boundaries[unit + 1U]);
    }
    if (covered_length < length) result.residuals.emplace_back(covered_length, length);
    if (!boundary_scores.empty()) {
        result.periodicity_score = std::accumulate(boundary_scores.begin(),
                                                   boundary_scores.end(), 0.0) /
                                   static_cast<double>(boundary_scores.size());
    }
    return result;
}

unsigned parse_color_channel(const std::string& text, const std::string& context) {
    const uint64_t value = parse_coordinate(text, context);
    if (value > 255U) {
        throw std::runtime_error("RGB channel exceeds 255 in " + context);
    }
    return static_cast<unsigned>(value);
}

Rgb parse_rgb(const std::vector<std::string>& fields, const std::string& context) {
    if (fields.size() == 4U) {
        return {parse_color_channel(fields[1], context),
                parse_color_channel(fields[2], context),
                parse_color_channel(fields[3], context)};
    }
    if (fields.size() == 2U) {
        std::vector<std::string> channels;
        size_t begin = 0;
        while (true) {
            const size_t comma = fields[1].find(',', begin);
            channels.emplace_back(fields[1].substr(begin, comma - begin));
            if (comma == std::string::npos) break;
            begin = comma + 1U;
        }
        if (channels.size() == 3U) {
            return {parse_color_channel(channels[0], context),
                    parse_color_channel(channels[1], context),
                    parse_color_channel(channels[2], context)};
        }
    }
    throw std::runtime_error(
        "expected component<TAB>R,G,B or component<TAB>R<TAB>G<TAB>B in " + context);
}

Rgb generated_color(const std::string& name) {
    uint64_t hash = 1469598103934665603ULL;
    for (unsigned char c : name) {
        hash ^= c;
        hash *= 1099511628211ULL;
    }
    hash = mix64(hash);
    return {40U + static_cast<unsigned>(hash & 0xffU) % 176U,
            40U + static_cast<unsigned>((hash >> 8U) & 0xffU) % 176U,
            40U + static_cast<unsigned>((hash >> 16U) & 0xffU) % 176U};
}

std::unordered_map<std::string, Rgb> build_colors(const ReferenceIndex& index,
                                                   const std::string& override_path) {
    std::unordered_map<std::string, Rgb> colors;
    colors.reserve(index.references.size() + 1U);
    for (const auto& reference : index.references) {
        colors.emplace(reference.name, generated_color(reference.name));
    }

    // High-contrast defaults for the core ACRE components. Every other supplied
    // reference receives a stable name-derived color.
    const std::array<std::pair<const char*, Rgb>, 11> builtins = {{
        {"ACRO1", {255, 0, 0}},
        {"COMP-subunit_ACRO_rnd-5_family-38", {0, 0, 255}},
        {"COMP-subunit_ACRO_rnd-5_family-37", {0, 166, 81}},
        {"COMP-subunit_ACRO_rnd-5_family-1624", {255, 127, 0}},
        {"COMP-subunit_ACRO_rnd-5_family-1625", {128, 0, 180}},
        {"MER21A", {166, 86, 40}},
        {"MER21B", {210, 145, 30}},
        {"MER21C", {0, 153, 153}},
        {"MER21-int", {102, 102, 0}},
        {"RMER21A", {216, 27, 96}},
        {"RMER21B", {121, 85, 72}},
    }};
    for (const auto& item : builtins) {
        const auto found = colors.find(item.first);
        if (found != colors.end()) found->second = item.second;
    }
    colors.emplace("unclassified", Rgb{128, 128, 128});

    if (!override_path.empty()) {
        std::ifstream input(override_path);
        if (!input) throw std::runtime_error("cannot open color map: " + override_path);
        std::string line;
        size_t line_number = 0;
        while (std::getline(input, line)) {
            ++line_number;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty() || line.front() == '#') continue;
            const auto fields = split_tabs(line);
            if (fields[0] == "component" || fields[0] == "name") continue;
            const auto found = colors.find(fields[0]);
            if (found == colors.end()) {
                throw std::runtime_error("unknown component in color map " + override_path +
                                         " line " + std::to_string(line_number) + ": " +
                                         fields[0]);
            }
            found->second = parse_rgb(fields, override_path + " line " +
                                      std::to_string(line_number));
        }
    }
    return colors;
}

std::vector<AnnotationSegment> annotate_sequence(const std::string& sequence,
                                                 const ReferenceIndex& index,
                                                 const Options& options) {
    std::vector<AnnotationSegment> segments;
    const size_t ref_count = index.references.size();
    std::vector<double> votes(ref_count);
    std::vector<double> plus(ref_count);
    std::vector<double> minus(ref_count);

    for (size_t begin = 0; begin < sequence.size(); begin += options.window) {
        const size_t end = std::min(sequence.size(), begin + options.window);
        std::fill(votes.begin(), votes.end(), 0.0);
        std::fill(plus.begin(), plus.end(), 0.0);
        std::fill(minus.begin(), minus.end(), 0.0);
        uint64_t valid = 0;
        const std::string_view window(sequence.data() + begin, end - begin);
        for_each_canonical_kmer(window, options.k,
            [&](size_t, uint64_t code, bool query_forward) {
                ++valid;
                const auto found = index.kmers.find(code);
                if (found == index.kmers.end()) return;
                const KmerEntry& entry = found->second;
                const double weight = 1.0 / static_cast<double>(popcount64(entry.members));
                uint64_t members = entry.members;
                while (members) {
                    const unsigned ref = trailing_zeroes64(members);
                    const uint64_t bit = uint64_t{1} << ref;
                    votes[ref] += weight;
                    const bool reference_forward = (entry.reference_forward & bit) != 0;
                    const bool reference_reverse = (entry.reference_reverse & bit) != 0;
                    if (reference_forward != reference_reverse) {
                        const bool same = query_forward == reference_forward;
                        (same ? plus[ref] : minus[ref]) += weight;
                    }
                    members &= members - 1U;
                }
            });

        int best_ref = -1;
        double best = 0;
        double second = 0;
        for (size_t ref = 0; ref < ref_count; ++ref) {
            if (votes[ref] > best) {
                second = best;
                best = votes[ref];
                best_ref = static_cast<int>(ref);
            } else if (votes[ref] > second) {
                second = votes[ref];
            }
        }
        const double density = valid ? best / static_cast<double>(valid) : 0.0;
        if (best < options.min_hits || density < options.min_density) best_ref = -1;
        const double confidence = best > 0 ? std::max(0.0, (best - second) / best) : 0.0;
        char strand = '.';
        if (best_ref >= 0) {
            if (plus[best_ref] > minus[best_ref]) strand = '+';
            else if (minus[best_ref] > plus[best_ref]) strand = '-';
        }

        if (!segments.empty() && segments.back().reference == best_ref &&
            segments.back().strand == strand && segments.back().relative_end == begin) {
            AnnotationSegment& segment = segments.back();
            segment.relative_end = end;
            segment.weighted_hits += best;
            segment.valid_kmers += valid;
            segment.confidence_sum += confidence * static_cast<double>(valid);
        } else {
            segments.push_back({begin, end, best_ref, strand, best, valid,
                                confidence * static_cast<double>(valid)});
        }
    }
    return segments;
}

std::vector<size_t> kmeans(const std::vector<std::vector<float>>& values,
                           size_t cluster_count, uint64_t seed) {
    const size_t n = values.size();
    if (n == 0) return {};
    const size_t dimensions = values.front().size();
    cluster_count = std::max<size_t>(1, std::min(cluster_count, n));

    // Keep standardization parameters and row norms instead of copying the full
    // N x D matrix. At pangenome scale that copy can otherwise consume gigabytes.
    std::vector<double> means(dimensions, 0.0);
    std::vector<double> scales(dimensions, 0.0);
    for (size_t col = 0; col < dimensions; ++col) {
        for (const auto& row : values) means[col] += row[col];
        means[col] /= static_cast<double>(n);
        double variance = 0;
        for (const auto& row : values) {
            const double delta = row[col] - means[col];
            variance += delta * delta;
        }
        scales[col] = n > 1 ? std::sqrt(variance / static_cast<double>(n - 1U)) : 0;
    }
    std::vector<double> norms(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        double norm = 0;
        for (size_t col = 0; col < dimensions; ++col) {
            const double standardized = scales[col] > 1e-12
                ? (values[i][col] - means[col]) / scales[col] : 0.0;
            norm += standardized * standardized;
        }
        norms[i] = std::sqrt(norm);
    }
    auto normalized_value = [&](size_t row, size_t col) {
        if (scales[col] <= 1e-12 || norms[row] <= 0.0) return 0.0;
        return ((values[row][col] - means[col]) / scales[col]) / norms[row];
    };
    auto normalized_row = [&](size_t row) {
        std::vector<double> result(dimensions, 0.0);
        for (size_t col = 0; col < dimensions; ++col) {
            result[col] = normalized_value(row, col);
        }
        return result;
    };
    auto distance_to_center = [&](size_t row, const std::vector<double>& center) {
        double result = 0;
        for (size_t j = 0; j < dimensions; ++j) {
            const double delta = normalized_value(row, j) - center[j];
            result += delta * delta;
        }
        return result;
    };

    std::mt19937_64 engine(seed);
    std::uniform_int_distribution<size_t> choose_first(0, n - 1U);
    std::vector<std::vector<double>> centers;
    centers.reserve(cluster_count);
    centers.push_back(normalized_row(choose_first(engine)));
    std::vector<double> nearest(n, std::numeric_limits<double>::infinity());
    while (centers.size() < cluster_count) {
        double total = 0;
        for (size_t i = 0; i < n; ++i) {
            nearest[i] = std::min(nearest[i], distance_to_center(i, centers.back()));
            total += nearest[i];
        }
        size_t selected = 0;
        if (total > 0) {
            std::uniform_real_distribution<double> draw(0, total);
            double target = draw(engine);
            for (; selected + 1U < n; ++selected) {
                target -= nearest[selected];
                if (target <= 0) break;
            }
        } else {
            selected = centers.size() % n;
        }
        centers.push_back(normalized_row(selected));
    }

    std::vector<size_t> assignments(n, cluster_count);
    for (unsigned iteration = 0; iteration < 100; ++iteration) {
        bool changed = false;
        for (size_t i = 0; i < n; ++i) {
            double best_distance = std::numeric_limits<double>::infinity();
            size_t best_cluster = 0;
            for (size_t cluster = 0; cluster < cluster_count; ++cluster) {
                const double d = distance_to_center(i, centers[cluster]);
                if (d < best_distance) {
                    best_distance = d;
                    best_cluster = cluster;
                }
            }
            if (assignments[i] != best_cluster) changed = true;
            assignments[i] = best_cluster;
        }
        std::vector<std::vector<double>> sums(cluster_count,
                                              std::vector<double>(dimensions, 0));
        std::vector<size_t> counts(cluster_count, 0);
        for (size_t i = 0; i < n; ++i) {
            ++counts[assignments[i]];
            for (size_t j = 0; j < dimensions; ++j) {
                sums[assignments[i]][j] += normalized_value(i, j);
            }
        }
        bool had_empty_cluster = false;
        for (size_t cluster = 0; cluster < cluster_count; ++cluster) {
            if (counts[cluster] == 0) {
                had_empty_cluster = true;
                size_t farthest = 0;
                double farthest_distance = -1;
                for (size_t i = 0; i < n; ++i) {
                    const double d = distance_to_center(i, centers[assignments[i]]);
                    if (d > farthest_distance) { farthest_distance = d; farthest = i; }
                }
                centers[cluster] = normalized_row(farthest);
            } else {
                for (size_t j = 0; j < dimensions; ++j) {
                    centers[cluster][j] = sums[cluster][j] /
                                          static_cast<double>(counts[cluster]);
                }
            }
        }
        if (!changed && !had_empty_cluster) break;
    }

    // Stabilize labels: cluster 1 is the one whose first input interval appears first.
    std::vector<size_t> first(cluster_count, n);
    for (size_t i = 0; i < n; ++i) first[assignments[i]] = std::min(first[assignments[i]], i);
    std::vector<size_t> order(cluster_count);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return first[a] < first[b];
    });
    std::vector<size_t> relabel(cluster_count);
    for (size_t i = 0; i < cluster_count; ++i) relabel[order[i]] = i;
    for (size_t& assignment : assignments) assignment = relabel[assignment];
    return assignments;
}

double median_value(std::vector<double> values) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    const size_t middle = values.size() / 2U;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if (values.size() % 2U != 0U) return upper;
    const double lower = *std::max_element(values.begin(), values.begin() + middle);
    return (lower + upper) / 2.0;
}

double unit_sketch_ani(const KmerSketch& first, const KmerSketch& second, unsigned k) {
    return mash_ani(sketch_jaccard(first, second), k);
}

SwitchThresholds switch_thresholds(const Options& options) {
    SwitchThresholds thresholds;
    thresholds.flank_units = options.switch_flank;
    thresholds.sketch_size = options.sketch_size;
    thresholds.min_coherence = options.switch_min_coherence;
    thresholds.max_cross_ani = options.switch_max_cross_ani;
    thresholds.min_ani_drop = options.switch_min_drop;
    thresholds.max_adjacent_ani = options.switch_max_adjacent_ani;
    return thresholds;
}

DeepBoundaryRecord make_boundary_record(const BedRecord& array,
                                        const ArraySegmentation& segmentation,
                                        const SwitchBoundary& boundary) {
    if (boundary.boundary_index == 0U ||
        boundary.boundary_index >= segmentation.units.size()) {
        throw std::runtime_error("internal switch boundary index is out of range");
    }
    DeepBoundaryRecord record;
    record.array_id = array.id;
    record.source_name = array.source_name;
    record.chrom = array.chrom;
    record.array_start = array.start;
    record.array_end = array.end;
    record.boundary_after_unit = boundary.boundary_index;
    record.boundary_before_unit = boundary.boundary_index + 1U;
    record.boundary_position =
        array.start + segmentation.units[boundary.boundary_index].first;
    record.metrics = boundary;
    return record;
}

const char* deep_metric_header() {
    return "array_id\tsource_name\tchrom\tboundary_after_unit\tboundary_before_unit"
           "\tboundary_position\tadjacent_phase_sketch_identity"
           "\tleft_flank_phase_sketch_identity\tright_flank_phase_sketch_identity"
           "\tminimum_flank_phase_sketch_identity\tcross_flank_phase_sketch_identity"
           "\tphase_sketch_identity_drop\tpasses_thresholds";
}

void write_boundary_record(std::ostream& output, const DeepBoundaryRecord& record) {
    const auto& boundary = record.metrics;
    output << record.array_id << '\t' << record.source_name << '\t' << record.chrom << '\t'
           << record.boundary_after_unit << '\t' << record.boundary_before_unit << '\t'
           << record.boundary_position << '\t' << boundary.adjacent_ani << '\t'
           << boundary.left_flank_ani << '\t' << boundary.right_flank_ani << '\t'
           << boundary.minimum_flank_coherence << '\t' << boundary.cross_flank_ani << '\t'
           << boundary.ani_drop << '\t' << (boundary.candidate ? "True" : "False");
}

int switch_bed_score(double identity_drop) {
    return static_cast<int>(std::llround(
        std::min(1.0, std::max(0.0, identity_drop)) * 1000.0));
}

void write_switch_coordinate_beds(const Options& options,
                                  const std::vector<DeepSwitchEvent>& events) {
    const std::string switches_bed_path = options.output_prefix + ".switches.bed";
    const std::string arrays_bed_path = options.output_prefix + ".switch_arrays.bed";
    std::ofstream switches_bed(switches_bed_path);
    std::ofstream arrays_bed(arrays_bed_path);
    if (!switches_bed || !arrays_bed) {
        throw std::runtime_error("cannot create switch BED outputs with prefix: " +
                                 options.output_prefix);
    }
    switches_bed << "#chrom\tstart\tend\tname\tscore\tstrand\tarray_id"
                    "\tboundary_after_unit\tboundary_before_unit"
                    "\tphase_sketch_identity_drop\tsource_name\n"
                 << std::setprecision(9);
    arrays_bed << "#chrom\tstart\tend\tname\tscore\tstrand\tswitch_count"
                  "\tswitch_ids\tmax_phase_sketch_identity_drop\tsource_name\n"
               << std::setprecision(9);

    for (const auto& event : events) {
        const auto& record = event.record;
        const uint64_t marker_end = record.boundary_position + 1U;
        switches_bed << record.chrom << '\t' << record.boundary_position << '\t'
                     << marker_end << '\t' << event.id << '\t'
                     << switch_bed_score(record.metrics.ani_drop) << "\t.\t"
                     << record.array_id << '\t' << record.boundary_after_unit << '\t'
                     << record.boundary_before_unit << '\t' << record.metrics.ani_drop
                     << '\t' << record.source_name << '\n';
    }

    size_t begin = 0;
    while (begin < events.size()) {
        size_t end = begin + 1U;
        double maximum_drop = events[begin].record.metrics.ani_drop;
        while (end < events.size() &&
               events[end].record.array_id == events[begin].record.array_id) {
            maximum_drop = std::max(maximum_drop, events[end].record.metrics.ani_drop);
            ++end;
        }
        const auto& record = events[begin].record;
        arrays_bed << record.chrom << '\t' << record.array_start << '\t'
                   << record.array_end << '\t' << record.array_id << '\t'
                   << switch_bed_score(maximum_drop) << "\t.\t" << end - begin << '\t';
        for (size_t index = begin; index < end; ++index) {
            if (index != begin) arrays_bed << ',';
            arrays_bed << events[index].id;
        }
        arrays_bed << '\t' << maximum_drop << '\t' << record.source_name << '\n';
        begin = end;
    }
}

std::pair<size_t, size_t> event_context(const DeepSwitchEvent& event, size_t radius) {
    const size_t left_available = event.boundary_global - event.array_begin;
    const size_t right_available = event.array_end - event.boundary_global;
    const size_t begin = event.boundary_global - std::min(radius, left_available);
    const size_t end = event.boundary_global + std::min(radius, right_available);
    return {begin, end};
}

void write_switch_context_subfamilies(const Options& options,
                                      const std::vector<AcreUnit>& units,
                                      const std::vector<KmerSketch>& sketches,
                                      const std::vector<DeepSwitchEvent>& events,
                                      std::ostream& bed_output,
                                      std::ostream& summary_output) {
    if (units.size() != sketches.size()) {
        throw std::runtime_error("internal deep-subfamily sketch count mismatch");
    }

    // Collect only the fixed-radius context around called switches. Global
    // subfamilies are discovered across this sparse subset, not all pangenome
    // units, and candidate lookup uses an inverted sketch-hash index.
    std::vector<size_t> selected;
    for (const auto& event : events) {
        const auto range = event_context(event, options.switch_radius);
        for (size_t index = range.first; index < range.second; ++index) {
            selected.push_back(index);
        }
    }
    std::sort(selected.begin(), selected.end(), [&](size_t left, size_t right) {
        const auto left_key = std::tie(units[left].chrom, units[left].start, units[left].end,
                                       units[left].array_id, units[left].unit_number);
        const auto right_key = std::tie(units[right].chrom, units[right].start, units[right].end,
                                        units[right].array_id, units[right].unit_number);
        return left_key < right_key;
    });
    selected.erase(std::unique(selected.begin(), selected.end()), selected.end());

    const size_t unassigned = std::numeric_limits<size_t>::max();
    std::vector<size_t> subfamily_by_unit(units.size(), unassigned);
    std::vector<double> representative_ani(units.size(), 0.0);
    std::vector<size_t> representatives;
    std::unordered_map<uint64_t, std::vector<size_t>> postings;
    postings.reserve(selected.size() * 8U + 1U);
    for (size_t unit_index : selected) {
        std::unordered_map<size_t, size_t> overlap_counts;
        for (uint64_t hash : sketches[unit_index].hashes) {
            const auto found = postings.find(hash);
            if (found == postings.end()) continue;
            for (size_t cluster : found->second) ++overlap_counts[cluster];
        }
        for (uint64_t hash : sketches[unit_index].reverse_hashes) {
            const auto found = postings.find(hash);
            if (found == postings.end()) continue;
            for (size_t cluster : found->second) ++overlap_counts[cluster];
        }
        std::vector<std::pair<size_t, size_t>> candidates(overlap_counts.begin(),
                                                          overlap_counts.end());
        std::sort(candidates.begin(), candidates.end(), [](const auto& left, const auto& right) {
            if (left.second != right.second) return left.second > right.second;
            return left.first < right.first;
        });
        size_t best_cluster = unassigned;
        double best_ani = -1.0;
        for (const auto& candidate : candidates) {
            const size_t cluster = candidate.first;
            const double ani = unit_sketch_ani(sketches[unit_index],
                                                sketches[representatives[cluster]], options.k);
            if (ani >= options.subfamily_min_ani &&
                (ani > best_ani + 1e-12 ||
                 (std::abs(ani - best_ani) <= 1e-12 && cluster < best_cluster))) {
                best_cluster = cluster;
                best_ani = ani;
            }
        }
        if (best_cluster == unassigned) {
            best_cluster = representatives.size();
            representatives.push_back(unit_index);
            for (uint64_t hash : sketches[unit_index].hashes) {
                postings[hash].push_back(best_cluster);
            }
            for (uint64_t hash : sketches[unit_index].reverse_hashes) {
                postings[hash].push_back(best_cluster);
            }
            best_ani = 1.0;
        }
        subfamily_by_unit[unit_index] = best_cluster;
        representative_ani[unit_index] = best_ani;
    }

    bed_output << "#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd"
                  "\titemRgb\tacre_id\tarray_id\tunit_in_array\tswitch_id"
                  "\tunits_from_boundary\tswitch_state_local\tswitch_subfamily_global"
                  "\tphase_sketch_identity_to_left\tphase_sketch_identity_to_right"
                  "\tright_minus_left_phase_sketch_identity"
                  "\tphase_sketch_identity_to_representative"
                  "\tvalid_kmers\tdistinct_sketch_tokens"
                  "\tsource_name\n" << std::setprecision(9);
    std::vector<size_t> membership_counts(representatives.size(), 0U);
    std::vector<std::unordered_set<std::string>> cluster_arrays(representatives.size());
    std::vector<std::unordered_set<std::string>> cluster_events(representatives.size());
    for (const auto& event : events) {
        const auto range = event_context(event, options.switch_radius);
        for (size_t index = range.first; index < range.second; ++index) {
            std::vector<double> left_values;
            std::vector<double> right_values;
            for (size_t flank_index = event.boundary_global - options.switch_flank;
                 flank_index < event.boundary_global; ++flank_index) {
                if (flank_index != index) {
                    left_values.push_back(unit_sketch_ani(sketches[index],
                                                          sketches[flank_index], options.k));
                }
            }
            for (size_t flank_index = event.boundary_global;
                 flank_index < event.boundary_global + options.switch_flank; ++flank_index) {
                if (flank_index != index) {
                    right_values.push_back(unit_sketch_ani(sketches[index],
                                                           sketches[flank_index], options.k));
                }
            }
            const double left_ani = median_value(std::move(left_values));
            const double right_ani = median_value(std::move(right_values));
            const double margin = right_ani - left_ani;
            const std::string local_state = margin > 0.01 ? "right" :
                                            margin < -0.01 ? "left" : "ambiguous";
            const size_t cluster = subfamily_by_unit[index];
            if (cluster == unassigned) {
                throw std::runtime_error("internal unassigned switch-region subfamily");
            }
            const std::string label = "switch_subfamily_" + std::to_string(cluster + 1U);
            const Rgb color = generated_color(label);
            const int score = static_cast<int>(std::llround(
                std::min(1.0, std::abs(margin) * 10.0) * 1000.0));
            const int64_t offset = static_cast<int64_t>(index) -
                                   static_cast<int64_t>(event.boundary_global);
            bed_output << units[index].chrom << '\t' << units[index].start << '\t'
                       << units[index].end << '\t' << label << '\t' << score
                       << "\t.\t" << units[index].start << '\t' << units[index].end << '\t'
                       << color.text() << '\t' << units[index].id << '\t'
                       << units[index].array_id << '\t' << units[index].unit_number << '\t'
                       << event.id << '\t' << offset << '\t' << local_state << '\t'
                       << cluster + 1U << '\t' << left_ani << '\t' << right_ani << '\t'
                       << margin << '\t' << representative_ani[index] << '\t'
                       << sketches[index].valid_kmers << '\t' << sketches[index].distinct_kmers
                       << '\t' << units[index].source_name << '\n';
            ++membership_counts[cluster];
            cluster_arrays[cluster].insert(units[index].array_id);
            cluster_events[cluster].insert(event.id);
        }
    }

    summary_output << "switch_subfamily_global\trepresentative_acre_id"
                      "\trepresentative_array_id\tunique_units\tswitch_region_memberships"
                      "\tarrays\tswitch_events\n";
    std::vector<size_t> unique_counts(representatives.size(), 0U);
    for (size_t index : selected) ++unique_counts[subfamily_by_unit[index]];
    for (size_t cluster = 0; cluster < representatives.size(); ++cluster) {
        const size_t representative = representatives[cluster];
        summary_output << cluster + 1U << '\t' << units[representative].id << '\t'
                       << units[representative].array_id << '\t' << unique_counts[cluster]
                       << '\t' << membership_counts[cluster] << '\t'
                       << cluster_arrays[cluster].size() << '\t'
                       << cluster_events[cluster].size() << '\n';
    }
    std::cerr << "Deep subfamily analysis called " << events.size()
              << " nonredundant switch event(s) and assigned " << selected.size()
              << " switch-region ACRE(s) to " << representatives.size()
              << " global subfamily/subfamilies\n";
}

void write_deep_subfamily_analysis(const Options& options,
                                   const std::vector<AcreUnit>& units,
                                   const std::vector<KmerSketch>& sketches) {
    if (units.size() != sketches.size()) {
        throw std::runtime_error("internal deep-subfamily sketch count mismatch");
    }

    const std::string boundary_path = options.output_prefix + ".switch_boundaries.tsv";
    const std::string switch_path = options.output_prefix + ".switches.tsv";
    const std::string bed_path = options.output_prefix + ".switch_subfamilies.bed";
    const std::string summary_path = options.output_prefix + ".switch_subfamily_summary.tsv";
    std::ofstream boundary_output(boundary_path);
    std::ofstream switch_output(switch_path);
    std::ofstream bed_output(bed_path);
    std::ofstream summary_output(summary_path);
    if (!boundary_output || !switch_output || !bed_output || !summary_output) {
        throw std::runtime_error("cannot create deep-subfamily outputs with prefix: " +
                                 options.output_prefix);
    }

    const auto thresholds = switch_thresholds(options);
    boundary_output << deep_metric_header() << '\n' << std::setprecision(9);
    switch_output << "switch_id\t" << deep_metric_header() << '\n'
                  << std::setprecision(9);

    std::vector<DeepSwitchEvent> events;
    size_t array_begin = 0;
    size_t switch_number = 0;
    while (array_begin < units.size()) {
        size_t array_end = array_begin + 1U;
        while (array_end < units.size() &&
               units[array_end].array_id == units[array_begin].array_id) {
            ++array_end;
        }
        std::vector<KmerSketch> array_sketches(
            sketches.begin() + static_cast<std::ptrdiff_t>(array_begin),
            sketches.begin() + static_cast<std::ptrdiff_t>(array_end));
        const auto boundaries = detect_switch_boundaries(array_sketches, options.k,
                                                          thresholds);
        for (const auto& boundary : boundaries) {
            DeepBoundaryRecord record;
            record.array_id = units[array_begin].array_id;
            record.source_name = units[array_begin].source_name;
            record.chrom = units[array_begin].chrom;
            record.array_start = units[array_begin].array_start;
            record.array_end = units[array_begin].array_end;
            const size_t global = array_begin + boundary.boundary_index;
            record.boundary_after_unit = units[global - 1U].unit_number;
            record.boundary_before_unit = units[global].unit_number;
            record.boundary_position = units[global].start;
            record.metrics = boundary;
            write_boundary_record(boundary_output, record);
            boundary_output << '\n';
        }
        const auto collapsed = collapse_switch_boundaries(boundaries, thresholds);
        for (const auto& boundary : collapsed) {
            std::ostringstream id;
            id << "SWITCH_" << std::setw(6) << std::setfill('0') << ++switch_number;
            DeepSwitchEvent event;
            event.id = id.str();
            event.array_begin = array_begin;
            event.array_end = array_end;
            event.boundary_global = array_begin + boundary.boundary_index;
            event.record.array_id = units[array_begin].array_id;
            event.record.source_name = units[array_begin].source_name;
            event.record.chrom = units[array_begin].chrom;
            event.record.array_start = units[array_begin].array_start;
            event.record.array_end = units[array_begin].array_end;
            event.record.boundary_after_unit = units[event.boundary_global - 1U].unit_number;
            event.record.boundary_before_unit = units[event.boundary_global].unit_number;
            event.record.boundary_position = units[event.boundary_global].start;
            event.record.metrics = boundary;
            events.push_back(event);
            switch_output << event.id << '\t';
            write_boundary_record(switch_output, event.record);
            switch_output << '\n';
        }
        array_begin = array_end;
    }
    write_switch_coordinate_beds(options, events);
    write_switch_context_subfamilies(options, units, sketches, events,
                                     bed_output, summary_output);
}

std::string seconds_since(Clock::time_point start) {
    const double seconds = std::chrono::duration<double>(Clock::now() - start).count();
    std::ostringstream out;
    out << std::fixed << std::setprecision(3) << seconds << 's';
    return out.str();
}

std::string numbered_id(const char* prefix, size_t number) {
    std::ostringstream id;
    id << prefix << std::setw(6) << std::setfill('0') << number;
    return id.str();
}

int run_streaming_pangenome_analysis(const Options& options,
                                     const std::vector<BedRecord>& arrays,
                                     const std::shared_ptr<GenomeIndex>& genome,
                                     Clock::time_point overall_start) {
    const auto analysis_start = Clock::now();
    const auto thresholds = switch_thresholds(options);

    const std::string boundary_path = options.output_prefix + ".switch_boundaries.tsv";
    const std::string switch_path = options.output_prefix + ".switches.tsv";
    const std::string bed_path = options.output_prefix + ".switch_subfamilies.bed";
    const std::string summary_path = options.output_prefix + ".switch_subfamily_summary.tsv";
    const std::string acres_path = options.output_prefix + ".acres.bed";
    const std::string partials_path = options.output_prefix + ".partials.bed";
    const std::string arrays_path = options.output_prefix + ".arrays.tsv";
    std::ofstream boundary_output(boundary_path);
    std::ofstream switch_output(switch_path);
    std::ofstream bed_output(bed_path);
    std::ofstream summary_output(summary_path);
    std::ofstream acres_output(acres_path);
    std::ofstream partials_output(partials_path);
    std::ofstream arrays_output(arrays_path);
    if (!boundary_output || !switch_output || !bed_output || !summary_output ||
        !acres_output || !partials_output || !arrays_output) {
        throw std::runtime_error("cannot create pangenome outputs with prefix: " +
                                 options.output_prefix);
    }
    boundary_output << deep_metric_header() << '\n' << std::setprecision(9);
    switch_output << "switch_id\t" << deep_metric_header() << '\n'
                  << std::setprecision(9);
    arrays_output << "array_id\tchrom\tstart\tend\tinput_length\tcomplete_acres\t"
                     "complete_bp\tresidual_bp\tsegmentation_method\testimated_period\t"
                     "periodicity_score\tsource_name\n" << std::setprecision(9);
    acres_output << std::setprecision(9);

    size_t thread_count = options.threads;
    if (thread_count == 0U) {
        thread_count = std::thread::hardware_concurrency();
        if (thread_count == 0U) thread_count = 1U;
        thread_count = std::min<size_t>(thread_count, 8U);
    }
    thread_count = std::max<size_t>(1U, std::min(thread_count, arrays.size()));
    const size_t batch_size = thread_count > std::numeric_limits<size_t>::max() / 2U
        ? thread_count : thread_count * 2U;

    std::vector<AcreUnit> context_units;
    std::vector<KmerSketch> context_sketches;
    std::vector<DeepSwitchEvent> events;
    size_t acre_count = 0;
    size_t partial_count = 0;
    size_t switch_count = 0;
    uint64_t input_bases = 0;
    uint64_t complete_bases = 0;
    uint64_t residual_bases = 0;
    uint64_t retained_hashes = 0;

    for (size_t batch_begin = 0; batch_begin < arrays.size(); batch_begin += batch_size) {
        const size_t batch_end = std::min(arrays.size(), batch_begin + batch_size);
        const size_t batch_count = batch_end - batch_begin;
        std::vector<StreamingArrayAnalysis> results(batch_count);
        std::atomic<size_t> next{0U};
        std::atomic<bool> failed{false};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        std::vector<std::thread> workers;
        workers.reserve(std::min(thread_count, batch_count));

        for (size_t worker = 0; worker < std::min(thread_count, batch_count); ++worker) {
            workers.emplace_back([&] {
                try {
                    FastaReader reader(genome);
                    while (!failed.load(std::memory_order_relaxed)) {
                        const size_t offset = next.fetch_add(1U, std::memory_order_relaxed);
                        if (offset >= batch_count) break;
                        const size_t array_index = batch_begin + offset;
                        const BedRecord& array = arrays[array_index];
                        StreamingArrayAnalysis result;
                        std::string sequence = reader.fetch(array.chrom, array.start, array.end);
                        result.input_bases = sequence.size();
                        result.segmentation = segment_sequence(sequence, options);
                        result.complete_units = result.segmentation.units.size();

                        std::vector<KmerSketch> sketches;
                        sketches.reserve(result.complete_units);
                        for (const auto& interval : result.segmentation.units) {
                            const size_t length = interval.second - interval.first;
                            result.complete_bases += length;
                            const std::string_view view(sequence.data() + interval.first, length);
                            sketches.push_back(make_bottom_sketch(
                                view, options.k, options.sketch_size,
                                options.sketch_position_bins));
                        }
                        for (const auto& interval : result.segmentation.residuals) {
                            result.residual_bases += interval.second - interval.first;
                        }
                        if (result.complete_bases + result.residual_bases !=
                            result.input_bases) {
                            throw std::runtime_error(
                                "internal streaming ACRE segmentation coverage failure");
                        }

                        const auto boundaries = detect_switch_boundaries(
                            sketches, options.k, thresholds);
                        result.boundaries.reserve(boundaries.size());
                        for (const auto& boundary : boundaries) {
                            result.boundaries.push_back(make_boundary_record(
                                array, result.segmentation, boundary));
                        }
                        result.switches = collapse_switch_boundaries(boundaries, thresholds);

                        std::vector<bool> retain(result.complete_units, false);
                        const size_t radius = std::max(options.switch_radius,
                                                       options.switch_flank);
                        for (const auto& call : result.switches) {
                            const size_t begin = call.boundary_index -
                                std::min(radius, call.boundary_index);
                            const size_t available = result.complete_units -
                                                     call.boundary_index;
                            const size_t end = call.boundary_index +
                                               std::min(radius, available);
                            for (size_t index = begin; index < end; ++index) {
                                retain[index] = true;
                            }
                        }
                        for (size_t index = 0; index < result.complete_units; ++index) {
                            if (!retain[index]) continue;
                            const auto interval = result.segmentation.units[index];
                            AcreUnit unit;
                            unit.chrom = array.chrom;
                            unit.start = array.start + interval.first;
                            unit.end = array.start + interval.second;
                            unit.array_start = array.start;
                            unit.array_end = array.end;
                            unit.source_name = array.source_name;
                            unit.array_id = array.id;
                            unit.unit_number = index + 1U;
                            unit.estimated_period = result.segmentation.estimated_period;
                            unit.periodicity_score = result.segmentation.periodicity_score;
                            result.retained_local_indices.push_back(index);
                            result.retained_units.push_back(std::move(unit));
                            result.retained_hashes += sketches[index].hashes.size() +
                                                      sketches[index].reverse_hashes.size();
                            result.retained_sketches.push_back(std::move(sketches[index]));
                        }
                        results[offset] = std::move(result);
                    }
                } catch (...) {
                    failed.store(true, std::memory_order_relaxed);
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (!worker_error) worker_error = std::current_exception();
                }
            });
        }
        for (auto& worker : workers) worker.join();
        if (worker_error) std::rethrow_exception(worker_error);

        for (size_t offset = 0; offset < batch_count; ++offset) {
            const size_t array_index = batch_begin + offset;
            const BedRecord& array = arrays[array_index];
            auto& result = results[offset];
            input_bases += result.input_bases;
            complete_bases += result.complete_bases;
            residual_bases += result.residual_bases;
            retained_hashes += result.retained_hashes;

            for (const auto& record : result.boundaries) {
                write_boundary_record(boundary_output, record);
                boundary_output << '\n';
            }

            const size_t acre_base = acre_count;
            for (size_t index = 0; index < result.segmentation.units.size(); ++index) {
                const auto interval = result.segmentation.units[index];
                const std::string id = numbered_id("ACRE_", acre_base + index + 1U);
                acres_output << array.chrom << '\t' << array.start + interval.first << '\t'
                             << array.start + interval.second << '\t' << id
                             << "\t0\t.\t" << array.id << '\t' << index + 1U << '\t'
                             << result.segmentation.estimated_period << '\t'
                             << result.segmentation.periodicity_score << '\t'
                             << array.source_name << '\n';
            }
            acre_count += result.complete_units;
            for (const auto& interval : result.segmentation.residuals) {
                partials_output << array.chrom << '\t' << array.start + interval.first << '\t'
                                << array.start + interval.second << '\t'
                                << numbered_id("PARTIAL_", ++partial_count)
                                << "\t0\t.\t" << array.id << '\t'
                                << (result.segmentation.units.empty()
                                        ? "shorter_than_min_period" : "edge_residual")
                                << '\n';
            }
            arrays_output << array.id << '\t' << array.chrom << '\t' << array.start << '\t'
                          << array.end << '\t' << result.input_bases << '\t'
                          << result.complete_units << '\t' << result.complete_bases << '\t'
                          << result.residual_bases << '\t' << result.segmentation.method << '\t'
                          << result.segmentation.estimated_period << '\t'
                          << result.segmentation.periodicity_score << '\t'
                          << array.source_name << '\n';

            const size_t context_base = context_units.size();
            for (size_t index = 0; index < result.retained_units.size(); ++index) {
                auto unit = std::move(result.retained_units[index]);
                unit.id = numbered_id("ACRE_", acre_base + unit.unit_number);
                context_units.push_back(std::move(unit));
                context_sketches.push_back(std::move(result.retained_sketches[index]));
            }
            auto compact_index = [&](size_t local_index) -> size_t {
                const auto found = std::lower_bound(result.retained_local_indices.begin(),
                                                    result.retained_local_indices.end(),
                                                    local_index);
                if (found == result.retained_local_indices.end() || *found != local_index) {
                    throw std::runtime_error(
                        "internal streaming switch context is incomplete");
                }
                return context_base + static_cast<size_t>(
                    found - result.retained_local_indices.begin());
            };
            for (const auto& call : result.switches) {
                DeepSwitchEvent event;
                event.id = numbered_id("SWITCH_", ++switch_count);
                event.record = make_boundary_record(array, result.segmentation, call);
                const size_t report_begin = call.boundary_index -
                    std::min(options.switch_radius, call.boundary_index);
                const size_t available = result.complete_units - call.boundary_index;
                const size_t report_end = call.boundary_index +
                    std::min(options.switch_radius, available);
                event.array_begin = compact_index(report_begin);
                event.array_end = compact_index(report_end - 1U) + 1U;
                event.boundary_global = compact_index(call.boundary_index);
                switch_output << event.id << '\t';
                write_boundary_record(switch_output, event.record);
                switch_output << '\n';
                events.push_back(std::move(event));
            }
        }
    }

    if (acre_count == 0U) {
        throw std::runtime_error("no complete ACREs fit the configured period range");
    }
    write_switch_coordinate_beds(options, events);
    write_switch_context_subfamilies(options, context_units, context_sketches,
                                     events, bed_output, summary_output);
    std::cerr << "Streaming pangenome mode processed " << arrays.size() << " arrays ("
              << input_bases << " bp), emitted " << acre_count << " complete ACREs ("
              << complete_bases << " bp) and " << partial_count << " residuals ("
              << residual_bases << " bp), and retained " << context_units.size()
              << " switch-context ACRE sketches / " << retained_hashes << " hashes ("
              << retained_hashes * sizeof(uint64_t) << " bytes)\n"
              << "Finished N-TRIAGE pangenome analysis in "
              << seconds_since(overall_start) << " (streaming analysis "
              << seconds_since(analysis_start) << ")\n";
    return 0;
}

}  // namespace

std::vector<BedRecord> read_bed(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open BED: " + path);
    std::vector<BedRecord> records;
    std::unordered_set<std::string> used_ids;
    std::string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line.front() == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 3) {
            throw std::runtime_error("BED line " + std::to_string(line_number) +
                                     " has fewer than 3 columns");
        }
        BedRecord record;
        record.chrom = fields[0];
        const std::string context = path + " line " + std::to_string(line_number);
        record.start = parse_coordinate(fields[1], context);
        record.end = parse_coordinate(fields[2], context);
        if (record.end <= record.start) {
            throw std::runtime_error("empty or reversed interval in " + context);
        }
        record.source_name = fields.size() >= 4 ? fields[3] : ".";
        const bool stable_array_id = record.source_name.size() > 6U &&
            record.source_name.compare(0, 6, "ARRAY_") == 0 &&
            std::all_of(record.source_name.begin() + 6, record.source_name.end(),
                        [](unsigned char c) { return std::isdigit(c) != 0; });
        if (stable_array_id) {
            record.id = record.source_name;
        } else {
            std::ostringstream id;
            id << "ARRAY_" << std::setw(6) << std::setfill('0') << (records.size() + 1U);
            record.id = id.str();
        }
        if (!used_ids.insert(record.id).second) {
            throw std::runtime_error("duplicate resolved array ID in " + context + ": " +
                                     record.id);
        }
        records.push_back(std::move(record));
    }
    if (records.empty()) throw std::runtime_error("BED has no intervals: " + path);
    return records;
}

std::vector<Reference> read_fasta(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open reference FASTA: " + path);
    std::vector<Reference> references;
    std::unordered_set<std::string> names;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line.front() == '>') {
            const size_t end = line.find_first_of(" \t", 1);
            const std::string name = line.substr(1, end == std::string::npos ? end : end - 1U);
            if (name.empty()) throw std::runtime_error("empty FASTA reference name");
            if (!names.insert(name).second) throw std::runtime_error("duplicate reference: " + name);
            references.push_back({name, {}});
        } else if (!line.empty()) {
            if (references.empty()) throw std::runtime_error("sequence before FASTA header");
            for (char c : line) if (c != ' ' && c != '\t') references.back().sequence.push_back(c);
        }
    }
    if (references.empty()) throw std::runtime_error("reference FASTA is empty: " + path);
    for (const auto& reference : references) {
        if (reference.sequence.empty()) throw std::runtime_error("empty reference: " + reference.name);
    }
    return references;
}

ReferenceIndex build_reference_index(std::vector<Reference> references, unsigned k) {
    if (references.size() > 64U) {
        throw std::runtime_error("at most 64 component references are supported");
    }
    ReferenceIndex index;
    index.references = std::move(references);
    size_t total_bases = 0;
    for (const auto& reference : index.references) total_bases += reference.sequence.size();
    index.kmers.reserve(total_bases);
    for (size_t ref = 0; ref < index.references.size(); ++ref) {
        const uint64_t bit = uint64_t{1} << ref;
        for_each_canonical_kmer(index.references[ref].sequence, k,
            [&](size_t, uint64_t code, bool forward) {
                KmerEntry& entry = index.kmers[code];
                entry.members |= bit;
                (forward ? entry.reference_forward : entry.reference_reverse) |= bit;
            });
    }
    index.diagnostic_kmers.assign(index.references.size(), 0);
    for (const auto& item : index.kmers) {
        if (popcount64(item.second.members) == 1U) {
            ++index.diagnostic_kmers[trailing_zeroes64(item.second.members)];
        }
    }
    return index;
}

int run_pipeline(const Options& options) {
    const auto overall_start = Clock::now();
    auto arrays = read_bed(options.bed_path);
    auto genome = std::make_shared<GenomeIndex>(options.fasta_path);
    if (options.pangenome) {
        std::cerr << "Loaded " << arrays.size()
                  << " ACRE arrays for streaming pangenome switch analysis\n";
        return run_streaming_pangenome_analysis(options, arrays, genome, overall_start);
    }
    auto reference_index = build_reference_index(read_fasta(options.refs_path), options.k);
    const auto colors = build_colors(reference_index, options.colors_path);
    std::cerr << "Loaded " << arrays.size() << " ACRE arrays, "
              << reference_index.references.size() << " references, and "
              << reference_index.kmers.size() << " distinct reference k-mers\n";

    size_t thread_count = options.threads;
    if (thread_count == 0) {
        thread_count = std::thread::hardware_concurrency();
        if (thread_count == 0) thread_count = 1;
        thread_count = std::min<size_t>(thread_count, 8);
    }
    thread_count = std::max<size_t>(1, std::min(thread_count, arrays.size()));
    std::vector<std::string> array_sequences(arrays.size());
    std::atomic<size_t> next{0};
    std::atomic<bool> failed{false};
    std::exception_ptr worker_error;
    std::mutex error_mutex;
    const auto extraction_start = Clock::now();
    std::vector<std::thread> workers;
    workers.reserve(thread_count);
    for (size_t worker = 0; worker < thread_count; ++worker) {
        workers.emplace_back([&] {
            try {
                FastaReader reader(genome);
                while (!failed.load(std::memory_order_relaxed)) {
                    const size_t i = next.fetch_add(1, std::memory_order_relaxed);
                    if (i >= arrays.size()) break;
                    array_sequences[i] = reader.fetch(arrays[i].chrom, arrays[i].start,
                                                      arrays[i].end);
                }
            } catch (...) {
                failed.store(true, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
        });
    }
    for (auto& worker : workers) worker.join();
    if (worker_error) std::rethrow_exception(worker_error);
    uint64_t input_bases = 0;
    for (const auto& sequence : array_sequences) input_bases += sequence.size();
    std::cerr << "Fetched " << input_bases << " bp with " << thread_count
              << " thread(s) in " << seconds_since(extraction_start) << '\n';

    const auto segmentation_start = Clock::now();
    std::vector<ArraySegmentation> segmentations(arrays.size());
    std::vector<AcreUnit> units;
    std::vector<ResidualFragment> residuals;
    units.reserve(static_cast<size_t>(input_bases / options.period_target) + arrays.size());
    size_t acre_number = 0;
    uint64_t complete_bases = 0;
    uint64_t residual_bases = 0;
    for (size_t array_index = 0; array_index < arrays.size(); ++array_index) {
        const BedRecord& array = arrays[array_index];
        ArraySegmentation segmentation = segment_sequence(array_sequences[array_index], options);
        for (size_t unit_index = 0; unit_index < segmentation.units.size(); ++unit_index) {
            const auto interval = segmentation.units[unit_index];
            AcreUnit unit;
            unit.chrom = array.chrom;
            unit.start = array.start + interval.first;
            unit.end = array.start + interval.second;
            unit.array_start = array.start;
            unit.array_end = array.end;
            unit.source_name = array.source_name;
            unit.array_id = array.id;
            unit.unit_number = unit_index + 1U;
            unit.estimated_period = segmentation.estimated_period;
            unit.periodicity_score = segmentation.periodicity_score;
            std::ostringstream id;
            id << "ACRE_" << std::setw(6) << std::setfill('0') << ++acre_number;
            unit.id = id.str();
            unit.sequence = array_sequences[array_index].substr(
                interval.first, interval.second - interval.first);
            complete_bases += unit.sequence.size();
            units.push_back(std::move(unit));
        }
        for (const auto& interval : segmentation.residuals) {
            const uint64_t length = interval.second - interval.first;
            residual_bases += length;
            residuals.push_back({array.chrom, array.start + interval.first,
                                 array.start + interval.second, array.id,
                                 segmentation.units.empty() ? "shorter_than_min_period"
                                                            : "edge_residual"});
        }
        segmentations[array_index] = std::move(segmentation);
    }
    if (complete_bases + residual_bases != input_bases) {
        throw std::runtime_error("internal ACRE segmentation coverage failure");
    }
    if (units.empty()) {
        throw std::runtime_error("no complete ACREs fit the configured period range");
    }
    std::cerr << "Segmented " << arrays.size() << " arrays into " << units.size()
              << " complete ACREs (" << complete_bases << " bp) and "
              << residuals.size() << " edge/short fragments (" << residual_bases
              << " bp) in " << seconds_since(segmentation_start) << '\n';

    // Complete ACRE sequences are now owned by `units`; release the duplicate
    // per-array strings before building k-mer structures.
    std::vector<std::string>().swap(array_sequences);

    // owner >= 0 means this canonical k-mer occurs in only that ACRE; -1 means
    // it has occurred in at least two segmented ACRE units.
    const auto feature_start = Clock::now();
    std::unordered_map<uint64_t, int32_t> owners;
    if (!options.pangenome) {
        owners.reserve(static_cast<size_t>(std::min<uint64_t>(complete_bases, 50000000ULL)));
        for (size_t i = 0; i < units.size(); ++i) {
            for_each_canonical_kmer(units[i].sequence, options.k,
                [&](size_t, uint64_t code, bool) {
                    const auto inserted = owners.emplace(code, static_cast<int32_t>(i));
                    if (!inserted.second && inserted.first->second != static_cast<int32_t>(i)) {
                        inserted.first->second = -1;
                    }
                });
        }
    }

    const size_t ref_count = reference_index.references.size();
    const size_t dimensions = ref_count * 2U + options.feature_bins + 2U;
    std::vector<std::vector<float>> features(units.size(),
                                             std::vector<float>(dimensions, 0));
    std::vector<RegionStats> stats(units.size());
    std::vector<KmerSketch> unit_sketches(options.deep_subfamilies ? units.size() : 0U);
    std::vector<std::vector<AnnotationSegment>> annotations(units.size());
    std::atomic<size_t> feature_next{0};
    std::atomic<bool> feature_failed{false};
    std::exception_ptr feature_error;
    std::mutex feature_error_mutex;
    auto compute_unit = [&](size_t i) {
        std::unordered_set<uint64_t> distinct;
        distinct.reserve(units[i].sequence.size() / 2U + 1U);
        std::vector<uint64_t> diagnostic_counts(ref_count, 0);
        for_each_canonical_kmer(units[i].sequence, options.k,
            [&](size_t, uint64_t code, bool) {
                ++stats[i].valid_kmers;
                distinct.insert(code);
                ++features[i][ref_count * 2U + mix64(code) % options.feature_bins];
                const auto match = reference_index.kmers.find(code);
                if (match != reference_index.kmers.end() &&
                    popcount64(match->second.members) == 1U) {
                    ++diagnostic_counts[trailing_zeroes64(match->second.members)];
                }
            });
        stats[i].distinct_kmers = distinct.size();
        if (!options.pangenome) {
            for (uint64_t code : distinct) {
                const auto owner = owners.find(code);
                if (owner != owners.end() && owner->second == static_cast<int32_t>(i)) {
                    ++stats[i].exclusive_kmers;
                }
            }
        }
        const double valid_denominator = stats[i].valid_kmers
                                             ? static_cast<double>(stats[i].valid_kmers)
                                             : 1.0;
        for (size_t ref = 0; ref < ref_count; ++ref) {
            features[i][ref] = static_cast<float>(diagnostic_counts[ref] / valid_denominator);
            features[i][ref_count + ref] = diagnostic_counts[ref] ? 1.0 : 0.0;
        }
        for (size_t bin = 0; bin < options.feature_bins; ++bin) {
            features[i][ref_count * 2U + bin] = static_cast<float>(
                features[i][ref_count * 2U + bin] / valid_denominator);
        }
        const double exclusive_fraction = !options.pangenome && stats[i].distinct_kmers
                                              ? static_cast<double>(stats[i].exclusive_kmers) /
                                                    static_cast<double>(stats[i].distinct_kmers)
                                              : 0.0;
        features[i][dimensions - 2U] = static_cast<float>(exclusive_fraction);
        features[i][dimensions - 1U] = !options.pangenome && stats[i].exclusive_kmers ? 1.0 : 0.0;
        if (options.deep_subfamilies) {
            unit_sketches[i] = make_bottom_sketch(units[i].sequence, options.k,
                                                   options.sketch_size,
                                                   options.sketch_position_bins);
        }
        annotations[i] = annotate_sequence(units[i].sequence, reference_index, options);
        std::vector<uint64_t> covered(ref_count, 0);
        for (const auto& segment : annotations[i]) {
            if (segment.reference >= 0) {
                covered[static_cast<size_t>(segment.reference)] +=
                    segment.relative_end - segment.relative_start;
            }
        }
        if (!covered.empty()) {
            const auto best = std::max_element(covered.begin(), covered.end());
            if (*best > 0) {
                const size_t ref = static_cast<size_t>(best - covered.begin());
                stats[i].top_component = reference_index.references[ref].name;
                stats[i].component_coverage = static_cast<double>(*best) /
                                              static_cast<double>(units[i].sequence.size());
            }
        }
    };
    std::vector<std::thread> feature_workers;
    size_t feature_thread_count = options.threads;
    if (feature_thread_count == 0U) {
        feature_thread_count = std::thread::hardware_concurrency();
        if (feature_thread_count == 0U) feature_thread_count = 1U;
        feature_thread_count = std::min<size_t>(feature_thread_count, 8U);
    }
    feature_thread_count = std::max<size_t>(1U, std::min(feature_thread_count, units.size()));
    feature_workers.reserve(feature_thread_count);
    for (size_t worker = 0; worker < feature_thread_count; ++worker) {
        feature_workers.emplace_back([&] {
            try {
                while (!feature_failed.load(std::memory_order_relaxed)) {
                    const size_t i = feature_next.fetch_add(1, std::memory_order_relaxed);
                    if (i >= units.size()) break;
                    compute_unit(i);
                }
            } catch (...) {
                feature_failed.store(true, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(feature_error_mutex);
                if (!feature_error) feature_error = std::current_exception();
            }
        });
    }
    for (auto& worker : feature_workers) worker.join();
    if (feature_error) std::rethrow_exception(feature_error);

    size_t cluster_count = options.clusters;
    if (cluster_count == 0) {
        cluster_count = static_cast<size_t>(std::llround(std::sqrt(units.size() / 2.0)));
        cluster_count = std::max<size_t>(2, std::min<size_t>(8, cluster_count));
        cluster_count = std::min(cluster_count, units.size());
    }
    const auto assignments = kmeans(features, cluster_count, options.seed);

    std::cerr << "Computed features and annotations in " << seconds_since(feature_start) << '\n';
    if (options.deep_subfamilies) {
        uint64_t retained_hashes = 0;
        for (const auto& sketch : unit_sketches) {
            retained_hashes += sketch.hashes.size() + sketch.reverse_hashes.size();
        }
        std::cerr << "Retained " << retained_hashes
                  << " phase-aware sketch hashes ("
                  << retained_hashes * sizeof(uint64_t) << " bytes) across "
                  << unit_sketches.size() << " complete ACREs\n";
        const auto deep_start = Clock::now();
        write_deep_subfamily_analysis(options, units, unit_sketches);
        std::cerr << "Completed switch-region deep subfamily analysis in "
                  << seconds_since(deep_start) << '\n';
    }

    const std::string annotations_path = options.output_prefix + ".annotations.bed";
    const std::string features_path = options.output_prefix + ".features.tsv";
    const std::string clusters_path = options.output_prefix + ".clusters.tsv";
    const std::string acres_path = options.output_prefix + ".acres.bed";
    const std::string partials_path = options.output_prefix + ".partials.bed";
    const std::string arrays_path = options.output_prefix + ".arrays.tsv";
    const std::string colors_path = options.output_prefix + ".colors.tsv";
    std::ofstream annotation_output(annotations_path);
    std::ofstream feature_output(features_path);
    std::ofstream cluster_output(clusters_path);
    std::ofstream acres_output(acres_path);
    std::ofstream partials_output(partials_path);
    std::ofstream arrays_output(arrays_path);
    std::ofstream colors_output(colors_path);
    if (!annotation_output || !feature_output || !cluster_output || !acres_output ||
        !partials_output || !arrays_output || !colors_output) {
        throw std::runtime_error("cannot create output files with prefix: " +
                                 options.output_prefix);
    }
    annotation_output << std::setprecision(6);
    for (size_t i = 0; i < units.size(); ++i) {
        for (const auto& segment : annotations[i]) {
            const std::string name = segment.reference >= 0
                ? reference_index.references[static_cast<size_t>(segment.reference)].name
                : "unclassified";
            const uint64_t hit_start = units[i].start + segment.relative_start;
            const uint64_t hit_end = units[i].start + segment.relative_end;
            const double density = segment.valid_kmers
                ? segment.weighted_hits / static_cast<double>(segment.valid_kmers) : 0;
            const double confidence = segment.valid_kmers
                ? segment.confidence_sum / static_cast<double>(segment.valid_kmers) : 0;
            const int bed_score = static_cast<int>(std::llround(
                std::min(1.0, std::max(0.0, density)) * 1000.0));
            annotation_output << units[i].chrom << '\t' << hit_start << '\t' << hit_end << '\t'
                              << name << '\t' << bed_score << '\t' << segment.strand << '\t'
                              << hit_start << '\t' << hit_end << '\t' << colors.at(name).text()
                              << '\t' << units[i].id << '\t' << units[i].array_id << '\t'
                              << segment.weighted_hits << '\t' << density << '\t' << confidence
                              << '\n';
        }
    }

    feature_output << "acre_id\tarray_id\tcluster";
    for (const auto& reference : reference_index.references) {
        feature_output << "\tdiagnostic_frequency:" << reference.name;
    }
    for (const auto& reference : reference_index.references) {
        feature_output << "\tdiagnostic_present:" << reference.name;
    }
    for (size_t bin = 0; bin < options.feature_bins; ++bin) {
        feature_output << "\thash_frequency:" << bin;
    }
    feature_output << "\texclusive_fraction\texclusive_present\n" << std::setprecision(9);
    for (size_t i = 0; i < units.size(); ++i) {
        feature_output << units[i].id << '\t' << units[i].array_id << '\t'
                       << assignments[i] + 1U;
        for (size_t dimension = 0; dimension < features[i].size(); ++dimension) {
            if (options.pangenome && dimension + 2U >= features[i].size()) {
                feature_output << "\tNA";
            } else {
                feature_output << '\t' << features[i][dimension];
            }
        }
        feature_output << '\n';
    }

    cluster_output << "acre_id\tarray_id\tunit_in_array\tchrom\tstart\tend\tlength\t"
                      "estimated_period\tperiodicity_score\tcluster\ttop_component\t"
                      "component_coverage\tvalid_kmers\tdistinct_kmers\texclusive_kmers\t"
                      "exclusive_fraction\tsource_name\n" << std::setprecision(9);
    for (size_t i = 0; i < units.size(); ++i) {
        const double fraction = !options.pangenome && stats[i].distinct_kmers
            ? static_cast<double>(stats[i].exclusive_kmers) /
                  static_cast<double>(stats[i].distinct_kmers) : 0;
        cluster_output << units[i].id << '\t' << units[i].array_id << '\t'
                       << units[i].unit_number << '\t' << units[i].chrom << '\t'
                       << units[i].start << '\t' << units[i].end << '\t'
                       << units[i].end - units[i].start << '\t' << units[i].estimated_period
                       << '\t' << units[i].periodicity_score << '\t' << assignments[i] + 1U
                       << '\t' << stats[i].top_component << '\t' << stats[i].component_coverage
                       << '\t' << stats[i].valid_kmers << '\t' << stats[i].distinct_kmers << '\t';
        if (options.pangenome) cluster_output << "NA\tNA\t";
        else cluster_output << stats[i].exclusive_kmers << '\t' << fraction << '\t';
        cluster_output << units[i].source_name << '\n';
    }

    acres_output << std::setprecision(9);
    for (const auto& unit : units) {
        acres_output << unit.chrom << '\t' << unit.start << '\t' << unit.end << '\t'
                     << unit.id << "\t0\t.\t" << unit.array_id << '\t' << unit.unit_number
                     << '\t' << unit.estimated_period << '\t' << unit.periodicity_score << '\t'
                     << unit.source_name << '\n';
    }
    for (size_t i = 0; i < residuals.size(); ++i) {
        partials_output << residuals[i].chrom << '\t' << residuals[i].start << '\t'
                        << residuals[i].end << "\tPARTIAL_" << std::setw(6)
                        << std::setfill('0') << i + 1U << std::setfill(' ') << "\t0\t.\t"
                        << residuals[i].array_id << '\t' << residuals[i].reason << '\n';
    }
    arrays_output << "array_id\tchrom\tstart\tend\tinput_length\tcomplete_acres\t"
                     "complete_bp\tresidual_bp\tsegmentation_method\testimated_period\t"
                     "periodicity_score\tsource_name\n" << std::setprecision(9);
    for (size_t i = 0; i < arrays.size(); ++i) {
        uint64_t complete = 0;
        uint64_t residual = 0;
        for (const auto& interval : segmentations[i].units) complete += interval.second - interval.first;
        for (const auto& interval : segmentations[i].residuals) residual += interval.second - interval.first;
        arrays_output << arrays[i].id << '\t' << arrays[i].chrom << '\t' << arrays[i].start
                      << '\t' << arrays[i].end << '\t' << arrays[i].end - arrays[i].start
                      << '\t' << segmentations[i].units.size() << '\t' << complete << '\t'
                      << residual << '\t' << segmentations[i].method << '\t'
                      << segmentations[i].estimated_period << '\t'
                      << segmentations[i].periodicity_score << '\t'
                      << arrays[i].source_name << '\n';
    }
    colors_output << "component\titemRgb\n";
    for (const auto& reference : reference_index.references) {
        colors_output << reference.name << '\t' << colors.at(reference.name).text() << '\n';
    }
    colors_output << "unclassified\t" << colors.at("unclassified").text() << '\n';

    if (!options.acro1_monomers_path.empty()) {
        write_acro1_subfamilies(options, reference_index, units);
    }

    if (options.dump_exclusive) {
        const std::string path = options.output_prefix + ".exclusive_kmers.tsv";
        std::ofstream output(path);
        if (!output) throw std::runtime_error("cannot create output: " + path);
        std::vector<std::pair<uint64_t, int32_t>> exclusive;
        for (const auto& item : owners) if (item.second >= 0) exclusive.push_back(item);
        std::sort(exclusive.begin(), exclusive.end(), [](const auto& a, const auto& b) {
            return a.second < b.second || (a.second == b.second && a.first < b.first);
        });
        output << "acre_id\tkmer\n";
        for (const auto& item : exclusive) {
            output << units[static_cast<size_t>(item.second)].id << '\t'
                   << decode_kmer(item.first, options.k) << '\n';
        }
    }

    if (options.write_fasta) {
        const std::string path = options.output_prefix + ".acres.fasta";
        std::ofstream output(path);
        if (!output) throw std::runtime_error("cannot create output: " + path);
        for (const auto& unit : units) {
            output << '>' << unit.id << "|array=" << unit.array_id << '|'
                   << unit.chrom << ':' << unit.start << '-' << unit.end << '\n';
            for (size_t offset = 0; offset < unit.sequence.size(); offset += 80U) {
                output.write(unit.sequence.data() + offset,
                             static_cast<std::streamsize>(
                                 std::min<size_t>(80U, unit.sequence.size() - offset)));
                output << '\n';
            }
        }
    }

    std::cerr << "Wrote per-ACRE annotations, features, clusters, unit BED, array QC, "
                 "partials, colors, and optional ACRO1 subfamilies (total "
              << seconds_since(overall_start) << ")\n";
    return 0;
}

int self_test() {
    try {
        std::vector<std::pair<uint64_t, bool>> forward;
        std::vector<std::pair<uint64_t, bool>> reverse;
        for_each_canonical_kmer("AACCGTT", 5,
            [&](size_t, uint64_t code, bool orientation) {
                forward.emplace_back(code, orientation);
            });
        for_each_canonical_kmer("AACGGTT", 5,
            [&](size_t, uint64_t code, bool orientation) {
                reverse.emplace_back(code, orientation);
            });
        if (forward.size() != 3 || reverse.size() != 3 ||
            forward[0].first != reverse[2].first ||
            forward[1].first != reverse[1].first ||
            forward[2].first != reverse[0].first) {
            throw std::runtime_error("canonical reverse-complement test failed");
        }
        size_t valid = 0;
        for_each_canonical_kmer("AAAANAAAA", 4,
            [&](size_t, uint64_t, bool) { ++valid; });
        if (valid != 2) throw std::runtime_error("ambiguous-base reset test failed");

        auto index = build_reference_index({{"one", "AAAACCCCGGGG"},
                                            {"two", "TTTTCCCCAAAA"}}, 5);
        if (index.references.size() != 2 || index.kmers.empty() ||
            index.diagnostic_kmers[0] == 0 || index.diagnostic_kmers[1] == 0) {
            throw std::runtime_error("reference index test failed");
        }
        const auto clustered = kmeans({{0, 0}, {0.01, 0}, {10, 10}, {10.01, 10}}, 2, 1);
        if (clustered[0] != clustered[1] || clustered[2] != clustered[3] ||
            clustered[0] == clustered[2]) {
            throw std::runtime_error("clustering test failed");
        }
        Options segmentation_options;
        segmentation_options.k = 3;
        segmentation_options.min_period = 4;
        segmentation_options.max_period = 4;
        segmentation_options.period_target = 4;
        segmentation_options.boundary_span = 4;
        const auto segmented = segment_sequence("ACGTACGTACGT", segmentation_options);
        if (segmented.units.size() != 3U || !segmented.residuals.empty() ||
            segmented.units.front() != std::make_pair<size_t, size_t>(0, 4) ||
            segmented.units.back() != std::make_pair<size_t, size_t>(8, 12)) {
            throw std::runtime_error("ACRE segmentation test failed");
        }
        std::string acro1_test_reference;
        const char bases[] = {'A', 'C', 'G', 'T'};
        for (size_t i = 0; i < 147U; ++i) {
            acro1_test_reference.push_back(bases[mix64(i + 17U) % 4U]);
        }
        const std::string internal = acro1_test_reference.substr(0, 5) +
                                     acro1_test_reference.substr(8);
        const std::string terminal = acro1_test_reference.substr(0, 144);
        size_t internal_mismatches = 0;
        size_t terminal_mismatches = 0;
        size_t compared = 0;
        if (classify_acro1_144(internal, internal, terminal, internal_mismatches,
                               terminal_mismatches, compared) != 2U ||
            classify_acro1_144(terminal, internal, terminal, internal_mismatches,
                               terminal_mismatches, compared) != 3U) {
            throw std::runtime_error("ACRO1 subfamily classifier test failed");
        }
        KmerSketch first_sketch{3, 3, {1, 2, 3}, 0, {}};
        KmerSketch second_sketch{3, 3, {2, 3, 4}, 0, {}};
        if (std::abs(sketch_jaccard(first_sketch, second_sketch) - 0.5) > 1e-12) {
            throw std::runtime_error("bottom-k Jaccard test failed");
        }
        const std::string phase_test = "ACGTCAGTACGATCGATTGGCCTATGCAAGCT";
        const auto phase_forward = make_bottom_sketch(phase_test, 5, 64, 8);
        const auto phase_reverse = make_bottom_sketch(reverse_complement(phase_test), 5, 64, 8);
        if (std::abs(sketch_jaccard(phase_forward, phase_reverse) - 1.0) > 1e-12) {
            throw std::runtime_error("phase-aware reverse-complement sketch test failed");
        }
        std::string truncation_test;
        for (size_t i = 0; i < 512U; ++i) {
            truncation_test.push_back(bases[mix64(i * 31U + 7U) % 4U]);
        }
        const auto truncated = make_bottom_sketch(truncation_test, 7, 32, 16);
        if (truncated.hashes.size() != 32U || truncated.reverse_hashes.size() != 32U ||
            truncated.distinct_kmers <= truncated.hashes.size() ||
            truncated.reverse_distinct_kmers <= truncated.reverse_hashes.size()) {
            throw std::runtime_error("fixed-size truncated sketch test failed");
        }
        std::vector<KmerSketch> switch_sketches;
        for (size_t i = 0; i < 3U; ++i) {
            switch_sketches.push_back({5, 5, {1, 2, 3, 4, 5}, 0, {}});
        }
        for (size_t i = 0; i < 3U; ++i) {
            switch_sketches.push_back({5, 5, {101, 102, 103, 104, 105}, 0, {}});
        }
        SwitchThresholds switch_thresholds;
        switch_thresholds.sketch_size = 5;
        const auto switch_calls = detect_switch_boundaries(switch_sketches, 5,
                                                            switch_thresholds);
        const auto collapsed_switches = collapse_switch_boundaries(switch_calls,
                                                                    switch_thresholds);
        if (switch_calls.size() != 1U || !switch_calls.front().candidate ||
            collapsed_switches.size() != 1U ||
            collapsed_switches.front().boundary_index != 3U) {
            throw std::runtime_error("deep switch detection test failed");
        }
        std::vector<KmerSketch> no_switch_sketches(6U, switch_sketches.front());
        const auto no_switch_calls = detect_switch_boundaries(no_switch_sketches, 5,
                                                               switch_thresholds);
        if (no_switch_calls.size() != 1U || no_switch_calls.front().candidate ||
            !collapse_switch_boundaries(no_switch_calls, switch_thresholds).empty()) {
            throw std::runtime_error("deep no-switch test failed");
        }
        SwitchBoundary weaker;
        weaker.boundary_index = 4U;
        weaker.ani_drop = 0.1;
        weaker.candidate = true;
        SwitchBoundary stronger = weaker;
        stronger.boundary_index = 3U;
        stronger.ani_drop = 0.2;
        const auto nearby = collapse_switch_boundaries({weaker, stronger}, 3U);
        if (nearby.size() != 1U || nearby.front().boundary_index != 3U) {
            throw std::runtime_error("nearby switch collapse test failed");
        }
        std::cout << "All N-TRIAGE self-tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Self-test failed: " << error.what() << '\n';
        return 1;
    }
}

}  // namespace ntriage
