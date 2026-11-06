#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace ntriage {

struct BedRecord {
    std::string chrom;
    uint64_t start = 0;
    uint64_t end = 0;
    std::string source_name;
    std::string id;
};

struct Reference {
    std::string name;
    std::string sequence;
};

struct KmerEntry {
    uint64_t members = 0;
    uint64_t reference_forward = 0;
    uint64_t reference_reverse = 0;
};

struct ReferenceIndex {
    std::vector<Reference> references;
    std::unordered_map<uint64_t, KmerEntry> kmers;
    std::vector<uint64_t> diagnostic_kmers;
};

struct Options {
    std::string bed_path;
    std::string fasta_path;
    std::string refs_path;
    std::string colors_path;
    std::string acro1_monomers_path;
    std::string output_prefix;
    unsigned k = 15;
    size_t feature_bins = 256;
    size_t window = 200;
    size_t threads = 0;
    size_t clusters = 0;
    double min_hits = 3.0;
    double min_density = 0.02;
    size_t min_period = 6200;
    size_t max_period = 7100;
    size_t period_target = 6650;
    size_t boundary_span = 1024;
    uint64_t seed = 1;
    bool dump_exclusive = false;
    bool write_fasta = false;

    // Deep switch/subfamily analysis uses fixed-size bottom-k sketches. Standard
    // runs retain one per ACRE; pangenome mode streams bounded array batches and
    // keeps only the sketches required by called switch contexts.
    bool deep_subfamilies = false;
    bool pangenome = false;
    size_t sketch_size = 512;
    size_t sketch_position_bins = 16;
    size_t switch_flank = 3;
    size_t switch_radius = 3;
    double switch_min_coherence = 0.90;
    double switch_max_cross_ani = 0.92;
    double switch_min_drop = 0.04;
    double switch_max_adjacent_ani = 0.92;
    double subfamily_min_ani = 0.94;
};

std::vector<BedRecord> read_bed(const std::string& path);
std::vector<Reference> read_fasta(const std::string& path);
ReferenceIndex build_reference_index(std::vector<Reference> references, unsigned k);
int run_pipeline(const Options& options);
int self_test();

}  // namespace ntriage
