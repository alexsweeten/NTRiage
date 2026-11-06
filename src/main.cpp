#include "ntriage/core.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef NTRIAGE_VERSION
#define NTRIAGE_VERSION "development"
#endif

namespace {

void usage(std::ostream& out) {
    out <<
        "N-TRIAGE - nested tandem repeat inference & analysis of genomic evolution\n"
        "Fast k-mer annotation, switch detection, and subfamily analysis of ACRE arrays\n\n"
        "Usage:\n"
        "  ntriage --bed regions.bed --fasta genome.fa[.gz] --out PREFIX \\\n"
        "          [--refs components.fa] [options]\n\n"
        "Required:\n"
        "  --bed FILE             BED array intervals (0-based, half-open)\n"
        "  --fasta FILE           indexed FASTA; BGZF requires .fai and .gzi\n"
        "  --refs FILE            ACRE component reference FASTA (non-pangenome runs)\n"
        "  -o, --out PREFIX       output file prefix\n\n"
        "Options:\n"
        "  -k, --kmer INT         canonical k-mer size [15]\n"
        "  --bins INT             hashed frequency features [256]\n"
        "  --window INT           annotation window in bp [200]\n"
        "  --min-period INT       minimum complete ACRE length [6200]\n"
        "  --max-period INT       maximum complete ACRE length [7100]\n"
        "  --period-target INT    preferred period when evidence is weak [6650]\n"
        "  --boundary-span INT    bp compared to refine boundaries [1024]\n"
        "  --colors FILE          optional component-to-RGB TSV overrides\n"
        "  --acro1-monomers FILE  ACRO1 monomer BED; emit subfamily 1/2/3 calls\n"
        "  -t, --threads INT      extraction threads [auto, capped at 8]\n"
        "  -c, --clusters INT     k-means clusters [auto heuristic]\n"
        "  --min-hits FLOAT       minimum weighted hits/window [3]\n"
        "  --min-density FLOAT    minimum weighted hit density [0.02]\n"
        "  --seed INT             deterministic clustering seed [1]\n"
        "  --deep-subfamilies     detect intra-array switches and infer local subfamilies\n"
        "  --sketch-size INT      bottom-k hashes per orientation/complete ACRE [512]\n"
        "  --sketch-position-bins INT  coarse phase bins per complete ACRE [16]\n"
        "  --switch-flank INT     complete ACREs on each side of a switch [3]\n"
        "  --switch-radius INT    complete ACREs reported around each switch [3]\n"
        "  --switch-coherence F   minimum within-flank sketch identity [0.90]\n"
        "  --switch-cross-ani F   maximum cross-flank sketch identity [0.92]\n"
        "  --switch-drop F        minimum coherence-minus-cross identity [0.04]\n"
        "  --switch-adjacent F    maximum adjacent-unit sketch identity [0.92]\n"
        "  --subfamily-ani F      minimum sketch identity for one global subfamily [0.94]\n"
        "  --pangenome            streaming switch mode; retain only called contexts\n"
        "  --dump-exclusive       write every ACRE-exclusive k-mer\n"
        "  --write-fasta          write complete segmented ACRE sequences\n"
        "  --self-test            run built-in unit checks\n"
        "  --version              show the N-TRIAGE version\n"
        "  -h, --help             show this help\n";
}

uint64_t parse_u64(const std::string& value, const std::string& option) {
    size_t used = 0;
    uint64_t parsed = 0;
    try {
        parsed = std::stoull(value, &used);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid value for " + option + ": " + value);
    }
    if (used != value.size()) throw std::runtime_error("invalid value for " + option);
    return parsed;
}

double parse_double(const std::string& value, const std::string& option) {
    size_t used = 0;
    double parsed = 0;
    try {
        parsed = std::stod(value, &used);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid value for " + option + ": " + value);
    }
    if (used != value.size()) throw std::runtime_error("invalid value for " + option);
    return parsed;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        ntriage::Options options;
        bool do_self_test = false;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            auto next = [&]() -> std::string {
                if (++i >= argc) throw std::runtime_error("missing value after " + arg);
                return argv[i];
            };
            if (arg == "--bed") options.bed_path = next();
            else if (arg == "--fasta") options.fasta_path = next();
            else if (arg == "--refs") options.refs_path = next();
            else if (arg == "--colors") options.colors_path = next();
            else if (arg == "--acro1-monomers") options.acro1_monomers_path = next();
            else if (arg == "-o" || arg == "--out") options.output_prefix = next();
            else if (arg == "-k" || arg == "--kmer") options.k = parse_u64(next(), arg);
            else if (arg == "--bins") options.feature_bins = parse_u64(next(), arg);
            else if (arg == "--window") options.window = parse_u64(next(), arg);
            else if (arg == "--min-period") options.min_period = parse_u64(next(), arg);
            else if (arg == "--max-period") options.max_period = parse_u64(next(), arg);
            else if (arg == "--period-target") options.period_target = parse_u64(next(), arg);
            else if (arg == "--boundary-span") options.boundary_span = parse_u64(next(), arg);
            else if (arg == "-t" || arg == "--threads") options.threads = parse_u64(next(), arg);
            else if (arg == "-c" || arg == "--clusters") options.clusters = parse_u64(next(), arg);
            else if (arg == "--min-hits") options.min_hits = parse_double(next(), arg);
            else if (arg == "--min-density") options.min_density = parse_double(next(), arg);
            else if (arg == "--seed") options.seed = parse_u64(next(), arg);
            else if (arg == "--deep-subfamilies") options.deep_subfamilies = true;
            else if (arg == "--pangenome") options.pangenome = true;
            else if (arg == "--sketch-size") options.sketch_size = parse_u64(next(), arg);
            else if (arg == "--sketch-position-bins") {
                options.sketch_position_bins = parse_u64(next(), arg);
            }
            else if (arg == "--switch-flank") options.switch_flank = parse_u64(next(), arg);
            else if (arg == "--switch-radius") options.switch_radius = parse_u64(next(), arg);
            else if (arg == "--switch-coherence") {
                options.switch_min_coherence = parse_double(next(), arg);
            }
            else if (arg == "--switch-cross-ani") {
                options.switch_max_cross_ani = parse_double(next(), arg);
            }
            else if (arg == "--switch-drop") {
                options.switch_min_drop = parse_double(next(), arg);
            }
            else if (arg == "--switch-adjacent") {
                options.switch_max_adjacent_ani = parse_double(next(), arg);
            }
            else if (arg == "--subfamily-ani") {
                options.subfamily_min_ani = parse_double(next(), arg);
            }
            else if (arg == "--dump-exclusive") options.dump_exclusive = true;
            else if (arg == "--write-fasta") options.write_fasta = true;
            else if (arg == "--self-test") do_self_test = true;
            else if (arg == "--version") {
                std::cout << "N-TRIAGE " << NTRIAGE_VERSION << '\n';
                return 0;
            }
            else if (arg == "-h" || arg == "--help") {
                usage(std::cout);
                return 0;
            } else {
                throw std::runtime_error("unknown option: " + arg);
            }
        }
        if (do_self_test) return ntriage::self_test();
        if (options.bed_path.empty() || options.fasta_path.empty() ||
            options.output_prefix.empty() ||
            (!options.pangenome && options.refs_path.empty())) {
            usage(std::cerr);
            throw std::runtime_error(
                "--bed, --fasta, and --out are required; --refs is required outside "
                "--pangenome mode");
        }
        if (options.k == 0 || options.k > 31) throw std::runtime_error("--kmer must be 1..31");
        if (options.feature_bins == 0) throw std::runtime_error("--bins must be positive");
        if (options.window < options.k) throw std::runtime_error("--window must be >= k");
        if (options.min_period < options.k || options.max_period < options.min_period) {
            throw std::runtime_error("period range must satisfy k <= min-period <= max-period");
        }
        if (options.period_target < options.min_period ||
            options.period_target > options.max_period) {
            throw std::runtime_error("--period-target must be inside the period range");
        }
        if (options.boundary_span < options.k) {
            throw std::runtime_error("--boundary-span must be >= k");
        }
        if (options.min_hits < 0 || options.min_density < 0) {
            throw std::runtime_error("minimum thresholds cannot be negative");
        }
        if (options.sketch_size < 32U) {
            throw std::runtime_error("--sketch-size must be at least 32");
        }
        if (options.sketch_position_bins == 0U) {
            throw std::runtime_error("--sketch-position-bins must be positive");
        }
        if (options.switch_flank < 2U) {
            throw std::runtime_error("--switch-flank must be at least 2");
        }
        if (options.switch_radius == 0U) {
            throw std::runtime_error("--switch-radius must be positive");
        }
        const auto probability = [](double value) { return value >= 0.0 && value <= 1.0; };
        if (!probability(options.switch_min_coherence) ||
            !probability(options.switch_max_cross_ani) ||
            !probability(options.switch_min_drop) ||
            !probability(options.switch_max_adjacent_ani) ||
            !probability(options.subfamily_min_ani)) {
            throw std::runtime_error("identity and switch thresholds must be between 0 and 1");
        }
        if (options.pangenome && options.dump_exclusive) {
            throw std::runtime_error("--pangenome and --dump-exclusive are incompatible");
        }
        if (options.pangenome && !options.deep_subfamilies) {
            throw std::runtime_error("--pangenome requires --deep-subfamilies");
        }
        if (options.pangenome && options.write_fasta) {
            throw std::runtime_error("--pangenome and --write-fasta are incompatible");
        }
        if (options.pangenome && !options.acro1_monomers_path.empty()) {
            throw std::runtime_error(
                "--pangenome and --acro1-monomers are incompatible; run ACRO1 "
                "classification separately");
        }
        return ntriage::run_pipeline(options);
    } catch (const std::exception& error) {
        std::cerr << "ntriage: error: " << error.what() << '\n';
        return 2;
    }
}
