#pragma once

#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace ntriage {

struct FaiRecord {
    uint64_t length = 0;
    uint64_t offset = 0;
    uint64_t line_bases = 0;
    uint64_t line_width = 0;
};

struct GziPoint {
    uint64_t compressed = 0;
    uint64_t uncompressed = 0;
};

struct GenomeIndex {
    std::string fasta_path;
    bool compressed = false;
    std::unordered_map<std::string, FaiRecord> fai;
    std::vector<GziPoint> gzi;

    explicit GenomeIndex(std::string path);
};

// One reader should be used by one thread. GenomeIndex itself is immutable and
// may be shared by all worker threads.
class FastaReader {
public:
    explicit FastaReader(std::shared_ptr<const GenomeIndex> index);
    ~FastaReader();
    FastaReader(const FastaReader&) = delete;
    FastaReader& operator=(const FastaReader&) = delete;

    std::string fetch(const std::string& contig, uint64_t start, uint64_t end);

private:
    struct Block {
        std::vector<unsigned char> data;
        uint64_t compressed_size = 0;
    };

    std::shared_ptr<const GenomeIndex> index_;
    std::FILE* file_ = nullptr;

    Block read_bgzf_block(uint64_t compressed_offset);
    std::vector<unsigned char> read_uncompressed(uint64_t offset, size_t length);
};

}  // namespace ntriage
