#include "ntriage/fasta.hpp"

#include <zlib.h>

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace ntriage {
namespace {

uint64_t read_le64(std::istream& input) {
    unsigned char bytes[8]{};
    input.read(reinterpret_cast<char*>(bytes), 8);
    if (!input) throw std::runtime_error("truncated .gzi index");
    uint64_t value = 0;
    for (unsigned i = 0; i < 8; ++i) value |= uint64_t{bytes[i]} << (8U * i);
    return value;
}

uint16_t read_le16(const unsigned char* bytes) {
    return static_cast<uint16_t>(bytes[0]) |
           static_cast<uint16_t>(bytes[1] << 8U);
}

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

uint64_t parse_u64(const std::string& text, const std::string& context) {
    size_t used = 0;
    uint64_t value = 0;
    try {
        value = std::stoull(text, &used);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid integer in " + context + ": " + text);
    }
    if (used != text.size()) {
        throw std::runtime_error("invalid integer in " + context + ": " + text);
    }
    return value;
}

}  // namespace

GenomeIndex::GenomeIndex(std::string path) : fasta_path(std::move(path)) {
    std::ifstream fasta(fasta_path, std::ios::binary);
    if (!fasta) throw std::runtime_error("cannot open FASTA: " + fasta_path);
    unsigned char magic[2]{};
    fasta.read(reinterpret_cast<char*>(magic), 2);
    compressed = fasta.gcount() == 2 && magic[0] == 0x1f && magic[1] == 0x8b;

    std::ifstream fai_input(fasta_path + ".fai");
    if (!fai_input) {
        throw std::runtime_error("missing FASTA index: " + fasta_path + ".fai");
    }
    std::string line;
    size_t line_number = 0;
    while (std::getline(fai_input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 5) {
            throw std::runtime_error("invalid .fai line " + std::to_string(line_number));
        }
        FaiRecord record;
        const std::string context = fasta_path + ".fai line " + std::to_string(line_number);
        record.length = parse_u64(fields[1], context);
        record.offset = parse_u64(fields[2], context);
        record.line_bases = parse_u64(fields[3], context);
        record.line_width = parse_u64(fields[4], context);
        if (record.line_bases == 0 || record.line_width < record.line_bases) {
            throw std::runtime_error("invalid line geometry in " + context);
        }
        fai.emplace(fields[0], record);
    }
    if (fai.empty()) throw std::runtime_error("empty FASTA index: " + fasta_path + ".fai");

    if (compressed) {
        std::ifstream gzi_input(fasta_path + ".gzi", std::ios::binary);
        if (!gzi_input) {
            throw std::runtime_error("compressed FASTA requires BGZF index: " +
                                     fasta_path + ".gzi");
        }
        const uint64_t count = read_le64(gzi_input);
        if (count > 100000000ULL) throw std::runtime_error("implausible .gzi entry count");
        gzi.reserve(static_cast<size_t>(count) + 1U);
        gzi.push_back({0, 0});
        for (uint64_t i = 0; i < count; ++i) {
            GziPoint point{read_le64(gzi_input), read_le64(gzi_input)};
            if (point.compressed == 0 && point.uncompressed == 0) continue;
            if (point.compressed <= gzi.back().compressed ||
                point.uncompressed <= gzi.back().uncompressed) {
                throw std::runtime_error("unsorted .gzi index: " + fasta_path + ".gzi");
            }
            gzi.push_back(point);
        }
    }
}

FastaReader::FastaReader(std::shared_ptr<const GenomeIndex> index)
    : index_(std::move(index)) {
    file_ = std::fopen(index_->fasta_path.c_str(), "rb");
    if (!file_) {
        throw std::runtime_error("cannot open FASTA: " + index_->fasta_path +
                                 ": " + std::strerror(errno));
    }
}

FastaReader::~FastaReader() {
    if (file_) std::fclose(file_);
}

FastaReader::Block FastaReader::read_bgzf_block(uint64_t compressed_offset) {
    if (compressed_offset > static_cast<uint64_t>(std::numeric_limits<off_t>::max()) ||
        fseeko(file_, static_cast<off_t>(compressed_offset), SEEK_SET) != 0) {
        throw std::runtime_error("cannot seek BGZF FASTA");
    }
    unsigned char fixed[12]{};
    if (std::fread(fixed, 1, sizeof(fixed), file_) != sizeof(fixed)) {
        throw std::runtime_error("truncated BGZF block header");
    }
    if (fixed[0] != 0x1f || fixed[1] != 0x8b || fixed[2] != 8 || !(fixed[3] & 4)) {
        throw std::runtime_error("FASTA is gzip but not valid BGZF");
    }
    const uint16_t extra_length = read_le16(fixed + 10);
    std::vector<unsigned char> extra(extra_length);
    if (std::fread(extra.data(), 1, extra.size(), file_) != extra.size()) {
        throw std::runtime_error("truncated BGZF extra header");
    }
    uint64_t block_size = 0;
    for (size_t i = 0; i + 4 <= extra.size();) {
        const uint16_t sub_length = read_le16(extra.data() + i + 2);
        if (i + 4U + sub_length > extra.size()) break;
        if (extra[i] == 'B' && extra[i + 1] == 'C' && sub_length == 2) {
            block_size = static_cast<uint64_t>(read_le16(extra.data() + i + 4)) + 1U;
            break;
        }
        i += 4U + sub_length;
    }
    const size_t header_size = sizeof(fixed) + extra.size();
    if (block_size < header_size + 8U || block_size > 65536U) {
        throw std::runtime_error("invalid BGZF block size");
    }
    std::vector<unsigned char> compressed(static_cast<size_t>(block_size));
    std::copy(std::begin(fixed), std::end(fixed), compressed.begin());
    std::copy(extra.begin(), extra.end(), compressed.begin() + sizeof(fixed));
    const size_t remainder = compressed.size() - header_size;
    if (std::fread(compressed.data() + header_size, 1, remainder, file_) != remainder) {
        throw std::runtime_error("truncated BGZF block");
    }

    Block block;
    block.data.resize(65536U);
    z_stream stream{};
    stream.next_in = compressed.data();
    stream.avail_in = static_cast<uInt>(compressed.size());
    stream.next_out = block.data.data();
    stream.avail_out = static_cast<uInt>(block.data.size());
    if (inflateInit2(&stream, 47) != Z_OK) throw std::runtime_error("zlib init failed");
    const int status = inflate(&stream, Z_FINISH);
    inflateEnd(&stream);
    if (status != Z_STREAM_END) throw std::runtime_error("cannot inflate BGZF block");
    block.data.resize(stream.total_out);
    block.compressed_size = block_size;
    return block;
}

std::vector<unsigned char> FastaReader::read_uncompressed(uint64_t offset, size_t length) {
    std::vector<unsigned char> result;
    result.reserve(length);
    if (!index_->compressed) {
        if (offset > static_cast<uint64_t>(std::numeric_limits<off_t>::max()) ||
            fseeko(file_, static_cast<off_t>(offset), SEEK_SET) != 0) {
            throw std::runtime_error("cannot seek FASTA");
        }
        result.resize(length);
        const size_t got = std::fread(result.data(), 1, length, file_);
        result.resize(got);
        return result;
    }

    const auto it = std::upper_bound(
        index_->gzi.begin(), index_->gzi.end(), offset,
        [](uint64_t value, const GziPoint& point) { return value < point.uncompressed; });
    const GziPoint point = (it == index_->gzi.begin()) ? index_->gzi.front() : *(it - 1);
    uint64_t compressed_offset = point.compressed;
    uint64_t uncompressed_offset = point.uncompressed;
    while (result.size() < length) {
        Block block = read_bgzf_block(compressed_offset);
        if (block.data.empty()) break;
        const uint64_t block_end = uncompressed_offset + block.data.size();
        if (offset < block_end) {
            const size_t within = offset > uncompressed_offset
                                      ? static_cast<size_t>(offset - uncompressed_offset)
                                      : 0U;
            const size_t wanted = length - result.size();
            const size_t available = block.data.size() - within;
            const size_t take = std::min(wanted, available);
            result.insert(result.end(), block.data.begin() + within,
                          block.data.begin() + within + take);
            offset += take;
        }
        compressed_offset += block.compressed_size;
        uncompressed_offset = block_end;
    }
    return result;
}

std::string FastaReader::fetch(const std::string& contig, uint64_t start, uint64_t end) {
    const auto found = index_->fai.find(contig);
    if (found == index_->fai.end()) {
        throw std::runtime_error("BED contig not found in FASTA index: " + contig);
    }
    const FaiRecord& record = found->second;
    if (start > end || end > record.length) {
        throw std::runtime_error("interval outside FASTA contig " + contig + ":" +
                                 std::to_string(start) + "-" + std::to_string(end));
    }
    if (start == end) return {};
    const uint64_t raw_start = record.offset +
                               (start / record.line_bases) * record.line_width +
                               (start % record.line_bases);
    const uint64_t last = end - 1U;
    const uint64_t raw_end = record.offset +
                             (last / record.line_bases) * record.line_width +
                             (last % record.line_bases) + 1U;
    if (raw_end < raw_start || raw_end - raw_start > std::numeric_limits<size_t>::max()) {
        throw std::runtime_error("interval too large to load");
    }
    const auto raw = read_uncompressed(raw_start, static_cast<size_t>(raw_end - raw_start));
    std::string sequence;
    sequence.reserve(static_cast<size_t>(end - start));
    for (unsigned char c : raw) {
        if (c != '\n' && c != '\r') sequence.push_back(static_cast<char>(c));
    }
    if (sequence.size() != end - start) {
        throw std::runtime_error("short FASTA read for " + contig + ":" +
                                 std::to_string(start) + "-" + std::to_string(end));
    }
    return sequence;
}

}  // namespace ntriage
