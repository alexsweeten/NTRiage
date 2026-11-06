#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>

namespace ntriage {

inline int base_bits(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

// Calls fn(start, canonical_code, forward_is_canonical) for every valid k-mer.
// Ambiguous bases break the rolling k-mer and are otherwise ignored.
template <typename Fn>
inline void for_each_canonical_kmer(std::string_view sequence, unsigned k, Fn&& fn) {
    if (k == 0 || k > 31) {
        throw std::invalid_argument("k-mer size must be in [1, 31]");
    }
    const uint64_t mask = (uint64_t{1} << (2U * k)) - 1U;
    const unsigned top_shift = 2U * (k - 1U);
    uint64_t forward = 0;
    uint64_t reverse = 0;
    unsigned valid = 0;
    for (size_t i = 0; i < sequence.size(); ++i) {
        const int bits = base_bits(sequence[i]);
        if (bits < 0) {
            forward = reverse = 0;
            valid = 0;
            continue;
        }
        forward = ((forward << 2U) | static_cast<uint64_t>(bits)) & mask;
        reverse = (reverse >> 2U) |
                  (static_cast<uint64_t>(3 - bits) << top_shift);
        if (++valid >= k) {
            const bool forward_is_canonical = forward <= reverse;
            fn(i + 1U - k,
               forward_is_canonical ? forward : reverse,
               forward_is_canonical);
        }
    }
}

inline std::string decode_kmer(uint64_t code, unsigned k) {
    static constexpr char bases[] = "ACGT";
    std::string result(k, 'A');
    for (unsigned i = 0; i < k; ++i) {
        result[k - 1U - i] = bases[code & 3U];
        code >>= 2U;
    }
    return result;
}

}  // namespace ntriage
