// del3d - the pseudo-random generator behind the spatial sort's shuffle.
//
// spatial_sort() shuffles the point indices before Hilbert-sorting them, so
// this generator fixes the insertion order and therefore the triangulation's
// cell order. It has to be specified down to the bit - substituting another
// "uniform" generator changes the output on degenerate input - which is why the
// engine, the uniform-integer reduction and the shuffle loop are all written
// out here instead of taken from <random>, whose distributions and shuffle are
// implementation-defined.
#ifndef DEL3D_RAND48_H
#define DEL3D_RAND48_H

#include <algorithm>   // std::iter_swap
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace del3d {
namespace detail {

/// A 48-bit truncated linear congruential engine:
///     x <- (0x5DEECE66D * x + 0xB) mod 2^48
/// returning the top 31 bits of the state. Seeding places the seed in the high
/// bits, x = (seed << 16) | 0x330E, the classic drand48 initialisation.
class Rand48 {
public:
    using result_type = std::uint32_t;

    Rand48() { seed(1); }
    explicit Rand48(result_type v) { seed(v); }

    void seed(result_type v) {
        x_ = ((static_cast<std::uint64_t>(v) << 16) | 0x330EULL) & kMask;
    }

    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return 0x7FFFFFFF; }

    result_type operator()() {
        x_ = (kA * x_ + kC) & kMask;
        return static_cast<result_type>(x_ >> 17);
    }

private:
    static constexpr std::uint64_t kA    = 0x5DEECE66DULL;
    static constexpr std::uint64_t kC    = 0xBULL;
    static constexpr std::uint64_t kMask = (1ULL << 48) - 1ULL;
    std::uint64_t x_ = 0;
};

/**
 * Maps engine output uniformly onto [min_value, max_value].
 *
 * Three cases, by how the engine's range compares with the requested one:
 *   * equal        - one draw, shifted into place;
 *   * engine range smaller - concatenate several draws in mixed radix,
 *     rejecting any combination that overflows or lands out of range;
 *   * engine range larger  - divide the draw into equal-sized buckets and
 *     reject the leftover partial bucket, which keeps the result unbiased.
 *
 * Only the last case is reachable with Rand48 and the ranges used here; the
 * others are kept so the reduction stays correct if the engine is changed.
 */
template <class Engine, class T>
T generate_uniform_int(Engine& eng, T min_value, T max_value) {
    using range_type = std::make_unsigned_t<T>;
    using base_unsigned = std::make_unsigned_t<typename Engine::result_type>;

    const range_type range = static_cast<range_type>(
        static_cast<range_type>(max_value) - static_cast<range_type>(min_value));
    const base_unsigned bmin = static_cast<base_unsigned>(Engine::min());
    const base_unsigned brange =
        static_cast<base_unsigned>(static_cast<base_unsigned>(Engine::max()) - bmin);

    if (range == 0) return min_value;

    if (static_cast<range_type>(brange) == range) {
        const base_unsigned v = static_cast<base_unsigned>(eng()) - bmin;
        return static_cast<T>(static_cast<range_type>(v) + static_cast<range_type>(min_value));
    }

    if (static_cast<range_type>(brange) < range) {
        // Concatenate several draws, rejecting out-of-range results.
        for (;;) {
            range_type limit;
            if (range == (std::numeric_limits<range_type>::max)()) {
                limit = range / (static_cast<range_type>(brange) + 1);
                if (range % (static_cast<range_type>(brange) + 1) == static_cast<range_type>(brange))
                    ++limit;
            } else {
                limit = (range + 1) / (static_cast<range_type>(brange) + 1);
            }
            range_type result = 0, mult = 1;
            while (mult <= limit) {
                result += static_cast<range_type>(
                    static_cast<range_type>(static_cast<base_unsigned>(eng()) - bmin) * mult);
                if (mult * static_cast<range_type>(brange) == range - mult + 1)
                    return static_cast<T>(result);
                mult *= static_cast<range_type>(brange) + 1;
            }
            range_type inc = generate_uniform_int(
                eng, static_cast<range_type>(0), static_cast<range_type>(range / mult));
            if ((std::numeric_limits<range_type>::max)() / mult < inc) continue;
            inc *= mult;
            result += inc;
            if (result < inc) continue;          // overflowed
            if (result > range) continue;        // out of range
            return static_cast<T>(result + static_cast<range_type>(min_value));
        }
    }

    // brange > range: bucket the engine's output.
    base_unsigned bucket;
    if (brange == (std::numeric_limits<base_unsigned>::max)()) {
        bucket = brange / (static_cast<base_unsigned>(range) + 1);
        if (brange % (static_cast<base_unsigned>(range) + 1) == static_cast<base_unsigned>(range))
            ++bucket;
    } else {
        bucket = (brange + 1) / (static_cast<base_unsigned>(range) + 1);
    }
    for (;;) {
        base_unsigned r = static_cast<base_unsigned>(eng()) - bmin;
        r /= bucket;
        if (r <= static_cast<base_unsigned>(range))
            return static_cast<T>(static_cast<range_type>(r) + static_cast<range_type>(min_value));
    }
}

/// Adaptor turning the engine into the callable a shuffle wants: g(n) returns a
/// value in [0, n).
class ShuffleRng {
public:
    explicit ShuffleRng(Rand48& g) : g_(&g) {}
    std::ptrdiff_t operator()(std::ptrdiff_t n) {
        return generate_uniform_int<Rand48, std::ptrdiff_t>(*g_, 0, n - 1);
    }
private:
    Rand48* g_;
};

/// Fisher-Yates shuffle, written out so the number of draws and the swap order
/// are fixed: element k (k = 1..n-1) is swapped with a uniformly chosen element
/// of [0, k]. std::shuffle is not a substitute - it consumes the engine in an
/// unspecified way, so it produces a different permutation per implementation.
template <class RandomAccessIterator, class Rng>
void random_shuffle(RandomAccessIterator begin, RandomAccessIterator end, Rng& rng) {
    if (begin == end) return;
    for (RandomAccessIterator it = begin + 1; it != end; ++it) {
        // rng(N) yields a value in [0, N), hence the +1 to include `it` itself.
        std::iter_swap(it, begin + rng((it - begin) + 1));
    }
}

} // namespace detail
} // namespace del3d

#endif // DEL3D_RAND48_H
