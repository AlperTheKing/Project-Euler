#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <functional>

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kTargetN = 100'000;
constexpr int kMaxMembership = 8; // Empirically 8 for N=100000.

struct Sequence {
    u64 s2;   // s^2 = a^2 + b^2
    int t;    // minimum odd |k|
    int kmax; // maximum odd |k| such that radius <= N
};

struct Occurrence {
    u64 radius2;
    std::uint16_t sequence_id;

    bool operator<(const Occurrence& other) const {
        if (radius2 != other.radius2) {
            return radius2 < other.radius2;
        }
        return sequence_id < other.sequence_id;
    }
};

struct IdSetKey {
    std::array<std::uint16_t, kMaxMembership> ids{};
    std::uint8_t len = 0;

    bool operator==(const IdSetKey& other) const {
        if (len != other.len) {
            return false;
        }
        for (std::uint8_t i = 0; i < len; ++i) {
            if (ids[i] != other.ids[i]) {
                return false;
            }
        }
        return true;
    }
};

struct IdSetKeyHash {
    std::size_t operator()(const IdSetKey& key) const {
        std::size_t h = 0x9e3779b97f4a7c15ULL ^ static_cast<std::size_t>(key.len);
        for (std::uint8_t i = 0; i < key.len; ++i) {
            const std::size_t v = static_cast<std::size_t>(key.ids[i]) +
                                  0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            h ^= v;
        }
        return h;
    }
};

inline u64 isqrt_u64(u64 n) {
    u64 x = static_cast<u64>(
        std::sqrt(static_cast<long double>(n)));
    while ((x + 1) <= n / (x + 1)) {
        ++x;
    }
    while (x > n / x) {
        --x;
    }
    return x;
}

struct EGResult {
    i64 x;
    i64 y;
    i64 g;
};

EGResult extended_gcd(i64 a, i64 b) {
    i64 old_r = a;
    i64 r = b;
    i64 old_s = 1;
    i64 s = 0;
    i64 old_t = 0;
    i64 t = 1;

    while (r != 0) {
        const i64 q = old_r / r;

        const i64 nr = old_r - q * r;
        old_r = r;
        r = nr;

        const i64 ns = old_s - q * s;
        old_s = s;
        s = ns;

        const i64 nt = old_t - q * t;
        old_t = t;
        t = nt;
    }

    return EGResult{old_s, old_t, old_r};
}

int compute_threshold_for_vector(int a, int b) {
    // For primitive odd (a,b), circles through two lattice intersection points are
    // indexed by odd k with radius^2 = s^2*(1+k^2)/4, s^2=a^2+b^2.
    // The no-interior-lattice-point condition for the lens formed by k=-m and k=n
    // is equivalent to m,n >= t(a,b), where t is the odd ceiling of
    // M = max(A - (A^2+1)/s^2), over A induced by -b*x + a*y = 1.

    const i64 aa = static_cast<i64>(a);
    const i64 bb = static_cast<i64>(b);
    const i64 s2 = aa * aa + bb * bb;

    const EGResult eg = extended_gcd(aa, bb); // aa*eg.x + bb*eg.y = 1

    const i64 x0 = -eg.y;
    const i64 y0 = eg.x; // -bb*x0 + aa*y0 = 1

    i64 A0 = aa * x0 + bb * y0;
    i64 r = A0 % s2;
    if (r < 0) {
        r += s2;
    }
    if (r > s2 / 2) {
        r = s2 - r;
    }

    const i64 quad = static_cast<i64>((static_cast<u128>(r) * r + 1) / s2);
    const i64 M = r - quad;

    i64 t = (M & 1LL) ? M : (M + 1);
    if (t < 1) {
        t = 1;
    }

    return static_cast<int>(t);
}

std::vector<Sequence> generate_sequences(u64 n_limit) {
    const u64 n2 = n_limit * n_limit;
    const int b_max = static_cast<int>(isqrt_u64(2 * n_limit));

    std::vector<Sequence> sequences;
    sequences.reserve(2048);

    std::unordered_set<u64> seen;
    seen.reserve(4096);

    for (int a = 1; a <= b_max; a += 2) {
        for (int b = a; b <= b_max; b += 2) {
            if (std::gcd(a, b) != 1) {
                continue;
            }

            const u64 s2 = static_cast<u64>(a) * a + static_cast<u64>(b) * b;
            const int t = compute_threshold_for_vector(a, b);

            const u64 ratio = (4 * n2) / s2;
            if (ratio <= 1) {
                continue;
            }

            int kmax = static_cast<int>(isqrt_u64(ratio - 1));
            if ((kmax & 1) == 0) {
                --kmax;
            }
            if (kmax < t) {
                continue;
            }

            const u64 key = (static_cast<u64>(s2) << 32) |
                            static_cast<u64>(static_cast<std::uint32_t>(t));
            if (!seen.insert(key).second) {
                continue;
            }

            sequences.push_back(Sequence{s2, t, kmax});
        }
    }

    std::sort(sequences.begin(), sequences.end(), [](const Sequence& lhs, const Sequence& rhs) {
        if (lhs.s2 != rhs.s2) {
            return lhs.s2 < rhs.s2;
        }
        return lhs.t < rhs.t;
    });

    return sequences;
}

u64 total_occurrences(const std::vector<Sequence>& sequences) {
    u64 total = 0;
    for (const Sequence& seq : sequences) {
        total += static_cast<u64>((seq.kmax - seq.t) / 2 + 1);
    }
    return total;
}

std::vector<Occurrence> build_occurrences(const std::vector<Sequence>& sequences,
                                          bool allow_multithreading,
                                          unsigned requested_threads) {
    const u64 total = total_occurrences(sequences);

    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 1;
        }
    }

    if (!allow_multithreading || sequences.size() < 128 || total < 200'000) {
        threads = 1;
    }

    if (threads > sequences.size()) {
        threads = static_cast<unsigned>(sequences.size());
    }
    if (threads == 0) {
        threads = 1;
    }

    std::vector<Occurrence> merged;
    merged.reserve(static_cast<std::size_t>(total));

    if (threads == 1) {
        for (std::size_t sid = 0; sid < sequences.size(); ++sid) {
            const Sequence& seq = sequences[sid];
            for (int k = seq.t; k <= seq.kmax; k += 2) {
                const u64 k2 = static_cast<u64>(k) * k;
                const u64 radius2 = static_cast<u64>(
                    (static_cast<u128>(seq.s2) * (1 + k2)) / 4);
                merged.push_back(Occurrence{radius2, static_cast<std::uint16_t>(sid)});
            }
        }
        return merged;
    }

    std::atomic<std::size_t> next_index{0};
    std::vector<std::vector<Occurrence>> locals(threads);
    for (auto& local : locals) {
        local.reserve(static_cast<std::size_t>(total / threads) + 1024);
    }

    std::vector<std::thread> pool;
    pool.reserve(threads);

    for (unsigned tid = 0; tid < threads; ++tid) {
        pool.emplace_back([&, tid]() {
            std::vector<Occurrence>& out = locals[tid];
            while (true) {
                const std::size_t sid = next_index.fetch_add(1, std::memory_order_relaxed);
                if (sid >= sequences.size()) {
                    break;
                }

                const Sequence& seq = sequences[sid];
                for (int k = seq.t; k <= seq.kmax; k += 2) {
                    const u64 k2 = static_cast<u64>(k) * k;
                    const u64 radius2 = static_cast<u64>(
                        (static_cast<u128>(seq.s2) * (1 + k2)) / 4);
                    out.push_back(Occurrence{radius2, static_cast<std::uint16_t>(sid)});
                }
            }
        });
    }

    for (std::thread& worker : pool) {
        worker.join();
    }

    for (auto& local : locals) {
        merged.insert(merged.end(), local.begin(), local.end());
    }

    return merged;
}

u64 count_distinct_lenticular_pairs(u64 n_limit,
                                    bool allow_multithreading,
                                    unsigned requested_threads = 0) {
    const std::vector<Sequence> sequences = generate_sequences(n_limit);
    if (sequences.size() > std::numeric_limits<std::uint16_t>::max()) {
        std::cerr << "Too many sequence classes (" << sequences.size()
                  << "). Increase sequence_id width." << '\n';
        std::exit(1);
    }

    std::vector<Occurrence> occurrences =
        build_occurrences(sequences, allow_multithreading, requested_threads);
    std::sort(occurrences.begin(), occurrences.end());

    std::unordered_map<IdSetKey, u32, IdSetKeyHash> full_set_frequency;
    full_set_frequency.reserve(8192);

    u64 distinct_radii = 0;

    std::size_t pos = 0;
    while (pos < occurrences.size()) {
        const u64 current_radius = occurrences[pos].radius2;
        IdSetKey key;

        while (pos < occurrences.size() && occurrences[pos].radius2 == current_radius) {
            const std::uint16_t sid = occurrences[pos].sequence_id;
            if (key.len == 0 || sid != key.ids[key.len - 1]) {
                if (key.len >= kMaxMembership) {
                    std::cerr << "Membership overflow for radius^2=" << current_radius
                              << ". Increase kMaxMembership." << '\n';
                    std::exit(1);
                }
                key.ids[key.len++] = sid;
            }
            ++pos;
        }

        ++distinct_radii;
        ++full_set_frequency[key];
    }

    std::unordered_map<IdSetKey, u64, IdSetKeyHash> subset_frequency;
    subset_frequency.reserve(131'072);

    for (const auto& entry : full_set_frequency) {
        const IdSetKey& full = entry.first;
        const u64 freq = entry.second;

        const int subset_count = 1 << full.len;
        for (int mask = 1; mask < subset_count; ++mask) {
            IdSetKey subset;
            for (std::uint8_t bit = 0; bit < full.len; ++bit) {
                if ((mask >> bit) & 1) {
                    subset.ids[subset.len++] = full.ids[bit];
                }
            }
            subset_frequency[subset] += freq;
        }
    }

    i64 intersecting_distinct_pairs = 0;

    for (const auto& entry : subset_frequency) {
        const IdSetKey& subset = entry.first;
        const u64 cnt = entry.second;
        if (cnt < 2) {
            continue;
        }

        const i64 ways = static_cast<i64>(cnt * (cnt - 1) / 2);
        if (subset.len & 1U) {
            intersecting_distinct_pairs += ways;
        } else {
            intersecting_distinct_pairs -= ways;
        }
    }

    // Include (r, r) pairs: every realized radius can be paired with itself.
    const u64 answer = distinct_radii + static_cast<u64>(intersecting_distinct_pairs);
    return answer;
}

bool run_validation_checkpoints() {
    struct Checkpoint {
        u64 n_limit;
        u64 expected;
    };

    const std::vector<Checkpoint> checkpoints = {
        {10, 30},
        {100, 3442},
    };

    for (const Checkpoint cp : checkpoints) {
        const u64 got = count_distinct_lenticular_pairs(cp.n_limit, false, 1);
        if (got != cp.expected) {
            std::cerr << "Checkpoint failed for N=" << cp.n_limit << ": got " << got
                      << ", expected " << cp.expected << '\n';
            return false;
        }
    }

    return true;
}

} // namespace

int main() {
    if (!run_validation_checkpoints()) {
        return 1;
    }

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 4;
    }

    const u64 answer = count_distinct_lenticular_pairs(kTargetN, true, threads);
    std::cout << answer << '\n';
    return 0;
}
