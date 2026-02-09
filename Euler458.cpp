#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u32 kMod = 1'000'000'000U;
constexpr int kLen = 6;
constexpr int kAlphabet = 7;

struct Options {
    u64 n = 1'000'000'000'000ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = static_cast<u64>(std::stoull(tail));
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

u32 mod_pow(u64 base, u64 exp) {
    u64 result = 1ULL;
    base %= kMod;
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = (result * base) % kMod;
        }
        base = (base * base) % kMod;
        exp >>= 1ULL;
    }
    return static_cast<u32>(result);
}

using Pattern = std::array<std::uint8_t, kLen>;

u32 encode_pattern(const Pattern& p) {
    u32 key = 0U;
    for (int i = 0; i < kLen; ++i) {
        key |= static_cast<u32>(p[static_cast<std::size_t>(i)]) << (3 * i);
    }
    return key;
}

Pattern canonicalize(const Pattern& p) {
    std::array<int, kAlphabet + 1> remap{};
    remap.fill(-1);
    int next = 0;
    Pattern out{};
    for (int i = 0; i < kLen; ++i) {
        const int v = p[static_cast<std::size_t>(i)];
        if (remap[static_cast<std::size_t>(v)] == -1) {
            remap[static_cast<std::size_t>(v)] = next++;
        }
        out[static_cast<std::size_t>(i)] =
            static_cast<std::uint8_t>(remap[static_cast<std::size_t>(v)]);
    }
    return out;
}

void generate_patterns_recursive(int pos,
                                 int current_max,
                                 Pattern& current,
                                 std::vector<Pattern>& out) {
    if (pos == kLen) {
        out.push_back(current);
        return;
    }
    for (int v = 0; v <= current_max + 1; ++v) {
        current[static_cast<std::size_t>(pos)] = static_cast<std::uint8_t>(v);
        generate_patterns_recursive(pos + 1, (v > current_max ? v : current_max), current, out);
    }
}

using Matrix = std::vector<std::vector<u32>>;

Matrix multiply_matrix(const Matrix& a, const Matrix& b) {
    const std::size_t n = a.size();
    Matrix c(n, std::vector<u32>(n, 0U));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            u128 sum = 0;
            for (std::size_t k = 0; k < n; ++k) {
                sum += static_cast<u128>(a[i][k]) * b[k][j];
            }
            c[i][j] = static_cast<u32>(sum % kMod);
        }
    }
    return c;
}

std::vector<u32> multiply_mat_vec(const Matrix& a, const std::vector<u32>& v) {
    const std::size_t n = a.size();
    std::vector<u32> out(n, 0U);
    for (std::size_t i = 0; i < n; ++i) {
        u128 sum = 0;
        for (std::size_t j = 0; j < n; ++j) {
            sum += static_cast<u128>(a[i][j]) * v[j];
        }
        out[i] = static_cast<u32>(sum % kMod);
    }
    return out;
}

u32 solve(const u64 n) {
    if (n == 0ULL) {
        return 1U;
    }
    if (n <= 6ULL) {
        return mod_pow(7U, n);
    }

    std::vector<Pattern> patterns;
    patterns.reserve(203);
    Pattern seed{};
    seed.fill(0U);
    seed[0] = 0U;
    generate_patterns_recursive(1, 0, seed, patterns);

    std::unordered_map<u32, u32> id_by_key;
    id_by_key.reserve(patterns.size() * 2U);
    for (u32 id = 0U; id < static_cast<u32>(patterns.size()); ++id) {
        id_by_key.emplace(encode_pattern(patterns[static_cast<std::size_t>(id)]), id);
    }

    const std::size_t state_count = patterns.size();
    Matrix trans(state_count, std::vector<u32>(state_count, 0U));

    for (u32 id = 0U; id < static_cast<u32>(state_count); ++id) {
        const Pattern& p = patterns[static_cast<std::size_t>(id)];
        int k = 0;
        for (int i = 0; i < kLen; ++i) {
            const int label = p[static_cast<std::size_t>(i)];
            if (label + 1 > k) {
                k = label + 1;
            }
        }

        for (int label = 0; label < k; ++label) {
            Pattern next{};
            for (int t = 0; t < kLen - 1; ++t) {
                next[static_cast<std::size_t>(t)] = p[static_cast<std::size_t>(t + 1)];
            }
            next[static_cast<std::size_t>(kLen - 1)] = static_cast<std::uint8_t>(label);
            next = canonicalize(next);
            const u32 nid = id_by_key[encode_pattern(next)];
            trans[static_cast<std::size_t>(nid)][static_cast<std::size_t>(id)] =
                (trans[static_cast<std::size_t>(nid)][static_cast<std::size_t>(id)] + 1U) % kMod;
        }

        const int new_count = kAlphabet - k;
        if (new_count > 0 && k < 6) {
            Pattern next{};
            for (int t = 0; t < kLen - 1; ++t) {
                next[static_cast<std::size_t>(t)] = p[static_cast<std::size_t>(t + 1)];
            }
            next[static_cast<std::size_t>(kLen - 1)] = static_cast<std::uint8_t>(k);
            next = canonicalize(next);
            const u32 nid = id_by_key[encode_pattern(next)];
            trans[static_cast<std::size_t>(nid)][static_cast<std::size_t>(id)] =
                (trans[static_cast<std::size_t>(nid)][static_cast<std::size_t>(id)] +
                 static_cast<u32>(new_count)) %
                kMod;
        }
    }

    std::vector<u32> vec(state_count, 0U);
    for (u32 id = 0U; id < static_cast<u32>(state_count); ++id) {
        const Pattern& p = patterns[static_cast<std::size_t>(id)];
        int k = 0;
        for (int i = 0; i < kLen; ++i) {
            const int label = p[static_cast<std::size_t>(i)];
            if (label + 1 > k) {
                k = label + 1;
            }
        }

        u32 ways = 1U;
        for (int t = 0; t < k; ++t) {
            ways = static_cast<u32>((static_cast<u64>(ways) * (kAlphabet - t)) % kMod);
        }
        vec[static_cast<std::size_t>(id)] = ways;
    }

    u64 steps = n - 6ULL;
    Matrix power = trans;
    while (steps > 0ULL) {
        if ((steps & 1ULL) != 0ULL) {
            vec = multiply_mat_vec(power, vec);
        }
        steps >>= 1ULL;
        if (steps > 0ULL) {
            power = multiply_matrix(power, power);
        }
    }

    u64 answer = 0ULL;
    for (const u32 v : vec) {
        answer += v;
        answer %= kMod;
    }
    return static_cast<u32>(answer);
}

bool run_checkpoints() {
    if (solve(7ULL) != 818'503U) {
        std::cerr << "Checkpoint failed: T(7)\n";
        return false;
    }
    if (solve(8ULL) != 5'699'281U) {
        std::cerr << "Checkpoint failed: T(8)\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    std::cout << solve(options.n) << '\n';
    return 0;
}
