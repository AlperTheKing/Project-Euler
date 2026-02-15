#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
constexpr u64 kMod = 1000000000000000000ULL;  // last 18 digits

struct Options {
    int k = 10000;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = std::stoi(tail);
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
        if (parse_int_after_prefix(arg, "--k=", options.k)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.k >= 1;
}

inline u64 add_mod(u64 a, u64 b) {
    u64 s = a + b;
    if (s >= kMod) {
        s -= kMod;
    }
    return s;
}

struct FlatMap {
    std::vector<std::uint32_t> keys;
    std::vector<u64> values;
    std::size_t used = 0;

    static constexpr std::uint32_t kEmpty = 0xFFFFFFFFu;

    static std::uint32_t hash(std::uint32_t x) {
        x ^= x >> 16;
        x *= 0x7feb352d;
        x ^= x >> 15;
        x *= 0x846ca68b;
        x ^= x >> 16;
        return x;
    }

    void rehash(std::size_t cap) {
        std::vector<std::uint32_t> old_keys = std::move(keys);
        std::vector<u64> old_values = std::move(values);
        keys.assign(cap, kEmpty);
        values.assign(cap, 0);
        used = 0;
        for (std::size_t i = 0; i < old_keys.size(); ++i) {
            if (old_keys[i] != kEmpty) {
                insert(old_keys[i], old_values[i]);
            }
        }
    }

    bool get(std::uint32_t key, u64& out) const {
        if (keys.empty()) {
            return false;
        }
        const std::size_t mask = keys.size() - 1;
        std::size_t idx = hash(key) & mask;
        while (true) {
            const std::uint32_t cur = keys[idx];
            if (cur == kEmpty) {
                return false;
            }
            if (cur == key) {
                out = values[idx];
                return true;
            }
            idx = (idx + 1) & mask;
        }
    }

    void insert(std::uint32_t key, u64 value) {
        if (keys.empty() || (used + 1) * 10 >= keys.size() * 7) {
            const std::size_t next = keys.empty() ? 64 : keys.size() * 2;
            rehash(next);
        }
        const std::size_t mask = keys.size() - 1;
        std::size_t idx = hash(key) & mask;
        while (true) {
            const std::uint32_t cur = keys[idx];
            if (cur == kEmpty) {
                keys[idx] = key;
                values[idx] = value;
                ++used;
                return;
            }
            if (cur == key) {
                values[idx] = value;
                return;
            }
            idx = (idx + 1) & mask;
        }
    }
};

struct Solver400 {
    int n = 0;
    std::vector<int> sg;
    std::vector<FlatMap> cache;
    u64 nodes = 0;

    explicit Solver400(int n) : n(n), sg(n + 1, 0), cache(n + 1) {
        if (n >= 1) {
            sg[1] = 0;
        }
        if (n >= 2) {
            sg[2] = 1;
        }
        for (int i = 3; i <= n; ++i) {
            sg[i] = (1 + sg[i - 1]) ^ (1 + sg[i - 2]);
        }
    }

    u64 prune(int tree_size, int target) {
        if (target == -1) {
            return 1;
        }
        if (target < -1) {
            return 0;
        }
        if (tree_size == 1) {
            return 0;
        }
        if (tree_size == 2) {
            return target == 0 ? 1 : 0;
        }

        FlatMap& memo = cache[tree_size];
        u64 cached = 0;
        if (memo.get(static_cast<std::uint32_t>(target), cached)) {
            return cached;
        }

        ++nodes;
        const int right_target = ((sg[tree_size - 1] + 1) ^ target) - 1;
        const int left_target = ((sg[tree_size - 2] + 1) ^ target) - 1;
        const u64 right = prune(tree_size - 2, right_target);
        const u64 left = prune(tree_size - 1, left_target);
        u64 res = right + left;
        if (res >= kMod) {
            res -= kMod;
        }
        memo.insert(static_cast<std::uint32_t>(target), res);
        return res;
    }
};

u64 f_mod(const int target_k) {
    Solver400 solver(target_k);
    return solver.prune(target_k, 0);
}

bool run_checkpoints() {
    if (f_mod(5) != 1ULL) {
        std::cerr << "Checkpoint failed: f(5)\n";
        return false;
    }
    if (f_mod(10) != 17ULL) {
        std::cerr << "Checkpoint failed: f(10)\n";
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

    const u64 answer = f_mod(options.k);
    std::cout << std::setfill('0') << std::setw(18) << answer << '\n';
    return 0;
}
