#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
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

u64 f_mod(const int target_k) {
    if (target_k == 1) {
        return 0;
    }

    // For a component tree U, SG(U) = (SG(left) xor SG(right)) + 1.
    // Let h(k)=SG(T(k)) for component play (root removable).
    int h_k_minus_2 = 0;  // h(0)
    int h_k_minus_1 = 1;  // h(1)

    std::unordered_map<int, u64> d_k_minus_2;  // D(0): no moves
    std::unordered_map<int, u64> d_k_minus_1;  // D(1): remove root -> 0
    d_k_minus_1.emplace(0, 1ULL);
    d_k_minus_2.reserve(1024);
    d_k_minus_1.reserve(1024);

    std::unordered_map<int, u64> d_cur;
    d_cur.reserve(2048);

    u64 f_value = 0ULL;

    for (int k = 2; k <= target_k; ++k) {
        // Full poisoned-root game on T(k) is a sum of components T(k-1) and T(k-2).
        // Winning first moves are those that leave xor = 0.
        const auto it_a = d_k_minus_1.find(h_k_minus_2);
        const auto it_b = d_k_minus_2.find(h_k_minus_1);
        const u64 part_a = (it_a == d_k_minus_1.end() ? 0ULL : it_a->second);
        const u64 part_b = (it_b == d_k_minus_2.end() ? 0ULL : it_b->second);
        f_value = add_mod(part_a, part_b);

        const int h_k = (h_k_minus_1 ^ h_k_minus_2) + 1;

        d_cur.clear();
        d_cur.reserve(d_k_minus_1.size() + d_k_minus_2.size() + 8);

        auto add_count = [&](const int idx, const u64 val) {
            const auto it = d_cur.find(idx);
            if (it == d_cur.end()) {
                d_cur.emplace(idx, val % kMod);
            } else {
                it->second = add_mod(it->second, val);
            }
        };

        // D(k): move-result nimber counts for component T(k).
        // 1) remove root -> empty (nimber 0)
        add_count(0, 1ULL);
        // 2) move in left child T(k-1)
        for (const auto& entry : d_k_minus_1) {
            const int v = entry.first;
            const int out = (v ^ h_k_minus_2) + 1;
            add_count(out, entry.second);
        }
        // 3) move in right child T(k-2)
        for (const auto& entry : d_k_minus_2) {
            const int v = entry.first;
            const int out = (h_k_minus_1 ^ v) + 1;
            add_count(out, entry.second);
        }

        d_k_minus_2.swap(d_k_minus_1);
        d_k_minus_1.swap(d_cur);

        h_k_minus_2 = h_k_minus_1;
        h_k_minus_1 = h_k;
    }

    return f_value;
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
