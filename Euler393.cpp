#include <cstdint>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using BigInt = boost::multiprecision::cpp_int;

struct Options {
    int n = 10;
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
    int parsed = 0;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1 && options.n <= 12;
}

enum class MoveType : std::uint8_t { Filled, Horizontal, Vertical };

struct Move {
    MoveType type;
    int occ_add;
    int out_add;
};

class Solver {
   public:
    explicit Solver(const int n) : n_(n), bit_mask_((1 << n) - 1) {}

    BigInt solve() {
        std::unordered_map<int, BigInt> dp;
        dp.reserve(1 << 12);
        dp.emplace(0, BigInt(1));

        for (int row = 0; row < n_; ++row) {
            const bool last_row = (row == n_ - 1);
            std::unordered_map<int, BigInt> next;
            next.reserve(dp.size() * 4 + 16);

            for (const auto& entry : dp) {
                const int state = entry.first;
                const BigInt& ways = entry.second;
                const int in1 = state & bit_mask_;
                const int in2 = (state >> n_) & bit_mask_;

                const auto& transitions = get_transitions(in1, in2, last_row);
                for (const auto& tr : transitions) {
                    next[tr.first] += ways * tr.second;
                }
            }

            dp.swap(next);
        }

        const auto it = dp.find(0);
        if (it == dp.end()) {
            return BigInt(0);
        }
        return it->second;
    }

   private:
    int n_;
    int bit_mask_;
    std::unordered_map<std::uint32_t, std::vector<std::pair<int, std::uint32_t>>> transition_cache_;

    std::vector<Move> options_for_cell(const int c, const int occ, const bool last_row) const {
        std::vector<Move> options;
        const int bit = 1 << c;
        if ((occ & bit) != 0) {
            options.push_back(Move{MoveType::Filled, 0, 0});
            return options;
        }

        if (c + 1 < n_ && (occ & (1 << (c + 1))) == 0) {
            options.push_back(Move{MoveType::Horizontal, bit | (1 << (c + 1)), 0});
        }
        if (!last_row) {
            options.push_back(Move{MoveType::Vertical, bit, bit});
        }

        return options;
    }

    const std::vector<std::pair<int, std::uint32_t>>& get_transitions(
        const int in1, const int in2, const bool last_row) {
        const std::uint32_t key = static_cast<std::uint32_t>(in1 | (in2 << n_) | (static_cast<int>(last_row) << 24));
        const auto cache_it = transition_cache_.find(key);
        if (cache_it != transition_cache_.end()) {
            return cache_it->second;
        }

        std::unordered_map<int, std::uint32_t> accum;
        accum.reserve(128);

        std::function<void(int, int, int, int, int)> dfs =
            [&](const int c, const int occ1, const int occ2, const int out1, const int out2) {
                if (c == n_) {
                    const int out_state = out1 | (out2 << n_);
                    ++accum[out_state];
                    return;
                }

                const std::vector<Move> ops1 = options_for_cell(c, occ1, last_row);
                const std::vector<Move> ops2 = options_for_cell(c, occ2, last_row);

                for (const Move& m1 : ops1) {
                    for (const Move& m2 : ops2) {
                        if (m1.type == MoveType::Horizontal && m2.type == MoveType::Horizontal) {
                            continue;  // Same undirected horizontal edge used by both matchings.
                        }
                        if (m1.type == MoveType::Vertical && m2.type == MoveType::Vertical) {
                            continue;  // Same undirected vertical edge used by both matchings.
                        }

                        const int next_occ1 = occ1 | m1.occ_add;
                        const int next_occ2 = occ2 | m2.occ_add;
                        const int next_out1 = out1 | m1.out_add;
                        const int next_out2 = out2 | m2.out_add;
                        dfs(c + 1, next_occ1, next_occ2, next_out1, next_out2);
                    }
                }
            };

        dfs(0, in1, in2, 0, 0);

        std::vector<std::pair<int, std::uint32_t>> packed;
        packed.reserve(accum.size());
        for (const auto& entry : accum) {
            packed.push_back(entry);
        }

        auto [it, inserted] = transition_cache_.emplace(key, std::move(packed));
        (void)inserted;
        return it->second;
    }
};

bool is_valid_assignment(const int n, const std::vector<int>& to) {
    const int total = n * n;
    std::vector<int> indeg(static_cast<std::size_t>(total), 0);
    auto rc = [n](const int id) {
        return std::pair<int, int>{id / n, id % n};
    };

    for (int u = 0; u < total; ++u) {
        const int v = to[static_cast<std::size_t>(u)];
        if (v < 0 || v >= total) {
            return false;
        }
        const auto [r1, c1] = rc(u);
        const auto [r2, c2] = rc(v);
        const int manhattan = std::abs(r1 - r2) + std::abs(c1 - c2);
        if (manhattan != 1) {
            return false;
        }
        ++indeg[static_cast<std::size_t>(v)];
        if (indeg[static_cast<std::size_t>(v)] > 1) {
            return false;
        }
    }

    for (int u = 0; u < total; ++u) {
        const int v = to[static_cast<std::size_t>(u)];
        if (to[static_cast<std::size_t>(v)] == u) {
            return false;  // Crossing the same edge in opposite directions.
        }
    }

    for (int v = 0; v < total; ++v) {
        if (indeg[static_cast<std::size_t>(v)] != 1) {
            return false;
        }
    }

    return true;
}

std::uint64_t brute_force_small(const int n) {
    const int total = n * n;
    std::vector<std::vector<int>> neighbors(static_cast<std::size_t>(total));
    auto id = [n](const int r, const int c) {
        return r * n + c;
    };

    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            const int u = id(r, c);
            if (r > 0) neighbors[static_cast<std::size_t>(u)].push_back(id(r - 1, c));
            if (r + 1 < n) neighbors[static_cast<std::size_t>(u)].push_back(id(r + 1, c));
            if (c > 0) neighbors[static_cast<std::size_t>(u)].push_back(id(r, c - 1));
            if (c + 1 < n) neighbors[static_cast<std::size_t>(u)].push_back(id(r, c + 1));
        }
    }

    std::vector<int> to(static_cast<std::size_t>(total), -1);
    std::uint64_t count = 0;

    std::function<void(int)> dfs = [&](const int u) {
        if (u == total) {
            if (is_valid_assignment(n, to)) {
                ++count;
            }
            return;
        }
        for (const int v : neighbors[static_cast<std::size_t>(u)]) {
            to[static_cast<std::size_t>(u)] = v;
            dfs(u + 1);
        }
    };

    dfs(0);
    return count;
}

bool run_checkpoints() {
    Solver s2(2);
    const BigInt f2 = s2.solve();
    if (f2 != 2) {
        std::cerr << "Checkpoint failed: f(2)\n";
        return false;
    }

    Solver s4(4);
    const BigInt f4 = s4.solve();
    if (f4 != 88) {
        std::cerr << "Checkpoint failed: f(4)\n";
        return false;
    }

    const std::uint64_t brute3 = brute_force_small(3);
    Solver s3(3);
    const BigInt f3 = s3.solve();
    if (f3 != brute3) {
        std::cerr << "Checkpoint failed: brute-force mismatch for n=3\n";
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

    Solver solver(options.n);
    std::cout << solver.solve() << '\n';
    return 0;
}
