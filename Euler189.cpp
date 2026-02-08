#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int size = 8;
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
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
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
        if (parse_int_after_prefix(arg, "--size=", options.size)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.size >= 1;
}

struct RowState {
    int odd_code = 0;
    int even_code = 0;
};

std::vector<RowState> build_row_states(const int length) {
    std::vector<RowState> states;
    std::vector<int> colors(static_cast<std::size_t>(length), 0);

    const auto dfs = [&](auto&& self, int idx) -> void {
        if (idx == length) {
            int odd_code = 0;
            int even_code = 0;
            int odd_pow = 1;
            int even_pow = 1;
            for (int i = 0; i < length; ++i) {
                if ((i & 1) == 0) {
                    even_code += colors[static_cast<std::size_t>(i)] * even_pow;
                    even_pow *= 3;
                } else {
                    odd_code += colors[static_cast<std::size_t>(i)] * odd_pow;
                    odd_pow *= 3;
                }
            }
            states.push_back({odd_code, even_code});
            return;
        }

        for (int c = 0; c < 3; ++c) {
            if (idx > 0 && colors[static_cast<std::size_t>(idx - 1)] == c) {
                continue;
            }
            colors[static_cast<std::size_t>(idx)] = c;
            self(self, idx + 1);
        }
    };
    dfs(dfs, 0);
    return states;
}

u64 solve(const int size) {
    std::vector<std::vector<RowState>> states_by_row(static_cast<std::size_t>(size));
    for (int r = 0; r < size; ++r) {
        states_by_row[static_cast<std::size_t>(r)] = build_row_states(2 * r + 1);
    }

    std::vector<u64> dp_prev(states_by_row[0].size(), 1ULL);
    for (int r = 1; r < size; ++r) {
        const int sig_count = [&]() {
            int v = 1;
            for (int i = 0; i < r; ++i) {
                v *= 3;
            }
            return v;
        }();

        std::vector<u64> agg(static_cast<std::size_t>(sig_count), 0ULL);
        for (std::size_t i = 0; i < dp_prev.size(); ++i) {
            const int sig = states_by_row[static_cast<std::size_t>(r - 1)][i].even_code;
            agg[static_cast<std::size_t>(sig)] += dp_prev[i];
        }

        std::vector<u64> memo(static_cast<std::size_t>(sig_count), std::numeric_limits<u64>::max());

        const auto compatible_sum = [&](int req_code) -> u64 {
            u64& cache = memo[static_cast<std::size_t>(req_code)];
            if (cache != std::numeric_limits<u64>::max()) {
                return cache;
            }

            std::array<int, 8> req_digits{};
            int tmp = req_code;
            for (int i = 0; i < r; ++i) {
                req_digits[static_cast<std::size_t>(i)] = tmp % 3;
                tmp /= 3;
            }

            u64 total = 0;
            const auto rec = [&](auto&& self, int idx, int code, int pow3) -> void {
                if (idx == r) {
                    total += agg[static_cast<std::size_t>(code)];
                    return;
                }
                for (int color = 0; color < 3; ++color) {
                    if (color == req_digits[static_cast<std::size_t>(idx)]) {
                        continue;
                    }
                    self(self, idx + 1, code + color * pow3, pow3 * 3);
                }
            };
            rec(rec, 0, 0, 1);
            cache = total;
            return total;
        };

        std::vector<u64> dp_curr(states_by_row[static_cast<std::size_t>(r)].size(), 0ULL);
        for (std::size_t i = 0; i < dp_curr.size(); ++i) {
            const int req = states_by_row[static_cast<std::size_t>(r)][i].odd_code;
            dp_curr[i] = compatible_sum(req);
        }

        dp_prev.swap(dp_curr);
    }

    u64 answer = 0;
    for (u64 value : dp_prev) {
        answer += value;
    }
    return answer;
}

u64 brute_small(const int size) {
    const auto id_of = [size](const int r, const int c) {
        return r * r + c;
    };
    const int nodes = size * size;
    std::vector<std::vector<int>> graph(static_cast<std::size_t>(nodes));

    for (int r = 0; r < size; ++r) {
        const int width = 2 * r + 1;
        for (int c = 0; c < width; ++c) {
            const int a = id_of(r, c);
            if (c + 1 < width) {
                const int b = id_of(r, c + 1);
                graph[static_cast<std::size_t>(a)].push_back(b);
                graph[static_cast<std::size_t>(b)].push_back(a);
            }
            if (r + 1 < size && (c % 2 == 0)) {
                const int b = id_of(r + 1, c + 1);
                graph[static_cast<std::size_t>(a)].push_back(b);
                graph[static_cast<std::size_t>(b)].push_back(a);
            }
        }
    }

    std::vector<int> color(static_cast<std::size_t>(nodes), -1);
    u64 answer = 0;

    const auto dfs = [&](auto&& self, int idx) -> void {
        if (idx == nodes) {
            ++answer;
            return;
        }
        for (int c = 0; c < 3; ++c) {
            bool ok = true;
            for (int nb : graph[static_cast<std::size_t>(idx)]) {
                if (nb < idx && color[static_cast<std::size_t>(nb)] == c) {
                    ok = false;
                    break;
                }
            }
            if (!ok) {
                continue;
            }
            color[static_cast<std::size_t>(idx)] = c;
            self(self, idx + 1);
            color[static_cast<std::size_t>(idx)] = -1;
        }
    };
    dfs(dfs, 0);
    return answer;
}

bool run_checkpoints() {
    for (int n = 1; n <= 4; ++n) {
        if (solve(n) != brute_small(n)) {
            std::cerr << "Checkpoint failed for brute cross-check at size " << n << '\n';
            return false;
        }
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
    std::cout << solve(options.size) << '\n';
    return 0;
}
