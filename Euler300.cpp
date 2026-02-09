#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    int n = 15;
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
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 2 && options.n <= 20;
}

struct ContactKey {
    u64 lo;
    u64 hi;

    bool operator==(const ContactKey& other) const {
        return lo == other.lo && hi == other.hi;
    }
};

struct ContactKeyHash {
    std::size_t operator()(const ContactKey& k) const {
        const u64 x = k.lo * 11400714819323198485ull;
        const u64 y = k.hi * 14029467366897019727ull;
        return static_cast<std::size_t>(x ^ (y >> 1));
    }
};

struct FoldingSolver {
    int n;
    int grid_side;
    int offset;

    std::vector<int> px;
    std::vector<int> py;
    std::vector<unsigned char> occupied;

    std::array<int, 4> dx{1, -1, 0, 0};
    std::array<int, 4> dy{0, 0, 1, -1};

    std::vector<std::pair<int, int>> pair_list;
    std::vector<std::vector<int>> pair_index;

    std::unordered_set<ContactKey, ContactKeyHash> contact_maps;

    explicit FoldingSolver(int n_) : n(n_) {
        grid_side = 2 * n + 1;
        offset = n;
        px.assign(static_cast<std::size_t>(n), 0);
        py.assign(static_cast<std::size_t>(n), 0);
        occupied.assign(static_cast<std::size_t>(grid_side * grid_side), 0);

        pair_index.assign(static_cast<std::size_t>(n), std::vector<int>(static_cast<std::size_t>(n), -1));
        int idx = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 2; j < n; ++j) {
                pair_index[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = idx;
                pair_list.push_back({i, j});
                ++idx;
            }
        }

        contact_maps.reserve(1 << 16);
    }

    int cell_id(int x, int y) const {
        return y * grid_side + x;
    }

    void add_current_contact_map() {
        ContactKey key{0, 0};

        for (int p = 0; p < static_cast<int>(pair_list.size()); ++p) {
            const int i = pair_list[static_cast<std::size_t>(p)].first;
            const int j = pair_list[static_cast<std::size_t>(p)].second;
            const int manhattan = std::abs(px[static_cast<std::size_t>(i)] - px[static_cast<std::size_t>(j)]) +
                                  std::abs(py[static_cast<std::size_t>(i)] - py[static_cast<std::size_t>(j)]);
            if (manhattan != 1) {
                continue;
            }

            if (p < 64) {
                key.lo |= (1ULL << p);
            } else {
                key.hi |= (1ULL << (p - 64));
            }
        }

        contact_maps.insert(key);
    }

    void dfs(int idx) {
        if (idx == n - 1) {
            add_current_contact_map();
            return;
        }

        const int cx = px[static_cast<std::size_t>(idx)];
        const int cy = py[static_cast<std::size_t>(idx)];

        for (int dir = 0; dir < 4; ++dir) {
            const int nx = cx + dx[static_cast<std::size_t>(dir)];
            const int ny = cy + dy[static_cast<std::size_t>(dir)];
            const int cid = cell_id(nx, ny);
            if (occupied[static_cast<std::size_t>(cid)] != 0) {
                continue;
            }
            occupied[static_cast<std::size_t>(cid)] = 1;
            px[static_cast<std::size_t>(idx + 1)] = nx;
            py[static_cast<std::size_t>(idx + 1)] = ny;
            dfs(idx + 1);
            occupied[static_cast<std::size_t>(cid)] = 0;
        }
    }

    std::vector<ContactKey> enumerate_maps() {
        const int x0 = offset;
        const int y0 = offset;
        const int x1 = offset + 1;
        const int y1 = offset;

        occupied[static_cast<std::size_t>(cell_id(x0, y0))] = 1;
        occupied[static_cast<std::size_t>(cell_id(x1, y1))] = 1;

        px[0] = x0;
        py[0] = y0;
        px[1] = x1;
        py[1] = y1;

        dfs(1);

        std::vector<ContactKey> maps;
        maps.reserve(contact_maps.size());
        for (const ContactKey& key : contact_maps) {
            maps.push_back(key);
        }
        return maps;
    }

    static std::string exact_average_decimal(u64 numerator, int n) {
        const u128 denominator = static_cast<u128>(1) << n;
        const u128 num = numerator;
        const u128 integer_part = num / denominator;
        u128 rem = num % denominator;

        std::string out = std::to_string(static_cast<u64>(integer_part));
        if (rem == 0) {
            return out;
        }

        out.push_back('.');
        while (rem != 0) {
            rem *= 10;
            const u64 digit = static_cast<u64>(rem / denominator);
            out.push_back(static_cast<char>('0' + digit));
            rem %= denominator;
        }
        return out;
    }

    u64 solve_numerator() {
        const std::vector<ContactKey> maps = enumerate_maps();

        std::vector<std::array<std::uint16_t, 20>> adj_masks;
        adj_masks.reserve(maps.size());

        for (const ContactKey& key : maps) {
            std::array<std::uint16_t, 20> adj{};
            adj.fill(0);

            for (int p = 0; p < static_cast<int>(pair_list.size()); ++p) {
                bool on = false;
                if (p < 64) {
                    on = ((key.lo >> p) & 1ULL) != 0ULL;
                } else {
                    on = ((key.hi >> (p - 64)) & 1ULL) != 0ULL;
                }
                if (!on) {
                    continue;
                }

                const int i = pair_list[static_cast<std::size_t>(p)].first;
                const int j = pair_list[static_cast<std::size_t>(p)].second;
                adj[static_cast<std::size_t>(i)] |= static_cast<std::uint16_t>(1U << j);
                adj[static_cast<std::size_t>(j)] |= static_cast<std::uint16_t>(1U << i);
            }

            adj_masks.push_back(adj);
        }

        const int mask_count = 1 << n;
        std::vector<std::uint8_t> best(static_cast<std::size_t>(mask_count), 0);
        std::vector<std::uint8_t> score(static_cast<std::size_t>(mask_count), 0);

        for (const auto& adj : adj_masks) {
            score[0] = 0;
            for (int mask = 1; mask < mask_count; ++mask) {
                const int lb = mask & -mask;
                const int i = __builtin_ctz(static_cast<unsigned>(lb));
                const int prev = mask ^ lb;
                const std::uint16_t m = adj[static_cast<std::size_t>(i)];
                const int add = __builtin_popcount(static_cast<unsigned>(prev & m));
                score[static_cast<std::size_t>(mask)] =
                    static_cast<std::uint8_t>(score[static_cast<std::size_t>(prev)] + add);

                if (score[static_cast<std::size_t>(mask)] > best[static_cast<std::size_t>(mask)]) {
                    best[static_cast<std::size_t>(mask)] = score[static_cast<std::size_t>(mask)];
                }
            }
        }

        u64 numerator = 0;
        for (int mask = 0; mask < mask_count; ++mask) {
            const int consecutive_hh = __builtin_popcount(static_cast<unsigned>(mask & (mask >> 1)));
            numerator += static_cast<u64>(best[static_cast<std::size_t>(mask)]) +
                         static_cast<u64>(consecutive_hh);
        }
        return numerator;
    }
};

std::string solve_average_decimal(int n) {
    FoldingSolver solver(n);
    const u64 numerator = solver.solve_numerator();
    return FoldingSolver::exact_average_decimal(numerator, n);
}

bool run_checkpoints() {
    const std::string avg8 = solve_average_decimal(8);
    if (avg8 != "3.3203125") {
        std::cerr << "Checkpoint failed for n=8: got " << avg8 << '\n';
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

    std::cout << solve_average_decimal(options.n) << '\n';
    return 0;
}
