#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>

namespace {

using i64 = std::int64_t;
using ld = long double;

constexpr int kSize = 5;
constexpr int kCells = kSize * kSize;
constexpr int kCols = 5;
constexpr int kMaskCount = 1 << kCols;
constexpr int kFullMask = kMaskCount - 1;
constexpr int kRhsCount = 1 + kCols;  // expected time + hit probabilities per column

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

int make_pos(const int x, const int y) {
    return y * kSize + x;
}

int top_pos(const int col) {
    return make_pos(col, 0);
}

int bottom_pos(const int col) {
    return make_pos(col, kSize - 1);
}

struct Grid {
    std::array<std::array<int, 4>, kCells> neighbors{};
    std::array<int, kCells> degree{};
};

Grid build_grid() {
    Grid grid;
    for (int y = 0; y < kSize; ++y) {
        for (int x = 0; x < kSize; ++x) {
            const int p = make_pos(x, y);
            int d = 0;
            if (x > 0) {
                grid.neighbors[p][d++] = make_pos(x - 1, y);
            }
            if (x + 1 < kSize) {
                grid.neighbors[p][d++] = make_pos(x + 1, y);
            }
            if (y > 0) {
                grid.neighbors[p][d++] = make_pos(x, y - 1);
            }
            if (y + 1 < kSize) {
                grid.neighbors[p][d++] = make_pos(x, y + 1);
            }
            grid.degree[p] = d;
        }
    }
    return grid;
}

bool solve_linear_system(
    std::array<std::array<ld, kCells>, kCells>& a,
    std::array<std::array<ld, kRhsCount>, kCells>& rhs
) {
    constexpr ld kEps = 1e-18L;

    for (int col = 0; col < kCells; ++col) {
        int pivot = col;
        ld best = std::fabsl(a[col][col]);
        for (int row = col + 1; row < kCells; ++row) {
            const ld cand = std::fabsl(a[row][col]);
            if (cand > best) {
                best = cand;
                pivot = row;
            }
        }
        if (best < kEps) {
            return false;
        }

        if (pivot != col) {
            std::swap(a[pivot], a[col]);
            std::swap(rhs[pivot], rhs[col]);
        }

        const ld inv_pivot = 1.0L / a[col][col];
        for (int j = col; j < kCells; ++j) {
            a[col][j] *= inv_pivot;
        }
        for (int r = 0; r < kRhsCount; ++r) {
            rhs[col][r] *= inv_pivot;
        }

        for (int row = 0; row < kCells; ++row) {
            if (row == col) {
                continue;
            }
            const ld factor = a[row][col];
            if (std::fabsl(factor) < kEps) {
                continue;
            }
            for (int j = col; j < kCells; ++j) {
                a[row][j] -= factor * a[col][j];
            }
            for (int r = 0; r < kRhsCount; ++r) {
                rhs[row][r] -= factor * rhs[col][r];
            }
        }
    }

    return true;
}

struct HittingData {
    std::array<ld, kCells> expect{};
    std::array<std::array<ld, kCols>, kCells> prob{};
};

HittingData solve_hitting_data(const Grid& grid, const int row, const int mask, bool& ok) {
    HittingData data;

    std::array<bool, kCells> is_target{};
    for (int col = 0; col < kCols; ++col) {
        if ((mask & (1 << col)) != 0) {
            is_target[static_cast<std::size_t>(make_pos(col, row))] = true;
        }
    }

    std::array<std::array<ld, kCells>, kCells> a{};
    std::array<std::array<ld, kRhsCount>, kCells> rhs{};

    for (int state = 0; state < kCells; ++state) {
        if (is_target[static_cast<std::size_t>(state)]) {
            a[static_cast<std::size_t>(state)][static_cast<std::size_t>(state)] = 1.0L;
            rhs[static_cast<std::size_t>(state)][0] = 0.0L;
            const int col = state % kSize;
            rhs[static_cast<std::size_t>(state)][static_cast<std::size_t>(1 + col)] = 1.0L;
            continue;
        }

        a[static_cast<std::size_t>(state)][static_cast<std::size_t>(state)] = 1.0L;
        const int deg = grid.degree[static_cast<std::size_t>(state)];
        const ld step = 1.0L / static_cast<ld>(deg);
        for (int i = 0; i < deg; ++i) {
            const int nb = grid.neighbors[static_cast<std::size_t>(state)][static_cast<std::size_t>(i)];
            a[static_cast<std::size_t>(state)][static_cast<std::size_t>(nb)] -= step;
        }
        rhs[static_cast<std::size_t>(state)][0] = 1.0L;
    }

    if (!solve_linear_system(a, rhs)) {
        ok = false;
        return data;
    }

    for (int state = 0; state < kCells; ++state) {
        data.expect[static_cast<std::size_t>(state)] = rhs[static_cast<std::size_t>(state)][0];
        for (int col = 0; col < kCols; ++col) {
            data.prob[static_cast<std::size_t>(state)][static_cast<std::size_t>(col)] =
                rhs[static_cast<std::size_t>(state)][static_cast<std::size_t>(1 + col)];
        }
    }

    return data;
}

struct SolutionData {
    std::array<HittingData, kMaskCount> hit_top{};
    std::array<HittingData, kMaskCount> hit_bottom{};
    std::array<std::array<std::array<ld, kCells>, kMaskCount>, kMaskCount> expected_nc{};
    std::array<std::array<std::array<ld, kCells>, kMaskCount>, kMaskCount> expected_c{};
};

bool build_solution_data(SolutionData& out) {
    const Grid grid = build_grid();

    bool ok = true;
    for (int mask = 1; mask < kMaskCount; ++mask) {
        out.hit_top[static_cast<std::size_t>(mask)] = solve_hitting_data(grid, 0, mask, ok);
        if (!ok) {
            return false;
        }
        out.hit_bottom[static_cast<std::size_t>(mask)] = solve_hitting_data(grid, kSize - 1, mask, ok);
        if (!ok) {
            return false;
        }
    }

    std::array<int, kMaskCount> popcount{};
    for (int mask = 0; mask < kMaskCount; ++mask) {
        popcount[static_cast<std::size_t>(mask)] = std::popcount(static_cast<unsigned>(mask));
    }

    for (int pos = 0; pos < kCells; ++pos) {
        out.expected_nc[0][kFullMask][static_cast<std::size_t>(pos)] = 0.0L;
    }

    for (int top_count = 4; top_count >= 0; --top_count) {
        // Carrying states at this top_count use non-carrying states at top_count + 1.
        for (int top_mask = 0; top_mask < kMaskCount; ++top_mask) {
            if (popcount[static_cast<std::size_t>(top_mask)] != top_count) {
                continue;
            }
            const int empty_top_mask = kFullMask ^ top_mask;
            const HittingData& drop = out.hit_top[static_cast<std::size_t>(empty_top_mask)];

            for (int bottom_mask = 0; bottom_mask < kMaskCount; ++bottom_mask) {
                if (popcount[static_cast<std::size_t>(bottom_mask)] != 4 - top_count) {
                    continue;
                }

                for (int pos = 0; pos < kCells; ++pos) {
                    ld value = drop.expect[static_cast<std::size_t>(pos)];

                    for (int col = 0; col < kCols; ++col) {
                        const int bit = 1 << col;
                        if ((empty_top_mask & bit) == 0) {
                            continue;
                        }
                        value +=
                            drop.prob[static_cast<std::size_t>(pos)][static_cast<std::size_t>(col)] *
                            out.expected_nc[static_cast<std::size_t>(bottom_mask)]
                                           [static_cast<std::size_t>(top_mask | bit)]
                                           [static_cast<std::size_t>(top_pos(col))];
                    }

                    out.expected_c[static_cast<std::size_t>(bottom_mask)]
                                  [static_cast<std::size_t>(top_mask)]
                                  [static_cast<std::size_t>(pos)] = value;
                }
            }
        }

        // Non-carrying states at this top_count use carrying states at same top_count.
        for (int top_mask = 0; top_mask < kMaskCount; ++top_mask) {
            if (popcount[static_cast<std::size_t>(top_mask)] != top_count) {
                continue;
            }
            for (int bottom_mask = 0; bottom_mask < kMaskCount; ++bottom_mask) {
                if (popcount[static_cast<std::size_t>(bottom_mask)] != 5 - top_count) {
                    continue;
                }

                const HittingData& pick = out.hit_bottom[static_cast<std::size_t>(bottom_mask)];
                for (int pos = 0; pos < kCells; ++pos) {
                    ld value = pick.expect[static_cast<std::size_t>(pos)];

                    for (int col = 0; col < kCols; ++col) {
                        const int bit = 1 << col;
                        if ((bottom_mask & bit) == 0) {
                            continue;
                        }
                        value +=
                            pick.prob[static_cast<std::size_t>(pos)][static_cast<std::size_t>(col)] *
                            out.expected_c[static_cast<std::size_t>(bottom_mask ^ bit)]
                                          [static_cast<std::size_t>(top_mask)]
                                          [static_cast<std::size_t>(bottom_pos(col))];
                    }

                    out.expected_nc[static_cast<std::size_t>(bottom_mask)]
                                   [static_cast<std::size_t>(top_mask)]
                                   [static_cast<std::size_t>(pos)] = value;
                }
            }
        }
    }

    return true;
}

bool run_checkpoints(const SolutionData& data) {
    constexpr ld kTol = 1e-10L;

    // Check hitting-system invariants.
    for (int mask = 1; mask < kMaskCount; ++mask) {
        for (const auto* table : {&data.hit_top, &data.hit_bottom}) {
            const HittingData& hit = (*table)[static_cast<std::size_t>(mask)];
            for (int pos = 0; pos < kCells; ++pos) {
                const int row = (table == &data.hit_top) ? 0 : (kSize - 1);
                const int col = pos % kSize;
                const bool is_target = (pos / kSize == row) && ((mask & (1 << col)) != 0);

                ld prob_sum = 0.0L;
                for (int c = 0; c < kCols; ++c) {
                    if ((mask & (1 << c)) != 0) {
                        prob_sum += hit.prob[static_cast<std::size_t>(pos)][static_cast<std::size_t>(c)];
                    }
                }

                if (std::fabsl(prob_sum - 1.0L) > kTol) {
                    std::cerr << "Probability sum checkpoint failed for mask=" << mask
                              << ", pos=" << pos << '\n';
                    return false;
                }
                if (is_target && std::fabsl(hit.expect[static_cast<std::size_t>(pos)]) > kTol) {
                    std::cerr << "Target expectation checkpoint failed for mask=" << mask
                              << ", pos=" << pos << '\n';
                    return false;
                }
            }
        }
    }

    // With one seed left and one top slot empty, the process decomposes into
    // deterministic pick/drop targets.
    for (int seed_col = 0; seed_col < kCols; ++seed_col) {
        const int bottom_mask = 1 << seed_col;
        const int bottom_start = bottom_pos(seed_col);

        for (int empty_col = 0; empty_col < kCols; ++empty_col) {
            const int empty_mask = 1 << empty_col;
            const int top_mask = kFullMask ^ empty_mask;

            const ld carry_from_bottom =
                data.hit_top[static_cast<std::size_t>(empty_mask)]
                    .expect[static_cast<std::size_t>(bottom_start)];

            for (int pos = 0; pos < kCells; ++pos) {
                const ld expected_nc =
                    data.hit_bottom[static_cast<std::size_t>(bottom_mask)]
                        .expect[static_cast<std::size_t>(pos)] +
                    carry_from_bottom;

                const ld actual_nc =
                    data.expected_nc[static_cast<std::size_t>(bottom_mask)]
                                    [static_cast<std::size_t>(top_mask)]
                                    [static_cast<std::size_t>(pos)];
                if (std::fabsl(actual_nc - expected_nc) > 5e-10L) {
                    std::cerr << "One-seed non-carrying checkpoint failed at pos=" << pos
                              << ", seed_col=" << seed_col << ", empty_col=" << empty_col << '\n';
                    return false;
                }

                const ld expected_carry =
                    data.hit_top[static_cast<std::size_t>(empty_mask)]
                        .expect[static_cast<std::size_t>(pos)];
                const ld actual_carry =
                    data.expected_c[0][static_cast<std::size_t>(top_mask)]
                                   [static_cast<std::size_t>(pos)];
                if (std::fabsl(actual_carry - expected_carry) > 5e-10L) {
                    std::cerr << "One-seed carrying checkpoint failed at pos=" << pos
                              << ", empty_col=" << empty_col << '\n';
                    return false;
                }
            }
        }
    }

    return true;
}

ld solve_expected_steps() {
    SolutionData data;
    if (!build_solution_data(data)) {
        return -1.0L;
    }

    if (!run_checkpoints(data)) {
        return -2.0L;
    }

    const int initial_pos = make_pos(2, 2);
    const int initial_bottom_mask = kFullMask;
    const int initial_top_mask = 0;
    return data.expected_nc[static_cast<std::size_t>(initial_bottom_mask)]
                           [static_cast<std::size_t>(initial_top_mask)]
                           [static_cast<std::size_t>(initial_pos)];
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    SolutionData data;
    if (!build_solution_data(data)) {
        std::cerr << "Failed to build linear systems." << '\n';
        return 2;
    }
    if (options.run_checkpoints && !run_checkpoints(data)) {
        return 3;
    }

    const int initial_pos = make_pos(2, 2);
    const ld answer = data.expected_nc[static_cast<std::size_t>(kFullMask)][0]
                                      [static_cast<std::size_t>(initial_pos)];

    std::cout << std::fixed << std::setprecision(6) << static_cast<double>(answer) << '\n';
    return 0;
}
