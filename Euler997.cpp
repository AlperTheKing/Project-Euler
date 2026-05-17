#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;

constexpr int TARGET_X = 9;
constexpr int TARGET_Y = 10;
constexpr int TARGET_Z = 11;

struct Orientation {
    std::array<int, 6> face{};
};

u64 pow2(const int exponent) {
    assert(0 <= exponent && exponent < 63);
    return 1ULL << exponent;
}

u64 f(const int x, const int y, const int z) {
    return 3ULL * (pow2(x) + pow2(y) + pow2(z) - 4ULL) * pow2(x + y + z - 1);
}

int parity(const std::array<int, 3>& p) {
    int inversions = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            if (p[i] > p[j]) {
                ++inversions;
            }
        }
    }
    return inversions % 2 == 0 ? 1 : -1;
}

int signed_axis(const int axis, const int sign) {
    return 2 * axis + (sign < 0 ? 1 : 0);
}

std::vector<Orientation> orientations() {
    std::vector<Orientation> result;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            for (int c = 0; c < 3; ++c) {
                if (a == b || a == c || b == c) {
                    continue;
                }
                const std::array<int, 3> perm{a, b, c};
                for (int mask = 0; mask < 8; ++mask) {
                    std::array<int, 3> sign{};
                    for (int i = 0; i < 3; ++i) {
                        sign[i] = (mask & (1 << i)) == 0 ? 1 : -1;
                    }
                    if (parity(perm) * sign[0] * sign[1] * sign[2] != 1) {
                        continue;
                    }

                    Orientation orientation;
                    for (int axis = 0; axis < 3; ++axis) {
                        orientation.face[2 * axis] = signed_axis(perm[axis], sign[axis]);
                        orientation.face[2 * axis + 1] = signed_axis(perm[axis], -sign[axis]);
                    }
                    result.push_back(orientation);
                }
            }
        }
    }
    return result;
}

u64 brute_count(const int x, const int y, const int z) {
    const std::vector<Orientation> all = orientations();
    std::vector<std::array<int, 3>> cells;
    for (int k = 0; k < z; ++k) {
        for (int j = 0; j < y; ++j) {
            for (int i = 0; i < x; ++i) {
                cells.push_back({i, j, k});
            }
        }
    }

    auto index = [=](const int i, const int j, const int k) {
        return (k * y + j) * x + i;
    };

    std::vector<int> grid(static_cast<std::size_t>(x * y * z), -1);
    u64 total = 0;

    auto search = [&](auto&& self, const int position) -> void {
        if (position == static_cast<int>(cells.size())) {
            ++total;
            return;
        }

        const auto [i, j, k] = cells[static_cast<std::size_t>(position)];
        for (int o = 0; o < static_cast<int>(all.size()); ++o) {
            const Orientation& current = all[static_cast<std::size_t>(o)];
            if (i > 0) {
                const Orientation& previous = all[static_cast<std::size_t>(grid[static_cast<std::size_t>(index(i - 1, j, k))])];
                if (previous.face[0] != current.face[1]) {
                    continue;
                }
            }
            if (j > 0) {
                const Orientation& previous = all[static_cast<std::size_t>(grid[static_cast<std::size_t>(index(i, j - 1, k))])];
                if (previous.face[2] != current.face[3]) {
                    continue;
                }
            }
            if (k > 0) {
                const Orientation& previous = all[static_cast<std::size_t>(grid[static_cast<std::size_t>(index(i, j, k - 1))])];
                if (previous.face[4] != current.face[5]) {
                    continue;
                }
            }

            grid[static_cast<std::size_t>(index(i, j, k))] = o;
            self(self, position + 1);
            grid[static_cast<std::size_t>(index(i, j, k))] = -1;
        }
    };

    search(search, 0);
    return total;
}

void run_checkpoints() {
    assert(orientations().size() == 24);
    assert(f(1, 1, 1) == 24);
    assert(f(2, 3, 4) == 18'432);
    assert(f(1, 1, 2) == brute_count(1, 1, 2));
    assert(f(1, 2, 2) == brute_count(1, 2, 2));
    assert(f(2, 2, 1) == brute_count(2, 2, 1));
    assert(f(2, 2, 2) == brute_count(2, 2, 2));
}

}  // namespace

int main() {
    run_checkpoints();
    std::cout << f(TARGET_X, TARGET_Y, TARGET_Z) << '\n';
    return 0;
}
