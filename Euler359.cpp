#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kTargetProduct = 71328803586048ULL;
constexpr u64 kMod = 100000000ULL;

u64 isqrt_u64(const u64 n) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1ULL) <= n / (r + 1ULL)) {
        ++r;
    }
    while (r > n / r) {
        --r;
    }
    return r;
}

u128 person_in_room(const u64 floor, const u64 room) {
    u128 first_person = 0;
    u64 first_square_root = 0;

    if (floor == 1ULL) {
        first_person = 1;
        first_square_root = 2;
    } else {
        first_person = (static_cast<u128>(floor) * floor) / 2U;
        first_square_root = (floor & 1ULL) ? floor : (floor + 1ULL);
    }

    if ((room & 1ULL) != 0ULL) {
        const u128 k = (room + 1ULL) / 2ULL;
        return first_person + (k - 1U) * (2U * static_cast<u128>(first_square_root) + 2U * k - 3U);
    }

    const u128 k = room / 2ULL;
    return static_cast<u128>(first_square_root) * first_square_root - first_person +
           (k - 1U) * (2U * static_cast<u128>(first_square_root) + 2U * k - 1U);
}

u64 sum_for_product_last8(const u64 product) {
    const u64 root = isqrt_u64(product);
    u128 total = 0;

    for (u64 d = 1; d <= root; ++d) {
        if (product % d != 0ULL) {
            continue;
        }
        const u64 q = product / d;
        total += person_in_room(d, q);
        if (d != q) {
            total += person_in_room(q, d);
        }
    }

    return static_cast<u64>(total % kMod);
}

bool validate_formula_against_simulation() {
    const u64 people_limit = 30000ULL;
    std::vector<u64> tops;
    std::vector<std::vector<u64>> floors;
    tops.reserve(512);
    floors.reserve(512);

    for (u64 person = 1; person <= people_limit; ++person) {
        bool placed = false;
        for (std::size_t i = 0; i < tops.size(); ++i) {
            const u64 s = tops[i] + person;
            const u64 r = isqrt_u64(s);
            if (r * r == s) {
                tops[i] = person;
                floors[i].push_back(person);
                placed = true;
                break;
            }
        }
        if (!placed) {
            tops.push_back(person);
            floors.push_back({person});
        }
    }

    for (u64 f = 1; f <= 60; ++f) {
        if (f > floors.size()) {
            break;
        }
        for (u64 r = 1; r <= 60; ++r) {
            if (floors[static_cast<std::size_t>(f - 1ULL)].size() < r) {
                break;
            }
            const u64 brute = floors[static_cast<std::size_t>(f - 1ULL)][static_cast<std::size_t>(r - 1ULL)];
            const u64 fast = static_cast<u64>(person_in_room(f, r));
            if (brute != fast) {
                std::cerr << "Checkpoint failed: mismatch at floor " << f << ", room " << r << '\n';
                return false;
            }
        }
    }
    return true;
}

bool run_checkpoints() {
    if (person_in_room(1ULL, 1ULL) != 1U || person_in_room(1ULL, 2ULL) != 3U || person_in_room(2ULL, 1ULL) != 2U) {
        std::cerr << "Checkpoint failed: base examples\n";
        return false;
    }

    if (person_in_room(10ULL, 20ULL) != 440U) {
        std::cerr << "Checkpoint failed: P(10,20)\n";
        return false;
    }
    if (person_in_room(25ULL, 75ULL) != 4863U) {
        std::cerr << "Checkpoint failed: P(25,75)\n";
        return false;
    }
    if (person_in_room(99ULL, 100ULL) != 19454U) {
        std::cerr << "Checkpoint failed: P(99,100)\n";
        return false;
    }

    if (sum_for_product_last8(20ULL) != 608ULL) {
        std::cerr << "Checkpoint failed: sum for product 20\n";
        return false;
    }

    if (!validate_formula_against_simulation()) {
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!skip_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const u64 answer = sum_for_product_last8(kTargetProduct);
    std::cout << answer << '\n';
    return 0;
}
