#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

// Project Euler 599: Distinct Colourings of a Rubik's Cube
//
// A 2x2x2 cube has 8 corner cubies, each with 3 stickers, for 24 stickers total.
// We colour each sticker with one of n colours (unlimited supply).
// Two colourings are equivalent if one can be transformed into the other by legal cube moves.
//
// By Burnside's lemma, the number of essentially distinct colourings is:
//   (1/|G|) * sum_{g in G} n^{c(g)}
// where c(g) is the number of cycles in the permutation of the 24 stickers induced by g.
//
// The move group on corners consists of all corner permutations (8!) and orientations (3^7,
// last orientation determined by the twist-sum constraint), so |G| = 8! * 3^7.
//
// We enumerate all group elements in corner representation and compute the cycle-count
// histogram for the induced 24-sticker permutation.

using boost::multiprecision::cpp_int;
using u64 = std::uint64_t;

static int cycles_in_perm(const std::array<std::uint8_t, 24>& p) {
    bool seen[24] = {false};
    int cycles = 0;
    for (int i = 0; i < 24; ++i) {
        if (seen[i]) continue;
        ++cycles;
        int x = i;
        while (!seen[x]) {
            seen[x] = true;
            x = p[(size_t)x];
        }
    }
    return cycles;
}

static std::array<u64, 25> cycle_histogram() {
    // Precompute all orientation vectors (3^7) with sum == 0 (mod 3).
    std::vector<std::array<std::uint8_t, 8>> oris;
    oris.reserve(2187);
    for (int mask = 0; mask < 2187; ++mask) {
        std::array<std::uint8_t, 8> ori{};
        int t = mask;
        int sum = 0;
        for (int i = 0; i < 7; ++i) {
            const int o = t % 3;
            t /= 3;
            ori[(size_t)i] = (std::uint8_t)o;
            sum += o;
        }
        ori[7] = (std::uint8_t)((3 - (sum % 3)) % 3);
        oris.push_back(ori);
    }

    std::array<u64, 25> freq{};
    freq.fill(0);

    std::array<std::uint8_t, 8> perm{};
    for (int i = 0; i < 8; ++i) perm[(size_t)i] = (std::uint8_t)i;

    do {
        for (const auto& ori : oris) {
            std::array<std::uint8_t, 24> p{};
            for (int pos = 0; pos < 8; ++pos) {
                const int c = perm[(size_t)pos];
                const int o = ori[(size_t)pos];
                for (int j = 0; j < 3; ++j) {
                    // Slot j at position pos receives cubie facelet k.
                    const int k = (j - o + 3) % 3;
                    p[(size_t)(3 * c + k)] = (std::uint8_t)(3 * pos + j);
                }
            }
            const int cyc = cycles_in_perm(p);
            ++freq[(size_t)cyc];
        }
    } while (std::next_permutation(perm.begin(), perm.end()));

    return freq;
}

static cpp_int orbits_from_hist(const std::array<u64, 25>& freq, int colours) {
    const u64 G = 88179840ULL; // 8! * 3^7

    cpp_int sum = 0;
    cpp_int pw = 1;
    for (int k = 0; k <= 24; ++k) {
        if (freq[(size_t)k]) sum += cpp_int(freq[(size_t)k]) * pw;
        pw *= colours;
    }
    return sum / G;
}

int main() {
    const auto freq = cycle_histogram();

    // Validation from the statement.
    if (orbits_from_hist(freq, 2) != 183) {
        std::cerr << "Validation failed for n=2\n";
        return 1;
    }

    std::cout << orbits_from_hist(freq, 10) << "\n";
    return 0;
}
