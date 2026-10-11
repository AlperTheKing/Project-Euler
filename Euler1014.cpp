#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

using Mask = std::uint32_t;
constexpr unsigned MEMBERS = 24;
constexpr Mask ALL = (Mask{1} << MEMBERS) - 1;
constexpr std::array<const char*, 11> GIVEN = {
    "ABCDKMPR", "ABJKLMNW", "ACIOQTUW", "ADFMPQTU", "AEFNPTVX", "AEFQRSTW",
    "AEGIMQSU", "BDEKORUX", "BEFGHIKU", "DEFHKOSW", "FHIJKRST"
};

void require(const bool condition, const std::string& description) {
    if (!condition) throw std::runtime_error("Check failed: " + description);
}

unsigned binomial(const unsigned n, const unsigned k) {
    if (k > n) return 0;
    unsigned result = 1;
    for (unsigned i = 1; i <= k; ++i) result = result * (n - i + 1) / i;
    return result;
}

unsigned weight(Mask mask) {
    unsigned count = 0;
    for (; mask; mask &= mask - 1) ++count;
    return count;
}

Mask encode(const std::string& names) {
    Mask mask = 0;
    for (const char name : names) {
        require(name >= 'A' && name <= 'X', "member name");
        const Mask bit = Mask{1} << (name - 'A');
        require((mask & bit) == 0, "distinct members");
        mask |= bit;
    }
    return mask;
}

std::string decode(const Mask mask) {
    std::string names;
    for (unsigned i = 0; i < MEMBERS; ++i) {
        if (mask & (Mask{1} << i)) names += static_cast<char>('A' + i);
    }
    return names;
}

std::array<unsigned, 9> intersection_counts() {
    std::array<unsigned, 9> counts{};
    counts[8] = 1;
    for (unsigned s = 5; s-- > 0;) {
        const unsigned numerator = binomial(24 - s, 5 - s);
        const unsigned denominator = binomial(8 - s, 5 - s);
        require(numerator % denominator == 0, "integral subset multiplicity");
        counts[s] = binomial(8, s) * (numerator / denominator);
        for (unsigned i = s + 1; i <= 8; ++i) {
            const unsigned contribution = binomial(i, s) * counts[i];
            require(counts[s] >= contribution, "nonnegative intersection count");
            counts[s] -= contribution;
        }
    }
    return counts;
}

Mask solve(const std::vector<Mask>& blocks, const Mask target) {
    require(weight(target) == 5, "five target members");
    const auto counts = intersection_counts();
    std::vector<Mask> remaining;
    for (unsigned i = 0; i < MEMBERS; ++i) {
        const Mask bit = Mask{1} << i;
        if (!(target & bit)) remaining.push_back(bit);
    }
    Mask answer = 0;
    for (std::size_t a = 0; a + 2 < remaining.size(); ++a) {
        for (std::size_t b = a + 1; b + 1 < remaining.size(); ++b) {
            for (std::size_t c = b + 1; c < remaining.size(); ++c) {
                const Mask candidate = target | remaining[a] | remaining[b] | remaining[c];
                const bool possible = std::all_of(blocks.begin(), blocks.end(), [&](const Mask block) {
                    return counts[weight(candidate & block)] != 0;
                });
                if (!possible) continue;
                require(answer == 0, "unique compatible night");
                answer = candidate;
            }
        }
    }
    require(answer != 0, "compatible night exists");
    return answer;
}

void run_tests(const std::vector<Mask>& given, const Mask target, const Mask answer) {
    const auto counts = intersection_counts();
    require(counts == std::array<unsigned, 9>{30, 0, 448, 0, 280, 0, 0, 0, 1},
        "intersection distribution");

    std::vector<Mask> generators = given;
    generators.push_back(ALL);
    std::vector<Mask> span{0};
    for (const Mask generator : generators) {
        const std::size_t size = span.size();
        for (std::size_t i = 0; i < size; ++i) span.push_back(span[i] ^ generator);
    }
    std::sort(span.begin(), span.end());
    require(std::adjacent_find(span.begin(), span.end()) == span.end(), "independent generators");

    std::vector<Mask> blocks;
    for (const Mask word : span) {
        if (weight(word) == 8) blocks.push_back(word);
    }
    const auto group_count = static_cast<std::size_t>(binomial(24, 5));
    require(blocks.size() * binomial(8, 5) == group_count, "total number of nights");
    for (const Mask block : given) {
        require(std::binary_search(blocks.begin(), blocks.end(), block), "given night retained");
    }

    std::unordered_set<Mask> covered;
    covered.reserve(group_count);
    Mask reconstructed = 0;
    for (const Mask block : blocks) {
        for (Mask subset = block; subset; subset = (subset - 1) & block) {
            if (weight(subset) == 5) require(covered.insert(subset).second, "five-member group occurs once");
        }
        if ((block & target) == target) {
            require(reconstructed == 0, "unique reconstructed night");
            reconstructed = block;
        }
    }
    require(covered.size() == group_count, "all five-member groups covered");
    require(reconstructed == answer, "agreement with full design reconstruction");
    std::cout << "All checks passed; " << blocks.size() << " nights cover " << covered.size()
              << " five-member groups exactly once.\n";
}

}

int main(int argc, char* argv[]) {
    try {
        const bool self_test = argc == 2 && std::string(argv[1]) == "--self-test";
        if (argc != 1 && !self_test) throw std::invalid_argument("Usage: Euler1014 [--self-test]");
        std::vector<Mask> blocks;
        for (const char* names : GIVEN) {
            blocks.push_back(encode(names));
            require(weight(blocks.back()) == 8, "eight members per given night");
        }
        const Mask target = encode("CGJLV");
        const Mask answer = solve(blocks, target);
        if (self_test) run_tests(blocks, target, answer);
        else std::cout << decode(answer) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
