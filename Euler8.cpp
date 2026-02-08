#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;

constexpr const char* DIGITS =
    "73167176531330624919225119674426574742355349194934"
    "96983520312774506326239578318016984801869478851843"
    "85861560789112949495459501737958331952853208805511"
    "12540698747158523863050715693290963295227443043557"
    "66896648950445244523161731856403098711121722383113"
    "62229893423380308135336276614282806444486645238749"
    "30358907296290491560440772390713810515859307960866"
    "70172427121883998797908792274921901699720888093776"
    "65727333001053367881220235421809751254540594752243"
    "52584907711670556013604839586446706324415722155397"
    "53697817977846174064955149290862569321978468622482"
    "83972241375657056057490261407972968652414535100474"
    "82166370484403199890008895243450658541227588666881"
    "16427171479924442928230863465674813919123162824586"
    "17866458359124566529476545682848912883142607690042"
    "24219022671055626321111109370544217506941658960408"
    "07198403850962455444362981230987879927244284909188"
    "84580156166097919133875499200524063689912560717606"
    "05886116467109405077541002256983155200055935729725"
    "71636269561882670428252483600823257530420752963450";

struct Options {
    int window = 13;
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
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }

    value = parsed;
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--window=", options.window)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.window > 0;
}

u64 max_adjacent_product(const std::string& digits, const int window) {
    if (window <= 0 || static_cast<std::size_t>(window) > digits.size()) {
        return 0ULL;
    }

    u64 best = 0ULL;
    for (std::size_t i = 0; i + static_cast<std::size_t>(window) <= digits.size(); ++i) {
        u64 prod = 1ULL;
        for (int j = 0; j < window; ++j) {
            const char ch = digits[i + static_cast<std::size_t>(j)];
            prod *= static_cast<u64>(ch - '0');
        }
        if (prod > best) {
            best = prod;
        }
    }

    return best;
}

u64 solve(const int window) {
    return max_adjacent_product(DIGITS, window);
}

bool run_checkpoints() {
    if (max_adjacent_product("123456789", 2) != 72ULL) {
        std::cerr << "Checkpoint failed for custom digits" << '\n';
        return false;
    }
    if (solve(4) != 5832ULL) {
        std::cerr << "Checkpoint failed for window=4" << '\n';
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

    std::cout << solve(options.window) << '\n';
    return 0;
}
