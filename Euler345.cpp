#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;

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

i64 matrix_sum(const std::vector<std::vector<int>>& a) {
    const int n = static_cast<int>(a.size());
    const int full = 1 << n;
    std::vector<i64> dp(static_cast<std::size_t>(full), -1);
    dp[0] = 0;

    for (int mask = 0; mask < full; ++mask) {
        const i64 base = dp[static_cast<std::size_t>(mask)];
        if (base < 0) {
            continue;
        }
        const int row = __builtin_popcount(static_cast<unsigned>(mask));
        if (row >= n) {
            continue;
        }
        for (int col = 0; col < n; ++col) {
            if ((mask & (1 << col)) != 0) {
                continue;
            }
            const int next = mask | (1 << col);
            dp[static_cast<std::size_t>(next)] =
                std::max(dp[static_cast<std::size_t>(next)],
                         base + static_cast<i64>(a[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]));
        }
    }
    return dp[static_cast<std::size_t>(full - 1)];
}

i64 solve() {
    const std::vector<std::vector<int>> m = {
        {7,53,183,439,863,497,383,563,79,973,287,63,343,169,583},
        {627,343,773,959,943,767,473,103,699,303,957,703,583,639,913},
        {447,283,463,29,23,487,463,993,119,883,327,493,423,159,743},
        {217,623,3,399,853,407,103,983,89,463,290,516,212,462,350},
        {960,376,682,962,300,780,486,502,912,800,250,346,172,812,350},
        {870,456,192,162,593,473,915,45,989,873,823,965,425,329,803},
        {973,965,905,919,133,673,665,235,509,613,673,815,165,992,326},
        {322,148,972,962,286,255,941,541,265,323,925,281,601,95,973},
        {445,721,11,525,473,65,511,164,138,672,18,428,154,448,848},
        {414,456,310,312,798,104,566,520,302,248,694,976,430,392,198},
        {184,829,373,181,631,101,969,613,840,740,778,458,284,760,390},
        {821,461,843,513,17,901,711,993,293,157,274,94,192,156,574},
        {34,124,4,878,450,476,712,914,838,669,875,299,823,329,699},
        {815,559,813,459,522,788,168,586,966,232,308,833,251,631,107},
        {813,883,451,509,615,77,281,613,459,205,380,274,302,35,805}
    };
    return matrix_sum(m);
}

bool run_checkpoints() {
    const std::vector<std::vector<int>> sample = {
        {7,53,183,439,863},
        {497,383,563,79,973},
        {287,63,343,169,583},
        {627,343,773,959,943},
        {767,473,103,699,303}
    };
    if (matrix_sum(sample) != 3315) {
        std::cerr << "Checkpoint failed for 5x5 sample matrix" << '\n';
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
    std::cout << solve() << '\n';
    return 0;
}
