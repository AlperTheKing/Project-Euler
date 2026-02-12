#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

using i64 = std::int64_t;

static constexpr i64 kMod = 1'005'075'251LL;

struct Challenge {
    i64 target;
    std::vector<i64> nums;
};

static std::vector<Challenge> load_challenges(const std::string& path) {
    std::ifstream fin(path);
    if (!fin) {
        throw std::runtime_error("cannot open challenge file");
    }

    std::vector<Challenge> out;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::size_t col = line.find(':');
        if (col == std::string::npos) continue;

        Challenge c;
        c.target = std::stoll(line.substr(0, col));

        std::string rest = line.substr(col + 1);
        std::stringstream ss(rest);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            c.nums.push_back(std::stoll(tok));
        }
        out.push_back(std::move(c));
    }

    return out;
}

static i64 min_score_for_target(const Challenge& c) {
    const int n = static_cast<int>(c.nums.size());
    const int all_masks = 1 << n;

    std::vector<i64> subset_sum(all_masks, 0);
    for (int mask = 1; mask < all_masks; ++mask) {
        int bit = __builtin_ctz(mask);
        int prev = mask & (mask - 1);
        subset_sum[mask] = subset_sum[prev] + c.nums[bit];
    }

    std::vector<std::vector<i64>> vals(all_masks);

    for (int mask = 1; mask < all_masks; ++mask) {
        if ((mask & (mask - 1)) == 0) {
            int bit = __builtin_ctz(mask);
            vals[mask].push_back(c.nums[bit]);
            continue;
        }

        std::unordered_set<i64> seen;
        for (int sub = (mask - 1) & mask; sub > 0; sub = (sub - 1) & mask) {
            int other = mask ^ sub;
            if (sub > other) continue;

            const auto& A = vals[sub];
            const auto& B = vals[other];

            for (i64 a : A) {
                for (i64 b : B) {
                    seen.insert(a + b);
                    seen.insert(a * b);

                    if (a > b) seen.insert(a - b);
                    if (b > a) seen.insert(b - a);

                    if (b != 0 && a % b == 0) seen.insert(a / b);
                    if (a != 0 && b % a == 0) seen.insert(b / a);
                }
            }
        }

        vals[mask].assign(seen.begin(), seen.end());
    }

    i64 best = std::numeric_limits<i64>::max();
    for (int mask = 1; mask < all_masks; ++mask) {
        const auto& v = vals[mask];
        if (std::find(v.begin(), v.end(), c.target) != v.end()) {
            best = std::min(best, subset_sum[mask]);
        }
    }

    if (best == std::numeric_limits<i64>::max()) {
        return 0;
    }
    return best;
}

int main() {
    const auto challenges = load_challenges("resources/documents/0828_number_challenges.txt");
    assert(challenges.size() == 200);

    const i64 s1 = min_score_for_target(challenges[0]);
    assert(s1 == 40);

    i64 ans = 0;
    i64 p3 = 1;

    for (std::size_t i = 0; i < challenges.size(); ++i) {
        p3 = (p3 * 3) % kMod;
        i64 s = min_score_for_target(challenges[i]);
        ans = (ans + (p3 * (s % kMod)) % kMod) % kMod;
    }

    std::cout << ans << '\n';
    return 0;
}
