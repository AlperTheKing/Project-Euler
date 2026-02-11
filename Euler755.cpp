#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;

u64 S(const u64 n) {
    std::vector<u64> fib = {1ULL, 2ULL};
    while (fib[fib.size() - 1] + fib[fib.size() - 2] <= n) {
        fib.push_back(fib[fib.size() - 1] + fib[fib.size() - 2]);
    }

    std::vector<u64> prefix(fib.size(), 0ULL);
    u64 running = 0ULL;
    for (std::size_t i = 0; i < fib.size(); ++i) {
        running += fib[i];
        prefix[i] = running;
    }

    std::vector<std::unordered_map<u64, u64>> memo(fib.size());

    const auto dfs = [&](auto&& self, const int i, const u64 limit) -> u64 {
        if (i < 0) {
            return 1ULL;
        }
        if (limit >= prefix[static_cast<std::size_t>(i)]) {
            return 1ULL << static_cast<unsigned>(i + 1);
        }

        auto& mp = memo[static_cast<std::size_t>(i)];
        const auto it = mp.find(limit);
        if (it != mp.end()) {
            return it->second;
        }

        u64 total = self(self, i - 1, limit);
        if (limit >= fib[static_cast<std::size_t>(i)]) {
            total += self(self, i - 1, limit - fib[static_cast<std::size_t>(i)]);
        }

        mp.emplace(limit, total);
        return total;
    };

    return dfs(dfs, static_cast<int>(fib.size()) - 1, n);
}

}  // namespace

int main() {
    assert(S(100ULL) == 415ULL);
    assert(S(10'000ULL) == 312'807ULL);

    std::cout << S(10'000'000'000'000ULL) << '\n';
    return 0;
}
