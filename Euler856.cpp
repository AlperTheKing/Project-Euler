#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <vector>

using u64 = std::uint64_t;

class Solver {
public:
    Solver(int ranks, int copies)
        : ranks_(ranks), copies_(copies), base_(ranks + 1) {}

    long double expected_draws() {
        if (ranks_ <= 0 || copies_ <= 0) return 0.0L;
        std::vector<int> a(copies_ + 1, 0);
        if (ranks_ >= 1) a[copies_] = ranks_ - 1;
        int m = copies_ - 1;
        memo_.clear();
        return 1.0L + dfs(m, a);
    }

private:
    int ranks_;
    int copies_;
    int base_;
    std::unordered_map<u64, long double> memo_;

    u64 encode(int m, const std::vector<int>& a) const {
        u64 key = static_cast<u64>(m);
        for (int v : a) key = key * static_cast<u64>(base_) + static_cast<u64>(v);
        return key;
    }

    long double dfs(int m, const std::vector<int>& a) {
        int remaining = m;
        for (int i = 1; i <= copies_; ++i) remaining += i * a[i];
        if (remaining == 0) return 0.0L;

        u64 key = encode(m, a);
        auto it = memo_.find(key);
        if (it != memo_.end()) return it->second;

        long double ans = 1.0L;
        for (int i = 1; i <= copies_; ++i) {
            if (a[i] == 0) continue;
            long double p = static_cast<long double>(i * a[i]) / static_cast<long double>(remaining);
            std::vector<int> b = a;
            --b[i];
            ++b[m];
            ans += p * dfs(i - 1, b);
        }

        memo_[key] = ans;
        return ans;
    }
};

int main() {
    assert(std::fabsl(3.0L / 51.0L - 1.0L / 17.0L) < 1e-18L);
    assert(std::fabsl(Solver(1, 4).expected_draws() - 2.0L) < 1e-18L);
    assert(std::fabsl(Solver(2, 2).expected_draws() - 3.0L) < 1e-18L);

    long double ans = Solver(13, 4).expected_draws();
    std::cout << std::fixed << std::setprecision(8) << ans << '\n';
    return 0;
}
