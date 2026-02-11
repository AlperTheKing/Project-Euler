#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;

std::vector<u64> fib_upto(u64 n) {
    std::vector<u64> f{0ULL, 1ULL, 2ULL};
    while (f.back() <= n) {
        f.push_back(f[f.size() - 1] + f[f.size() - 2]);
    }
    return f;
}

u64 H(u64 n, const std::vector<u64>& f) {
    u64 x = n;
    u64 best = n;

    for (int i = static_cast<int>(f.size()) - 1; i >= 1; --i) {
        if (f[static_cast<std::size_t>(i)] <= x) {
            x -= f[static_cast<std::size_t>(i)];
            if (f[static_cast<std::size_t>(i)] < best) {
                best = f[static_cast<std::size_t>(i)];
            }
        }
        if (x == 0ULL) {
            break;
        }
    }

    return best;
}

u64 G(u64 n) {
    const std::vector<u64> f = fib_upto(n);
    std::vector<u64> pref(f.size(), 0ULL);

    pref[1] = 0ULL;
    pref[2] = 1ULL;

    for (std::size_t m = 2; m + 1 < f.size(); ++m) {
        pref[m + 1] = pref[m] + f[m] + pref[m - 1];
    }

    u64 x = n;
    u64 out = 0ULL;

    while (x > 0ULL) {
        const std::size_t m = static_cast<std::size_t>(
            std::upper_bound(f.begin(), f.end(), x) - f.begin() - 1);
        out += pref[m] + f[m];
        x -= f[m];
    }

    return out;
}

}  // namespace

int main() {
    const std::vector<u64> f = fib_upto(1'000ULL);

    assert(H(1ULL, f) == 1ULL);
    assert(H(4ULL, f) == 1ULL);
    assert(H(17ULL, f) == 1ULL);
    assert(H(8ULL, f) == 8ULL);
    assert(H(18ULL, f) == 5ULL);

    assert(G(13ULL) == 43ULL);

    std::cout << G(23'416'728'348'467'685ULL) << "\n";
    return 0;
}
