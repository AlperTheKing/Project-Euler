#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Key {
    i64 x;
    i64 m;

    bool operator==(const Key& other) const noexcept { return x == other.x && m == other.m; }
};

struct KeyHash {
    std::size_t operator()(const Key& k) const noexcept {
        std::size_t h1 = std::hash<i64>{}(k.x);
        std::size_t h2 = std::hash<i64>{}(k.m);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6U) + (h1 >> 2U));
    }
};

class Solver702 {
  public:
    i64 solve(i64 n) {
        cache_.clear();
        const int D = bit_length(n);

        i64 ans = static_cast<i64>((static_cast<u128>(n) * (3 * static_cast<u128>(n) + 1) / 2) *
                                   static_cast<u128>(D + 1));

        for (int d = 2; d <= D; ++d) {
            ans -= g(n, 1LL << d);
        }
        ans += 2 * g(n, (1LL << D) - n);
        return ans;
    }

  private:
    std::unordered_map<Key, i64, KeyHash> cache_;

    static int bit_length(i64 x) {
        int len = 0;
        while (x > 0) {
            ++len;
            x >>= 1;
        }
        return len;
    }

    i64 f(i64 x, i64 m) {
        x %= m;
        if (x <= 1) return 0;

        const Key key{x, m};
        auto it = cache_.find(key);
        if (it != cache_.end()) return it->second;

        const i64 t = m / x;
        const i64 y = m % x;

        const i64 term = static_cast<i64>((static_cast<u128>(t) * (t + 1) * x * (x - 1)) / 4);
        const i64 res = term + (t + 1) * f(x, y) - t * f(x, x - y);

        cache_.emplace(key, res);
        return res;
    }

    i64 g(i64 x, i64 m) { return (m - 1) * (m - 2) - f(x, m); }
};

void run_validations() {
    Solver702 solver;

    struct Test {
        i64 N;
        i64 expected;
    };

    const Test tests[] = {
        {3, 42},
        {5, 126},
        {123, 167178},
        {12345, 3185041956LL},
    };

    for (const auto& t : tests) {
        i64 got = solver.solve(t.N);
        if (got != t.expected) {
            std::cerr << "Validation failed for N=" << t.N << ": got " << got
                      << ", expected " << t.expected << "\n";
            std::exit(1);
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    bool validate = false;
    i64 N = 123456789LL;

    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        if (s == "--validate") {
            validate = true;
        } else {
            N = std::stoll(s);
        }
    }

    if (validate) {
        run_validations();
    }

    Solver702 solver;
    std::cout << solver.solve(N) << '\n';
    return 0;
}
