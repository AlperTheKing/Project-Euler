#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

u64 solve_f(u64 n) {
    std::vector<u128> fib = {0, 1};
    fib.reserve(256);
    for (int i = 0; i < 240; ++i) {
        fib.push_back(fib[fib.size() - 1] + fib[fib.size() - 2]);
    }

    int depth = 0;
    while (fib[static_cast<std::size_t>(depth + 3)] < n) {
        ++depth;
    }

    const u64 kInf = std::numeric_limits<u64>::max();

    struct State {
        u64 x;
        u64 y;
        i64 p;
        i64 q;
        int rem;
    };

    while (true) {
        u64 best_m = kInf;
        std::vector<State> stack;
        stack.reserve(1 << 20);
        stack.push_back(State{2ULL, 1ULL, 0, 1, depth});

        while (!stack.empty()) {
            const State s = stack.back();
            stack.pop_back();

            if (s.y >= best_m) {
                continue;
            }
            if (s.x == n) {
                best_m = s.y;
                continue;
            }
            if (s.rem == 0 || s.x > n) {
                continue;
            }

            if (static_cast<u128>(s.x) + static_cast<u128>(s.rem) * s.y > n) {
                continue;
            }
            if (fib[static_cast<std::size_t>(s.rem + 1)] * s.x + fib[static_cast<std::size_t>(s.rem)] * s.y < n) {
                continue;
            }

            const i64 y_i = static_cast<i64>(s.y);
            const i64 p_mod = (s.p % y_i + y_i) % y_i;
            const i64 a = static_cast<i64>((static_cast<u128>(n) * static_cast<u64>(p_mod)) % s.y);
            if (static_cast<u128>(a) * s.x > n) {
                continue;
            }

            const u64 xp = s.x + s.y;
            const int rem_next = s.rem - 1;
            stack.push_back(State{xp, s.x, s.q, s.p - s.q, rem_next});
            stack.push_back(State{xp, s.y, s.p, s.q - s.p, rem_next});
        }

        if (best_m != kInf) {
            return best_m;
        }
        ++depth;
    }
}

void run_validations() {
    assert(solve_f(7ULL) == 2ULL);
    assert(solve_f(89ULL) == 34ULL);
    assert(solve_f(8'191ULL) == 1'856ULL);
}

}  // namespace

int main() {
    run_validations();
    std::cout << solve_f(1'000'000'000'039ULL) << '\n';
    return 0;
}
