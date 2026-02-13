#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;

class Solver {
public:
    Solver() {
        for (auto& row : memo_) {
            row.fill(-1);
        }
    }

    u64 W_pow10(int exp) {
        std::array<u64, 19> pow9{};
        pow9[0] = 1ULL;
        for (int i = 1; i <= 18; ++i) {
            pow9[static_cast<std::size_t>(i)] = pow9[static_cast<std::size_t>(i - 1)] * 9ULL;
        }

        u64 ans = 0ULL;
        for (int len = 1; len <= exp; ++len) {
            const u64 base = 1ULL << (len - 1);
            const u64 limit = 1ULL << (len - 1);
            for (u64 tail = 0ULL; tail < limit; ++tail) {
                const u64 mask = base | tail;
                if (is_winning(mask, len)) {
                    const int ones = std::popcount(mask);
                    ans += pow9[static_cast<std::size_t>(ones)];
                }
            }
        }
        return ans;
    }

private:
    std::array<std::array<int8_t, 1U << 18>, 19> memo_{};

    bool is_winning(u64 mask, int len) {
        if (len == 0) {
            return false;
        }
        int8_t& cell = memo_[static_cast<std::size_t>(len)][static_cast<std::size_t>(mask)];
        if (cell != -1) {
            return cell != 0;
        }

        for (int i = 0; i < len; ++i) {
            const int j = len - 1 - i;
            const u64 high = mask >> (j + 1);
            const u64 low = mask & ((1ULL << j) - 1ULL);
            u64 next_mask = (high << j) | low;
            int next_len = len - 1;

            while (next_len > 0 && ((next_mask >> (next_len - 1)) & 1ULL) == 0ULL) {
                --next_len;
            }

            if (!is_winning(next_mask, next_len)) {
                cell = 1;
                return true;
            }
        }

        cell = 0;
        return false;
    }
};

void run_validations() {
    Solver solver;
    assert(solver.W_pow10(2) == 18ULL);
    assert(solver.W_pow10(4) == 1'656ULL);
}

}  // namespace

int main() {
    run_validations();
    Solver solver;
    std::cout << solver.W_pow10(18) << '\n';
    return 0;
}
