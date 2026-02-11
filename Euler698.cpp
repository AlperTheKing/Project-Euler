#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kTargetN = 111'111'111'111'222'333ULL;
constexpr u64 kMod = 123'123'123ULL;
constexpr int kIs123Limit = 5000;

u64 sat_add(const u64 a, const u64 b, const u64 cap) {
    const u128 v = static_cast<u128>(a) + static_cast<u128>(b);
    return (v > static_cast<u128>(cap)) ? (cap + 1U) : static_cast<u64>(v);
}

u64 sat_mul(const u64 a, const u64 b, const u64 cap) {
    const u128 v = static_cast<u128>(a) * static_cast<u128>(b);
    return (v > static_cast<u128>(cap)) ? (cap + 1U) : static_cast<u64>(v);
}

std::vector<std::uint8_t> build_is123(const int limit) {
    std::vector<std::uint8_t> is123(static_cast<std::size_t>(limit) + 1U, 0U);
    is123[1] = 1U;

    for (int x = 2; x <= limit; ++x) {
        std::string s = std::to_string(x);
        bool digits_ok = true;
        int c1 = 0;
        int c2 = 0;
        int c3 = 0;
        for (const char ch : s) {
            if (ch == '1') {
                ++c1;
            } else if (ch == '2') {
                ++c2;
            } else if (ch == '3') {
                ++c3;
            } else {
                digits_ok = false;
                break;
            }
        }
        if (!digits_ok) {
            continue;
        }

        if ((c1 == 0 || is123[static_cast<std::size_t>(c1)]) &&
            (c2 == 0 || is123[static_cast<std::size_t>(c2)]) &&
            (c3 == 0 || is123[static_cast<std::size_t>(c3)])) {
            is123[static_cast<std::size_t>(x)] = 1U;
        }
    }

    return is123;
}

class Nth123Solver {
   public:
    Nth123Solver(const u64 cap, const std::vector<std::uint8_t>& is123)
        : cap_(cap), is123_(is123), choose_(1, std::vector<u64>(1, 1U)) {}

    std::string solve(const u64 n) {
        u64 cumulative = 0;
        int length = 0;

        while (true) {
            ++length;
            const u64 count = count_length(length);
            if (sat_add(cumulative, count, cap_) >= n) {
                length_ = length;
                break;
            }
            cumulative += count;
        }

        u64 rank_in_length = n - cumulative;
        memo_.clear();

        int used1 = 0;
        int used2 = 0;
        int used3 = 0;
        std::string out;
        out.reserve(static_cast<std::size_t>(length_));

        for (int pos = 0; pos < length_; ++pos) {
            for (int d = 1; d <= 3; ++d) {
                int n1 = used1;
                int n2 = used2;
                int n3 = used3;
                if (d == 1) {
                    ++n1;
                } else if (d == 2) {
                    ++n2;
                } else {
                    ++n3;
                }

                const u64 cnt = count_with_prefix(pos + 1, n1, n2, n3);
                if (rank_in_length > cnt) {
                    rank_in_length -= cnt;
                } else {
                    out.push_back(static_cast<char>('0' + d));
                    used1 = n1;
                    used2 = n2;
                    used3 = n3;
                    break;
                }
            }
        }

        return out;
    }

   private:
    void ensure_choose(const int n) {
        while (static_cast<int>(choose_.size()) <= n) {
            const int m = static_cast<int>(choose_.size());
            choose_.push_back(std::vector<u64>(static_cast<std::size_t>(m) + 1U, 0U));
            choose_[static_cast<std::size_t>(m)][0] = 1U;
            choose_[static_cast<std::size_t>(m)][static_cast<std::size_t>(m)] = 1U;
            for (int k = 1; k < m; ++k) {
                const u64 left = choose_[static_cast<std::size_t>(m - 1)][static_cast<std::size_t>(k - 1)];
                const u64 right = choose_[static_cast<std::size_t>(m - 1)][static_cast<std::size_t>(k)];
                choose_[static_cast<std::size_t>(m)][static_cast<std::size_t>(k)] =
                    sat_add(left, right, cap_);
            }
        }
    }

    u64 multinomial(const int n, const int a, const int b) {
        ensure_choose(n);
        const u64 first = choose_[static_cast<std::size_t>(n)][static_cast<std::size_t>(a)];
        const u64 second = choose_[static_cast<std::size_t>(n - a)][static_cast<std::size_t>(b)];
        return sat_mul(first, second, cap_);
    }

    u64 count_length(const int length) {
        u64 total = 0;
        for (int a = 0; a <= length; ++a) {
            if (a > 0 && !is123_[static_cast<std::size_t>(a)]) {
                continue;
            }
            for (int b = 0; b + a <= length; ++b) {
                const int c = length - a - b;
                if (b > 0 && !is123_[static_cast<std::size_t>(b)]) {
                    continue;
                }
                if (c > 0 && !is123_[static_cast<std::size_t>(c)]) {
                    continue;
                }
                if (a == 0 && b == 0 && c == 0) {
                    continue;
                }
                total = sat_add(total, multinomial(length, a, b), cap_);
                if (total > cap_) {
                    return cap_ + 1U;
                }
            }
        }
        return total;
    }

    u64 pack_state(const int pos, const int u1, const int u2, const int u3) const {
        return (static_cast<u64>(pos) << 48U) | (static_cast<u64>(u1) << 32U) |
               (static_cast<u64>(u2) << 16U) | static_cast<u64>(u3);
    }

    u64 count_with_prefix(const int pos, const int u1, const int u2, const int u3) {
        const u64 key = pack_state(pos, u1, u2, u3);
        const auto it = memo_.find(key);
        if (it != memo_.end()) {
            return it->second;
        }

        const int rem = length_ - pos;
        u64 total = 0;

        for (int A = u1; A <= length_; ++A) {
            if (A > 0 && !is123_[static_cast<std::size_t>(A)]) {
                continue;
            }
            const int max_b = length_ - A;
            if (u2 > max_b) {
                continue;
            }
            for (int B = u2; B <= max_b; ++B) {
                const int C = length_ - A - B;
                if (C < u3) {
                    continue;
                }
                if (B > 0 && !is123_[static_cast<std::size_t>(B)]) {
                    continue;
                }
                if (C > 0 && !is123_[static_cast<std::size_t>(C)]) {
                    continue;
                }

                const int a = A - u1;
                const int b = B - u2;
                const int c = C - u3;
                if (a + b + c != rem) {
                    continue;
                }

                total = sat_add(total, multinomial(rem, a, b), cap_);
                if (total > cap_) {
                    memo_.emplace(key, cap_ + 1U);
                    return cap_ + 1U;
                }
            }
        }

        memo_.emplace(key, total);
        return total;
    }

    const u64 cap_;
    const std::vector<std::uint8_t>& is123_;
    std::vector<std::vector<u64>> choose_;
    int length_ = 0;
    std::unordered_map<u64, u64> memo_;
};

u64 mod_of_decimal_string(const std::string& s, const u64 mod) {
    u64 value = 0;
    for (const char ch : s) {
        value = (value * 10U + static_cast<u64>(ch - '0')) % mod;
    }
    return value;
}

std::string nth_123_number(const u64 n, const std::vector<std::uint8_t>& is123) {
    Nth123Solver solver(n, is123);
    return solver.solve(n);
}

}  // namespace

int main() {
    const std::vector<std::uint8_t> is123 = build_is123(kIs123Limit);

    assert(nth_123_number(4U, is123) == "11");
    assert(nth_123_number(10U, is123) == "31");
    assert(nth_123_number(40U, is123) == "1112");
    assert(nth_123_number(1000U, is123) == "1223321");
    assert(nth_123_number(6000U, is123) == "2333333333323");

    const std::string value = nth_123_number(kTargetN, is123);
    std::cout << mod_of_decimal_string(value, kMod) << '\n';
    return 0;
}

