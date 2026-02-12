#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr u64 A = 25'214'903'917ULL;
static constexpr u64 C = 11ULL;
static constexpr u64 MOD = 1ULL << 48;
static constexpr u64 MASK = MOD - 1ULL;

static inline u64 mul_mod(u64 x, u64 y) {
    return static_cast<u64>((static_cast<u128>(x) * static_cast<u128>(y)) & static_cast<u128>(MASK));
}

static u64 pow_mod(u64 a, u64 e) {
    u64 r = 1ULL;
    while (e > 0) {
        if (e & 1ULL) {
            r = mul_mod(r, a);
        }
        a = mul_mod(a, a);
        e >>= 1ULL;
    }
    return r;
}

static u64 inv_odd(u64 a) {
    u64 x = 1ULL;
    for (int i = 0; i < 6; ++i) {
        const u64 ax = mul_mod(a, x);
        const u64 two_minus_ax = (2ULL + MOD - ax) & MASK;
        x = mul_mod(x, two_minus_ax);
    }
    return x;
}

static inline u64 next_state(u64 a) {
    return (mul_mod(A, a) + C) & MASK;
}

static inline int value_to_code(u64 state) {
    return static_cast<int>((state >> 16) % 52ULL);
}

static int char_code(char ch) {
    if (ch >= 'a' && ch <= 'z') {
        return ch - 'a';
    }
    return ch - 'A' + 26;
}

static std::vector<int> pattern_codes(const std::string& s) {
    std::vector<int> out;
    out.reserve(s.size());
    for (char ch : s) {
        out.push_back(char_code(ch));
    }
    return out;
}

static std::vector<u64> states_starting_with(const std::string& pattern) {
    const std::vector<int> t = pattern_codes(pattern);
    const int n = static_cast<int>(t.size());
    std::vector<int> req_mod4;
    req_mod4.reserve(n - 1);
    for (int i = 0; i + 1 < n; ++i) {
        req_mod4.push_back((t[i + 1] - t[i] + 4) & 3);
    }

    std::vector<u32> low16_candidates;
    low16_candidates.reserve(8);

    for (u32 low0 = 0; low0 < (1u << 16); ++low0) {
        u64 low = low0;
        bool ok = true;
        for (int i = 0; i + 1 < n; ++i) {
            const u64 x = A * low + C;
            const u32 carry = static_cast<u32>(x >> 16);
            if (static_cast<int>(carry & 3u) != req_mod4[i]) {
                ok = false;
                break;
            }
            low = x & 0xFFFFULL;
        }
        if (ok) {
            low16_candidates.push_back(low0);
        }
    }

    std::vector<u64> states;

    for (u32 low0 : low16_candidates) {
        std::vector<u32> carries(n - 1);
        u64 low = low0;
        for (int i = 0; i + 1 < n; ++i) {
            const u64 x = A * low + C;
            carries[i] = static_cast<u32>(x >> 16);
            low = x & 0xFFFFULL;
        }

        for (u64 m = 0;; ++m) {
            const u64 q0_64 = static_cast<u64>(t[0]) + 52ULL * m;
            if (q0_64 > 0xFFFFFFFFULL) {
                break;
            }

            u32 q = static_cast<u32>(q0_64);
            bool ok = true;
            for (int i = 0; i + 1 < n; ++i) {
                q = static_cast<u32>(static_cast<u64>(static_cast<u32>(A)) * static_cast<u64>(q) + static_cast<u64>(carries[i]));
                if (static_cast<int>(q % 52u) != t[i + 1]) {
                    ok = false;
                    break;
                }
            }

            if (ok) {
                states.push_back((q0_64 << 16) | static_cast<u64>(low0));
            }
        }
    }

    std::sort(states.begin(), states.end());
    states.erase(std::unique(states.begin(), states.end()), states.end());
    return states;
}

static u64 advance_state(u64 start, u64 steps) {
    u64 mul = 1ULL;
    u64 add = 0ULL;
    u64 base_mul = A & MASK;
    u64 base_add = C;

    while (steps > 0) {
        if (steps & 1ULL) {
            add = (mul_mod(base_mul, add) + base_add) & MASK;
            mul = mul_mod(base_mul, mul);
        }
        base_add = (mul_mod(base_mul, base_add) + base_add) & MASK;
        base_mul = mul_mod(base_mul, base_mul);
        steps >>= 1ULL;
    }

    return (mul_mod(mul, start) + add) & MASK;
}

static u64 discrete_log_base_A(u64 ratio) {
    static constexpr int ORDER_BITS = 46;
    const u64 inv_a = inv_odd(A & MASK);

    std::array<u64, ORDER_BITS> inv_pows{};
    inv_pows[0] = inv_a;
    for (int i = 1; i < ORDER_BITS; ++i) {
        inv_pows[i] = mul_mod(inv_pows[i - 1], inv_pows[i - 1]);
    }

    u64 cur = ratio;
    u64 exponent = 0ULL;

    for (int i = 0; i < ORDER_BITS; ++i) {
        const u64 test = pow_mod(cur, 1ULL << (ORDER_BITS - 1 - i));
        if (test != 1ULL) {
            exponent |= (1ULL << i);
            cur = mul_mod(cur, inv_pows[i]);
        }
    }

    assert(cur == 1ULL);
    return exponent;
}

static u64 distance_between_states(u64 start, u64 target) {
    const u64 b_start = (mul_mod((A - 1ULL) & MASK, start) + C) & MASK;
    const u64 b_target = (mul_mod((A - 1ULL) & MASK, target) + C) & MASK;

    const u64 inv_b_start = inv_odd(b_start);
    const u64 ratio = mul_mod(b_target, inv_b_start);
    const u64 base_distance = discrete_log_base_A(ratio);

    static constexpr u64 ORDER = 1ULL << 46;
    for (u64 k = 0; k < 4ULL; ++k) {
        const u64 candidate = base_distance + k * ORDER;
        if (advance_state(start, candidate) == target) {
            return candidate;
        }
    }

    assert(false);
    return 0ULL;
}

static u64 first_occurrence_bruteforce(u64 seed, const std::string& pattern, u64 limit) {
    const std::vector<int> t = pattern_codes(pattern);
    const int len = static_cast<int>(t.size());

    std::vector<int> window;
    window.reserve(len);

    u64 state = seed;
    for (u64 i = 0; i < limit; ++i) {
        const int c = value_to_code(state);
        state = next_state(state);

        if (static_cast<int>(window.size()) < len) {
            window.push_back(c);
        } else {
            for (int j = 0; j + 1 < len; ++j) {
                window[j] = window[j + 1];
            }
            window[len - 1] = c;
        }

        if (static_cast<int>(window.size()) == len) {
            bool ok = true;
            for (int j = 0; j < len; ++j) {
                if (window[j] != t[j]) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                return i + 1ULL - static_cast<u64>(len);
            }
        }
    }

    return std::numeric_limits<u64>::max();
}

int main() {
    {
        const auto states = states_starting_with("EULERcats");
        assert(states.size() == 1);
        assert(states[0] == 78'580'612'777'175ULL);
    }

    assert(first_occurrence_bruteforce(123456ULL, "RxqLBfWzv", 5000ULL) == 100ULL);

    const auto puzzle_states = states_starting_with("PuzzleOne");
    const auto lucky_states = states_starting_with("LuckyText");
    assert(puzzle_states.size() == 1);
    assert(lucky_states.size() == 1);

    const u64 answer = distance_between_states(puzzle_states[0], lucky_states[0]);
    std::cout << answer << '\n';
    return 0;
}
