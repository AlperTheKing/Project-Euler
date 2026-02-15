#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <unordered_set>
#include <unistd.h>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct State {
    u64 a;
    u64 b;
    u64 c;
    u32 mask;
};

struct Tail {
    int value;
    u32 mask;
};

struct TailTriple {
    int a;
    int b;
    int c;
};

bool feasible(u64 a, u64 b, u64 c, u64 p) {
    const u64 min_a = a * p;
    const u64 min_b = b * p;
    const u64 min_c = c * p;
    const u64 max_a = min_a + p - 1;
    const u64 max_b = min_b + p - 1;
    const u64 max_c = min_c + p - 1;

    if (min_c >= max_a + max_b) return false;

    const u128 max_left = static_cast<u128>(max_a) * max_a +
                          static_cast<u128>(max_b) * max_b +
                          static_cast<u128>(max_a) * max_b;
    const u128 min_left = static_cast<u128>(min_a) * min_a +
                          static_cast<u128>(min_b) * min_b +
                          static_cast<u128>(min_a) * min_b;
    const u128 min_right = static_cast<u128>(min_c) * min_c;
    const u128 max_right = static_cast<u128>(max_c) * max_c;

    if (max_left < min_right) return false;
    if (min_left > max_right) return false;
    return true;
}

std::vector<int> digits_from_mask(u32 mask, int base) {
    std::vector<int> digits;
    digits.reserve(static_cast<std::size_t>(base));
    for (int d = 0; d < base; ++d) {
        if ((mask & (1u << d)) == 0) digits.push_back(d);
    }
    return digits;
}

std::vector<State> initial_states(int base, u64 p) {
    std::vector<State> states;
    const int r = base % 3;

    if (r == 0) {
        for (int a = 1; a < base; ++a) {
            for (int b = a + 1; b < base; ++b) {
                for (int c = b + 1; c < base; ++c) {
                    u32 mask = (1u << a) | (1u << b) | (1u << c);
                    if (feasible(a, b, c, p)) {
                        states.push_back({static_cast<u64>(a), static_cast<u64>(b), static_cast<u64>(c), mask});
                    }
                }
            }
        }
        return states;
    }

    if (r == 1) {
        for (int a = 1; a < base; ++a) {
            for (int b = a + 1; b < base; ++b) {
                u32 mask_ab = (1u << a) | (1u << b);
                for (int c1 = 1; c1 < base; ++c1) {
                    if (mask_ab & (1u << c1)) continue;
                    for (int c0 = 0; c0 < base; ++c0) {
                        if (c0 == c1) continue;
                        if (mask_ab & (1u << c0)) continue;
                        const int c = c1 * base + c0;
                        const u32 mask = mask_ab | (1u << c1) | (1u << c0);
                        if (feasible(a, b, c, p)) {
                            states.push_back({static_cast<u64>(a), static_cast<u64>(b), static_cast<u64>(c), mask});
                        }
                    }
                }
            }
        }
        return states;
    }

    for (int a = 1; a < base; ++a) {
        const u32 mask_a = (1u << a);
        for (int b1 = 1; b1 < base; ++b1) {
            if (mask_a & (1u << b1)) continue;
            for (int b0 = 0; b0 < base; ++b0) {
                if (b0 == b1) continue;
                u32 mask_ab = mask_a | (1u << b1) | (1u << b0);
                if (__builtin_popcount(mask_ab) != 3) continue;
                const int b = b1 * base + b0;
                for (int c1 = 1; c1 < base; ++c1) {
                    if (mask_ab & (1u << c1)) continue;
                    for (int c0 = 0; c0 < base; ++c0) {
                        if (c0 == c1) continue;
                        if (mask_ab & (1u << c0)) continue;
                        const int c = c1 * base + c0;
                        if (b >= c) continue;
                        const u32 mask = mask_ab | (1u << c1) | (1u << c0);
                        if (feasible(a, b, c, p)) {
                            states.push_back({static_cast<u64>(a), static_cast<u64>(b), static_cast<u64>(c), mask});
                        }
                    }
                }
            }
        }
    }
    return states;
}

std::vector<State> extend_states(const std::vector<State>& states, int base, u64 p) {
    std::vector<State> next;
    for (const auto& s : states) {
        const auto digits = digits_from_mask(s.mask, base);
        const int m = static_cast<int>(digits.size());
        for (int ia = 0; ia < m; ++ia) {
            const int da = digits[ia];
            for (int ib = 0; ib < m; ++ib) {
                if (ib == ia) continue;
                const int db = digits[ib];
                for (int ic = 0; ic < m; ++ic) {
                    if (ic == ia || ic == ib) continue;
                    const int dc = digits[ic];
                    const u32 mask = s.mask | (1u << da) | (1u << db) | (1u << dc);
                    const u64 na = s.a * base + static_cast<u64>(da);
                    const u64 nb = s.b * base + static_cast<u64>(db);
                    const u64 nc = s.c * base + static_cast<u64>(dc);
                    if (feasible(na, nb, nc, p)) {
                        next.push_back({na, nb, nc, mask});
                    }
                }
            }
        }
    }
    return next;
}

u64 finalize_two_digits(const std::vector<State>& states, int base) {
    const u64 pow2 = static_cast<u64>(base) * static_cast<u64>(base);
    u64 sum = 0;
    for (const auto& s : states) {
        auto digits = digits_from_mask(s.mask, base);
        if (digits.size() != 6) continue;
        std::sort(digits.begin(), digits.end());
        do {
            const u64 a = s.a * pow2 + static_cast<u64>(digits[0] * base + digits[1]);
            const u64 b = s.b * pow2 + static_cast<u64>(digits[2] * base + digits[3]);
            const u64 c = s.c * pow2 + static_cast<u64>(digits[4] * base + digits[5]);
            const u128 left = static_cast<u128>(a) * a + static_cast<u128>(b) * b + static_cast<u128>(a) * b;
            const u128 right = static_cast<u128>(c) * c;
            if (left == right && c < a + b) sum += c;
        } while (std::next_permutation(digits.begin(), digits.end()));
    }
    return sum;
}

u64 finalize_three_digits(const std::vector<State>& states, int base) {
    const int mod = base * base * base;
    const u64 pow3 = static_cast<u64>(base) * base * base;
    std::vector<Tail> tails;
    tails.reserve(static_cast<std::size_t>(base * (base - 1) * (base - 2)));
    for (int d2 = 0; d2 < base; ++d2) {
        for (int d1 = 0; d1 < base; ++d1) {
            if (d1 == d2) continue;
            for (int d0 = 0; d0 < base; ++d0) {
                if (d0 == d1 || d0 == d2) continue;
                const int value = d2 * base * base + d1 * base + d0;
                const u32 mask = (1u << d2) | (1u << d1) | (1u << d0);
                tails.push_back({value, mask});
            }
        }
    }

    std::vector<std::vector<int>> c_by_res(static_cast<std::size_t>(mod));
    for (std::size_t i = 0; i < tails.size(); ++i) {
        const u64 v = static_cast<u64>(tails[i].value);
        const int res = static_cast<int>((v * v) % static_cast<u64>(mod));
        c_by_res[static_cast<std::size_t>(res)].push_back(static_cast<int>(i));
    }

    std::vector<std::vector<TailTriple>> by_mask(1u << base);
    for (const auto& a : tails) {
        for (const auto& b : tails) {
            if (a.mask & b.mask) continue;
            const u64 av = static_cast<u64>(a.value);
            const u64 bv = static_cast<u64>(b.value);
            const int res = static_cast<int>((av * av + bv * bv + av * bv) % static_cast<u64>(mod));
            const u32 ab_mask = a.mask | b.mask;
            for (int idx : c_by_res[static_cast<std::size_t>(res)]) {
                const auto& c = tails[static_cast<std::size_t>(idx)];
                if (ab_mask & c.mask) continue;
                const u32 mask = ab_mask | c.mask;
                by_mask[mask].push_back({a.value, b.value, c.value});
            }
        }
    }

    const u32 full_mask = (1u << base) - 1u;
    u64 sum = 0;
    for (const auto& s : states) {
        const u32 remaining = full_mask ^ s.mask;
        for (const auto& t : by_mask[remaining]) {
            const u64 a = s.a * pow3 + static_cast<u64>(t.a);
            const u64 b = s.b * pow3 + static_cast<u64>(t.b);
            const u64 c = s.c * pow3 + static_cast<u64>(t.c);
            const u128 left = static_cast<u128>(a) * a + static_cast<u128>(b) * b + static_cast<u128>(a) * b;
            const u128 right = static_cast<u128>(c) * c;
            if (left == right && c < a + b) sum += c;
        }
    }
    return sum;
}

u64 solve_base(int base) {
    const int k = base / 3;
    const int remaining = k - 1;
    std::vector<u64> pow_base(static_cast<std::size_t>(k + 1), 1);
    for (int i = 1; i <= k; ++i) pow_base[static_cast<std::size_t>(i)] = pow_base[static_cast<std::size_t>(i - 1)] * base;
    auto states = initial_states(base, pow_base[static_cast<std::size_t>(remaining)]);
    if (remaining > 3) {
        for (int rem = remaining - 1; rem >= 3; --rem) {
            states = extend_states(states, base, pow_base[static_cast<std::size_t>(rem)]);
        }
    }
    if (remaining == 2) return finalize_two_digits(states, base);
    return finalize_three_digits(states, base);
}

u64 brute_base9() {
    std::array<int, 9> digits{};
    for (int i = 0; i < 9; ++i) digits[static_cast<std::size_t>(i)] = i;
    struct TripleKey {
        u64 a;
        u64 b;
        u64 c;
    };
    struct Hash {
        std::size_t operator()(const TripleKey& t) const noexcept {
            std::size_t h = static_cast<std::size_t>(t.a);
            h ^= static_cast<std::size_t>(t.b) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            h ^= static_cast<std::size_t>(t.c) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            return h;
        }
    };
    struct Eq {
        bool operator()(const TripleKey& x, const TripleKey& y) const noexcept {
            return x.a == y.a && x.b == y.b && x.c == y.c;
        }
    };

    std::unordered_set<TripleKey, Hash, Eq> seen;
    u64 sum = 0;
    do {
        for (int i = 1; i <= 7; ++i) {
            if (digits[0] == 0) continue;
            for (int j = i + 1; j <= 8; ++j) {
                if (digits[static_cast<std::size_t>(i)] == 0 || digits[static_cast<std::size_t>(j)] == 0) continue;
                u64 a = 0;
                for (int x = 0; x < i; ++x) a = a * 9 + static_cast<u64>(digits[static_cast<std::size_t>(x)]);
                u64 b = 0;
                for (int x = i; x < j; ++x) b = b * 9 + static_cast<u64>(digits[static_cast<std::size_t>(x)]);
                u64 c = 0;
                for (int x = j; x < 9; ++x) c = c * 9 + static_cast<u64>(digits[static_cast<std::size_t>(x)]);

                std::array<u64, 3> s{a, b, c};
                std::sort(s.begin(), s.end());
                a = s[0];
                b = s[1];
                c = s[2];
                const u128 left = static_cast<u128>(a) * a + static_cast<u128>(b) * b + static_cast<u128>(a) * b;
                const u128 right = static_cast<u128>(c) * c;
                if (left != right || c >= a + b) continue;
                const TripleKey key{a, b, c};
                if (seen.insert(key).second) sum += c;
            }
        }
    } while (std::next_permutation(digits.begin(), digits.end()));
    return sum;
}

struct WorkerData {
    std::atomic<int>* next_base;
    std::vector<u64>* results;
};

void* worker_fn(void* arg) {
    auto* data = static_cast<WorkerData*>(arg);
    while (true) {
        const int base = data->next_base->fetch_add(1);
        if (base > 18) break;
        (*data->results)[static_cast<std::size_t>(base)] = solve_base(base);
    }
    return nullptr;
}

}  // namespace

int main() {
    std::vector<u64> results(19, 0);
    std::atomic<int> next_base{9};

    long threads = sysconf(_SC_NPROCESSORS_ONLN);
    if (threads < 1) threads = 1;
    if (threads > 10) threads = 10;

    std::vector<pthread_t> pool(static_cast<std::size_t>(threads));
    WorkerData data{&next_base, &results};
    for (long i = 0; i < threads; ++i) {
        pthread_create(&pool[static_cast<std::size_t>(i)], nullptr, worker_fn, &data);
    }
    for (long i = 0; i < threads; ++i) {
        pthread_join(pool[static_cast<std::size_t>(i)], nullptr);
    }

    const u64 brute9 = brute_base9();
    assert(results[9] == brute9);

    u64 sum = 0;
    for (int base = 9; base <= 18; ++base) sum += results[static_cast<std::size_t>(base)];
    std::cout << sum << '\n';
    return 0;
}
