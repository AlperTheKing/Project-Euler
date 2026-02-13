#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using boost::multiprecision::cpp_int;

constexpr u64 kMod = 1'000'000'007ULL;
constexpr u64 kTailMod = 100'000'000'000ULL;
constexpr u64 kTailCycleLen = 15'625'000ULL;  // 8 * 5^9
constexpr u64 kTargetN = 10'000'000'000'000'000ULL;

u64 add_mod(u64 x, u64 y) {
    x += y;
    if (x >= kMod) {
        x -= kMod;
    }
    return x;
}

u64 sub_mod(u64 x, u64 y) {
    return x >= y ? x - y : x + kMod - y;
}

u64 mul_mod(u64 x, u64 y, u64 mod) {
    return static_cast<u64>((static_cast<u128>(x) * y) % mod);
}

u64 step_mod(u64 x, u64 mod) {
    return (mul_mod(x, x, mod) + 2ULL) % mod;
}

u64 next_permutation_u64(u64 n) {
    std::string s = std::to_string(n);
    int i = static_cast<int>(s.size()) - 2;
    while (i >= 0 && s[i] >= s[i + 1]) {
        --i;
    }
    if (i < 0) {
        return 0;
    }
    int j = static_cast<int>(s.size()) - 1;
    while (s[j] <= s[i]) {
        --j;
    }
    std::swap(s[i], s[j]);
    std::reverse(s.begin() + i + 1, s.end());
    return static_cast<u64>(std::stoull(s));
}

std::pair<bool, u64> next_perm_delta_11(u64 tail) {
    std::array<int, 11> d{};
    u64 x = tail;
    for (int i = 10; i >= 0; --i) {
        d[i] = static_cast<int>(x % 10ULL);
        x /= 10ULL;
    }

    int i = 9;
    while (i >= 0 && d[i] >= d[i + 1]) {
        --i;
    }
    if (i < 0) {
        return {false, 0};
    }

    int j = 10;
    while (d[j] <= d[i]) {
        --j;
    }
    std::swap(d[i], d[j]);
    std::reverse(d.begin() + i + 1, d.end());

    u64 next_tail = 0;
    for (int digit : d) {
        next_tail = next_tail * 10ULL + static_cast<u64>(digit);
    }
    return {true, next_tail - tail};
}

cpp_int next_permutation_cpp(const cpp_int& value) {
    std::string s = value.convert_to<std::string>();
    int i = static_cast<int>(s.size()) - 2;
    while (i >= 0 && s[i] >= s[i + 1]) {
        --i;
    }
    if (i < 0) {
        return 0;
    }
    int j = static_cast<int>(s.size()) - 1;
    while (s[j] <= s[i]) {
        --j;
    }
    std::swap(s[i], s[j]);
    std::reverse(s.begin() + i + 1, s.end());
    return cpp_int(s);
}

u64 brute_u_mod(int n_limit) {
    cpp_int a = 0;
    cpp_int total = 0;
    for (int n = 1; n <= n_limit; ++n) {
        a = a * a + 2;
        total += next_permutation_cpp(a);
    }
    return static_cast<u64>(total % kMod);
}

struct CycleData {
    u64 mu = 0;
    u64 lambda = 0;
    std::vector<u64> values;
    std::vector<u64> pref_mod;
};

CycleData build_cycle_mod_p() {
    CycleData data;
    std::unordered_map<u64, u64> seen;
    seen.reserve(70000);
    seen.max_load_factor(0.7f);

    u64 x = 0;
    seen.emplace(x, 0);
    data.values.push_back(x);

    for (u64 idx = 1;; ++idx) {
        x = step_mod(x, kMod);
        auto it = seen.find(x);
        if (it != seen.end()) {
            data.mu = it->second;
            data.lambda = idx - it->second;
            break;
        }
        seen.emplace(x, idx);
        data.values.push_back(x);
    }

    data.pref_mod.assign(data.values.size() + 1, 0);
    for (std::size_t i = 0; i < data.values.size(); ++i) {
        data.pref_mod[i + 1] = add_mod(data.pref_mod[i], data.values[i]);
    }
    return data;
}

u64 prefix_sum_cycle(const CycleData& data, u64 n) {
    if (n + 1 <= data.values.size()) {
        return data.pref_mod[static_cast<std::size_t>(n + 1)];
    }

    const u64 mu = data.mu;
    const u64 lambda = data.lambda;
    const u64 sum_mu = data.pref_mod[static_cast<std::size_t>(mu)];
    const u64 cycle_sum = sub_mod(data.pref_mod[static_cast<std::size_t>(mu + lambda)],
                                  data.pref_mod[static_cast<std::size_t>(mu)]);

    const u64 remaining = n - mu + 1;
    const u64 full_cycles = remaining / lambda;
    const u64 rem = remaining % lambda;

    u64 total = sum_mu;
    total = add_mod(total, mul_mod(full_cycles % kMod, cycle_sum, kMod));
    const u64 partial = sub_mod(data.pref_mod[static_cast<std::size_t>(mu + rem)],
                                data.pref_mod[static_cast<std::size_t>(mu)]);
    total = add_mod(total, partial);
    return total;
}

u64 range_sum_cycle(const CycleData& data, u64 left, u64 right) {
    if (right < left) {
        return 0;
    }
    const u64 pref_r = prefix_sum_cycle(data, right);
    const u64 pref_l = (left == 0 ? 0 : prefix_sum_cycle(data, left - 1));
    return sub_mod(pref_r, pref_l);
}

u64 sum_small_exact(u64 n_limit) {
    u64 total = 0;
    u64 a = 0;
    const u64 upto = std::min<u64>(n_limit, 5);
    for (u64 n = 1; n <= upto; ++n) {
        a = a * a + 2;
        total = add_mod(total, next_permutation_u64(a) % kMod);
    }
    return total;
}

u64 sum_a_terms(u64 n_limit, const CycleData& mod_p_cycle) {
    if (n_limit < 6) {
        return 0;
    }
    return range_sum_cycle(mod_p_cycle, 6, n_limit);
}

u64 sum_delta_terms(u64 n_limit) {
    if (n_limit < 6) {
        return 0;
    }

    u64 a5_tail = 0;
    for (int i = 0; i < 5; ++i) {
        a5_tail = step_mod(a5_tail, kTailMod);
    }

    const u64 count = n_limit - 5;
    u64 cycle_sum = 0;
    u64 tail = a5_tail;

    for (u64 i = 0; i < kTailCycleLen; ++i) {
        tail = step_mod(tail, kTailMod);  // n = 6 + i
        const auto [ok, delta] = next_perm_delta_11(tail);
        assert(ok);
        cycle_sum = add_mod(cycle_sum, delta % kMod);
    }
    assert(tail == a5_tail);

    const u64 full_cycles = count / kTailCycleLen;
    const u64 rem = count % kTailCycleLen;

    u64 total = mul_mod(full_cycles % kMod, cycle_sum, kMod);

    tail = a5_tail;
    for (u64 i = 0; i < rem; ++i) {
        tail = step_mod(tail, kTailMod);
        const auto [ok, delta] = next_perm_delta_11(tail);
        assert(ok);
        total = add_mod(total, delta % kMod);
    }

    return total;
}

u64 solve(u64 n_limit) {
    const CycleData mod_p_cycle = build_cycle_mod_p();
    u64 answer = 0;
    answer = add_mod(answer, sum_small_exact(n_limit));
    answer = add_mod(answer, sum_a_terms(n_limit, mod_p_cycle));
    answer = add_mod(answer, sum_delta_terms(n_limit));
    return answer;
}

void validate() {
    assert(brute_u_mod(10) == 543870437ULL);

    u64 a_mod = 0;
    u64 a_tail = 0;
    cpp_int a_exact = 0;
    for (int n = 1; n <= 15; ++n) {
        a_mod = step_mod(a_mod, kMod);
        a_tail = step_mod(a_tail, kTailMod);
        a_exact = a_exact * a_exact + 2;

        const u64 exact_b_mod = static_cast<u64>(next_permutation_cpp(a_exact) % kMod);
        if (n <= 5) {
            assert(exact_b_mod == next_permutation_u64(static_cast<u64>(a_exact)) % kMod);
        } else {
            const auto [ok, delta] = next_perm_delta_11(a_tail);
            assert(ok);
            assert(exact_b_mod == add_mod(a_mod, delta % kMod));
        }
    }
}

}  // namespace

int main() {
    validate();
    std::cout << solve(kTargetN) << '\n';
    return 0;
}
