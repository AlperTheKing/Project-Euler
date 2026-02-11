#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;

struct ResidueEnumerator {
    int mod;
    std::vector<int> w;
    std::vector<int> nine_w;
    std::vector<std::uint8_t> digits;
    int res;
    bool done;

    ResidueEnumerator(int mod_value, const std::vector<int>& weights)
        : mod(mod_value), w(weights), nine_w(weights.size(), 0), digits(weights.size(), 0), res(0), done(false) {
        for (std::size_t i = 0; i < w.size(); ++i) {
            nine_w[i] = static_cast<int>((9LL * w[i]) % mod);
        }
        if (w.empty()) done = false;
    }

    bool next(int& out) {
        if (done) return false;
        out = res;

        if (w.empty()) {
            done = true;
            return true;
        }

        std::size_t p = 0;
        while (p < w.size()) {
            if (digits[p] < 9) {
                ++digits[p];
                res += w[p];
                if (res >= mod) res -= mod;
                return true;
            }
            digits[p] = 0;
            res -= nine_w[p];
            if (res < 0) res += mod;
            ++p;
        }

        done = true;
        return true;
    }
};

i64 count_len_pal_divisible(int len,
                            int mod,
                            const std::vector<int>& pow10_mod,
                            std::vector<int>& count_res,
                            std::vector<int>& touched,
                            std::vector<int>& rest_residues) {
    const int h = (len + 1) / 2;
    const int split = (h + 1) / 2;

    std::vector<int> weight(static_cast<std::size_t>(h), 0);
    for (int i = 0; i < h; ++i) {
        const int lo = i;
        const int hi = len - 1 - i;
        if (lo == hi) {
            weight[static_cast<std::size_t>(i)] = pow10_mod[static_cast<std::size_t>(lo)];
        } else {
            int v = pow10_mod[static_cast<std::size_t>(lo)] + pow10_mod[static_cast<std::size_t>(hi)];
            if (v >= mod) v -= mod;
            weight[static_cast<std::size_t>(i)] = v;
        }
    }

    const int left_w0 = weight[0];
    std::vector<int> left_rest_w;
    left_rest_w.reserve(static_cast<std::size_t>(std::max(0, split - 1)));
    for (int i = 1; i < split; ++i) left_rest_w.push_back(weight[static_cast<std::size_t>(i)]);

    rest_residues.clear();
    {
        ResidueEnumerator enum_rest(mod, left_rest_w);
        int r = 0;
        while (enum_rest.next(r)) rest_residues.push_back(r);
    }

    touched.clear();
    for (int lead = 1; lead <= 9; ++lead) {
        const int lead_contrib = static_cast<int>((1LL * lead * left_w0) % mod);
        for (int rr : rest_residues) {
            int v = lead_contrib + rr;
            if (v >= mod) v -= mod;
            if (count_res[static_cast<std::size_t>(v)] == 0) touched.push_back(v);
            ++count_res[static_cast<std::size_t>(v)];
        }
    }

    std::vector<int> right_w;
    right_w.reserve(static_cast<std::size_t>(h - split));
    for (int i = split; i < h; ++i) right_w.push_back(weight[static_cast<std::size_t>(i)]);

    i64 total = 0;
    {
        ResidueEnumerator enum_right(mod, right_w);
        int rr = 0;
        while (enum_right.next(rr)) {
            const int need = (rr == 0) ? 0 : (mod - rr);
            total += count_res[static_cast<std::size_t>(need)];
        }
    }

    for (int idx : touched) count_res[static_cast<std::size_t>(idx)] = 0;
    return total;
}

i64 solve_case(int mod, int max_len) {
    std::vector<int> pow10_mod(static_cast<std::size_t>(max_len), 0);
    pow10_mod[0] = 1 % mod;
    for (int i = 1; i < max_len; ++i) {
        pow10_mod[static_cast<std::size_t>(i)] =
            static_cast<int>((10LL * pow10_mod[static_cast<std::size_t>(i - 1)]) % mod);
    }

    std::vector<int> count_res(static_cast<std::size_t>(mod), 0);
    std::vector<int> touched;
    touched.reserve(static_cast<std::size_t>(mod));
    std::vector<int> rest_residues;

    i64 answer = 0;
    for (int len = 1; len <= max_len; ++len) {
        answer += count_len_pal_divisible(len, mod, pow10_mod, count_res, touched, rest_residues);
    }
    return answer;
}

}  // namespace

int main() {
    assert(solve_case(109, 5) == 9);
    std::cout << solve_case(10'000'019, 32) << "\n";
    return 0;
}
