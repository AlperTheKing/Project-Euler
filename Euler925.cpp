#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 1'000'000'007ULL;

u64 add_mod(u64 a, u64 b) {
    a += b;
    if (a >= kMod) {
        a -= kMod;
    }
    return a;
}

u64 sub_mod(u64 a, u64 b) {
    return (a >= b) ? (a - b) : (a + kMod - b);
}

u64 mul_mod(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * b) % kMod);
}

bool next_perm_digits(std::uint8_t* digits, int length) {
    int i = length - 2;
    while (i >= 0 && digits[i] >= digits[i + 1]) {
        --i;
    }
    if (i < 0) {
        return false;
    }
    int j = length - 1;
    while (digits[j] <= digits[i]) {
        --j;
    }
    std::swap(digits[i], digits[j]);
    std::reverse(digits + i + 1, digits + length);
    return true;
}

u64 b_of_square_mod(u64 n) {
    u128 x = static_cast<u128>(n) * n;
    if (x == 0) {
        return 0;
    }

    std::uint8_t rev[40];
    int length = 0;
    while (x > 0) {
        rev[length++] = static_cast<std::uint8_t>(x % 10);
        x /= 10;
    }

    std::uint8_t digits[40];
    for (int i = 0; i < length; ++i) {
        digits[i] = rev[length - 1 - i];
    }
    if (!next_perm_digits(digits, length)) {
        return 0;
    }

    u64 value_mod = 0;
    for (int i = 0; i < length; ++i) {
        value_mod = (value_mod * 10 + digits[i]) % kMod;
    }
    return value_mod;
}

u64 brute_sum(u64 n) {
    u64 total = 0;
    for (u64 i = 1; i <= n; ++i) {
        total = add_mod(total, b_of_square_mod(i));
    }
    return total;
}

struct Prefix {
    u64 a;
    u64 n_mod;
    u64 sq_low;
    int depth;
    int prev;
};

struct Frame {
    u64 a;
    u64 n_mod;
    u64 sq_low;
    int depth;
    int prev;
    int next_digit;
};

struct RangePre {
    bool split_by_a = false;
    u64 q_mod = 0;
    u64 q2_mod = 0;
    u64 cnt0_mod = 0;
    u64 cnt1_mod = 0;
    u64 s1_mod = 0;
    u64 s2_mod = 0;
};

class Solver {
  public:
    explicit Solver(int digits) : digits_(digits), split_(std::min(8, digits)) {
        pow10_.assign(digits_ + 1, 1);
        for (int i = 1; i <= digits_; ++i) {
            pow10_[i] = pow10_[i - 1] * 10ULL;
        }

        a_exp_ = (digits_ + 1) / 2;
        a_min_ = pow10_[a_exp_];

        pow10_mod_.assign(2 * digits_ + 2, 1);
        for (int i = 1; i < static_cast<int>(pow10_mod_.size()); ++i) {
            pow10_mod_[i] = mul_mod(pow10_mod_[i - 1], 10);
        }

        add_digit_mod_.assign(digits_ + 1, {});
        for (int k = 0; k <= digits_; ++k) {
            for (int d = 0; d <= 9; ++d) {
                add_digit_mod_[k][d] = mul_mod(static_cast<u64>(d), pow10_mod_[k]);
            }
        }

        pre_.assign(digits_ + 1, {});
        precompute_ranges();
    }

    u64 solve() {
        const u64 small = brute_small_parallel();
        build_prefixes();

        unsigned threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 8;
        }

        std::atomic<std::size_t> index{0};
        std::vector<u64> local(threads, 0);
        std::vector<std::thread> workers;
        workers.reserve(threads);

        for (unsigned tid = 0; tid < threads; ++tid) {
            workers.emplace_back([&, tid]() {
                u64 subtotal = 0;
                while (true) {
                    const std::size_t i = index.fetch_add(1, std::memory_order_relaxed);
                    if (i >= prefixes_.size()) {
                        break;
                    }
                    subtotal = add_mod(subtotal, process_prefix(prefixes_[i]));
                }
                local[tid] = subtotal;
            });
        }
        for (std::thread& t : workers) {
            t.join();
        }

        u64 answer = add_mod(small, base_sum_);
        for (u64 v : local) {
            answer = add_mod(answer, v);
        }
        return answer;
    }

  private:
    int digits_;
    int split_;
    int a_exp_;
    u64 a_min_;

    std::vector<u64> pow10_;
    std::vector<u64> pow10_mod_;
    std::vector<std::array<u64, 10>> add_digit_mod_;
    std::vector<RangePre> pre_;

    std::vector<Prefix> prefixes_;
    u64 base_sum_ = 0;

    u64 sum0(u64 n) const {
        const u64 n_mod = n % kMod;
        const u64 n1_mod = (n + 1) % kMod;
        constexpr u64 inv2 = (kMod + 1) / 2;
        return mul_mod(mul_mod(n_mod, n1_mod), inv2);
    }

    u64 sum2_0(u64 n) const {
        const u64 n_mod = n % kMod;
        const u64 n1_mod = (n + 1) % kMod;
        const u64 n2p1_mod = (2ULL * n_mod + 1) % kMod;
        constexpr u64 inv6 = 166666668ULL;
        return mul_mod(mul_mod(mul_mod(n_mod, n1_mod), n2p1_mod), inv6);
    }

    u64 sum_range_1(u64 l, u64 r) const {
        if (l > r) {
            return 0;
        }
        return sub_mod(sum0(r), (l == 0 ? 0 : sum0(l - 1)));
    }

    u64 sum_range_2(u64 l, u64 r) const {
        if (l > r) {
            return 0;
        }
        return sub_mod(sum2_0(r), (l == 0 ? 0 : sum2_0(l - 1)));
    }

    void precompute_ranges() {
        for (int r = 1; r <= digits_; ++r) {
            const u64 block_count = pow10_[digits_ - r];
            RangePre cur;
            cur.q_mod = pow10_mod_[r];
            cur.q2_mod = mul_mod(cur.q_mod, cur.q_mod);

            if (r <= a_exp_) {
                const u64 t0 = pow10_[a_exp_ - r];
                const u64 count = block_count - t0;
                cur.split_by_a = false;
                cur.cnt0_mod = count % kMod;
                cur.s1_mod = sum_range_1(t0, block_count - 1);
                cur.s2_mod = sum_range_2(t0, block_count - 1);
            } else {
                cur.split_by_a = true;
                cur.cnt0_mod = block_count % kMod;
                cur.cnt1_mod = (block_count - 1) % kMod;
                cur.s1_mod = sum0(block_count - 1);
                cur.s2_mod = sum2_0(block_count - 1);
            }

            pre_[r] = cur;
        }
    }

    u64 delta_fixed(u64 low, int r) const {
        const u64 original = low;
        std::uint8_t digits[40];
        u64 x = low;
        for (int i = r - 1; i >= 0; --i) {
            digits[i] = static_cast<std::uint8_t>(x % 10ULL);
            x /= 10ULL;
        }

        const bool ok = next_perm_digits(digits, r);
        assert(ok);

        u64 next_mod = 0;
        for (int i = 0; i < r; ++i) {
            next_mod = (next_mod * 10 + digits[i]) % kMod;
        }
        return sub_mod(next_mod, original % kMod);
    }

    u64 resolved_contrib(int r, u64 a, u64 a_mod, u64 low_r) const {
        const RangePre& pr = pre_[r];
        const u64 delta = delta_fixed(low_r, r);

        const u64 a2 = mul_mod(a_mod, a_mod);
        const u64 two_a_q = mul_mod((2ULL * a_mod) % kMod, pr.q_mod);

        auto build_sumsq = [&](u64 cnt_mod) {
            u64 res = 0;
            res = add_mod(res, mul_mod(cnt_mod, a2));
            res = add_mod(res, mul_mod(two_a_q, pr.s1_mod));
            res = add_mod(res, mul_mod(pr.q2_mod, pr.s2_mod));
            return res;
        };

        u64 cnt_mod = 0;
        u64 sumsq_mod = 0;

        if (!pr.split_by_a) {
            cnt_mod = pr.cnt0_mod;
            sumsq_mod = build_sumsq(cnt_mod);
        } else {
            const u64 full = build_sumsq(pr.cnt0_mod);
            if (a < a_min_) {
                cnt_mod = pr.cnt1_mod;
                sumsq_mod = sub_mod(full, a2);
            } else {
                cnt_mod = pr.cnt0_mod;
                sumsq_mod = full;
            }
        }

        return add_mod(sumsq_mod, mul_mod(cnt_mod, delta));
    }

    u64 leaf_exact(const Frame& f) const {
        if (f.a < a_min_) {
            return 0;
        }

        const u64 n2_mod = mul_mod(f.n_mod, f.n_mod);
        const u128 n2 = static_cast<u128>(f.a) * f.a;

        u128 high = n2 / pow10_[digits_];
        int prev_digit = f.prev;

        std::uint8_t extras[40];
        int extra_count = 0;
        bool found = false;

        while (high > 0) {
            const int digit = static_cast<int>(high % 10);
            high /= 10;

            extras[extra_count++] = static_cast<std::uint8_t>(digit);
            if (digit < prev_digit) {
                found = true;
                break;
            }
            prev_digit = digit;
        }

        if (!found) {
            return 0;
        }

        const int suffix_len = digits_ + extra_count;
        std::uint8_t lsd[40];

        u64 tmp = f.sq_low;
        for (int i = 0; i < digits_; ++i) {
            lsd[i] = static_cast<std::uint8_t>(tmp % 10ULL);
            tmp /= 10ULL;
        }
        for (int i = 0; i < extra_count; ++i) {
            lsd[digits_ + i] = extras[i];
        }

        u64 original_mod = f.sq_low % kMod;
        for (int i = 0; i < extra_count; ++i) {
            original_mod = add_mod(original_mod, mul_mod(extras[i], pow10_mod_[digits_ + i]));
        }

        std::uint8_t msd[40];
        for (int i = 0; i < suffix_len; ++i) {
            msd[i] = lsd[suffix_len - 1 - i];
        }

        const bool ok = next_perm_digits(msd, suffix_len);
        assert(ok);

        u64 perm_mod = 0;
        for (int i = 0; i < suffix_len; ++i) {
            perm_mod = (perm_mod * 10 + msd[i]) % kMod;
        }

        const u64 delta = sub_mod(perm_mod, original_mod);
        return add_mod(n2_mod, delta);
    }

    void build_prefixes() {
        prefixes_.clear();
        base_sum_ = 0;

        std::vector<Frame> st;
        st.reserve(64);
        st.push_back({0, 0, 0, 0, 0, 0});

        while (!st.empty()) {
            Frame& f = st.back();

            if (f.depth == split_) {
                prefixes_.push_back({f.a, f.n_mod, f.sq_low, f.depth, f.prev});
                st.pop_back();
                continue;
            }
            if (f.next_digit == 10) {
                st.pop_back();
                continue;
            }

            const int d = f.next_digit++;
            const u64 a2 = f.a + static_cast<u64>(d) * pow10_[f.depth];
            const u64 n_mod2 = add_mod(f.n_mod, add_digit_mod_[f.depth][d]);

            const u64 mod = pow10_[f.depth + 1];
            const u64 y = static_cast<u64>((static_cast<u128>(a2) * a2) % mod);
            const int dig = static_cast<int>(y / pow10_[f.depth]);

            if (f.depth == 0 || dig >= f.prev) {
                st.push_back({a2, n_mod2, y, f.depth + 1, dig, 0});
            } else {
                base_sum_ = add_mod(base_sum_, resolved_contrib(f.depth + 1, a2, n_mod2, y));
            }
        }
    }

    u64 process_prefix(const Prefix& p) const {
        std::vector<Frame> st;
        st.reserve(64);
        st.push_back({p.a, p.n_mod, p.sq_low, p.depth, p.prev, 0});

        u64 total = 0;

        while (!st.empty()) {
            Frame& f = st.back();

            if (f.depth == digits_) {
                total = add_mod(total, leaf_exact(f));
                st.pop_back();
                continue;
            }
            if (f.next_digit == 10) {
                st.pop_back();
                continue;
            }

            const int d = f.next_digit++;
            const u64 a2 = f.a + static_cast<u64>(d) * pow10_[f.depth];
            const u64 n_mod2 = add_mod(f.n_mod, add_digit_mod_[f.depth][d]);

            const u64 mod = pow10_[f.depth + 1];
            const u64 y = static_cast<u64>((static_cast<u128>(a2) * a2) % mod);
            const int dig = static_cast<int>(y / pow10_[f.depth]);

            if (dig >= f.prev) {
                st.push_back({a2, n_mod2, y, f.depth + 1, dig, 0});
            } else {
                total = add_mod(total, resolved_contrib(f.depth + 1, a2, n_mod2, y));
            }
        }

        return total;
    }

    u64 brute_small_parallel() const {
        const u64 upto = a_min_ - 1;
        if (upto == 0) {
            return 0;
        }

        unsigned threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 8;
        }

        const u64 chunk = (upto + threads - 1) / threads;
        std::vector<u64> local(threads, 0);
        std::vector<std::thread> th;
        th.reserve(threads);

        for (unsigned tid = 0; tid < threads; ++tid) {
            const u64 l = 1 + static_cast<u64>(tid) * chunk;
            const u64 r = std::min<u64>(upto, static_cast<u64>(tid + 1) * chunk);
            th.emplace_back([&, tid, l, r]() {
                if (l > r) {
                    return;
                }
                u64 subtotal = 0;
                for (u64 n = l; n <= r; ++n) {
                    subtotal = add_mod(subtotal, b_of_square_mod(n));
                }
                local[tid] = subtotal;
            });
        }
        for (std::thread& t : th) {
            t.join();
        }

        u64 total = 0;
        for (u64 v : local) {
            total = add_mod(total, v);
        }
        return total;
    }
};

u64 solve_power10(int d) {
    Solver solver(d);
    return solver.solve();
}

void validate() {
    assert(brute_sum(10) == 270);
    assert(brute_sum(100) == 335316);

    const u64 fast_d6 = solve_power10(6);
    const u64 brute_d6 = brute_sum(1'000'000);
    assert(fast_d6 == brute_d6);
}

}  // namespace

int main() {
    validate();
    std::cout << solve_power10(16) << '\n';
    return 0;
}
