#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using i64 = long long;

struct State {
    int sum;
    int block;
};

static bool can_reach_pow3_one_step_digits(const int *digits, int len, int max_target,
                                           const std::vector<uint8_t> &is_pow3) {
    static thread_local State buf1[1 << 14], buf2[1 << 14];
    State *cur = buf1;
    State *nxt = buf2;

    int curN = 1;
    cur[0] = {0, digits[0]};

    for (int pos = 1; pos < len; ++pos) {
        const int d = digits[pos];
        int nxtN = 0;
        for (int i = 0; i < curN; ++i) {
            const int s = cur[i].sum;
            const int b = cur[i].block;

            const int s_split = s + b;
            const int b_split = d;
            if (s_split + b_split <= max_target) nxt[nxtN++] = {s_split, b_split};

            const int b_cat = b * 3 + d;
            if (s + b_cat <= max_target) nxt[nxtN++] = {s, b_cat};
        }
        curN = nxtN;
        std::swap(cur, nxt);
    }

    for (int i = 0; i < curN; ++i) {
        const int total = cur[i].sum + cur[i].block;
        if (is_pow3[total]) return true;
    }
    return false;
}

static void init_digits(int x, int base, int *digits, int len, int &sum) {
    sum = 0;
    for (int i = 0; i < len; ++i) digits[i] = 0;
    int pos = 0;
    while (x > 0 && pos < len) {
        digits[pos] = x % base;
        sum += digits[pos];
        x /= base;
        ++pos;
    }
}

static inline void inc_digits(int base, int max_digit, int *digits, int len, int &sum) {
    int pos = 0;
    while (true) {
        if (digits[pos] < max_digit) {
            ++digits[pos];
            ++sum;
            return;
        }
        digits[pos] = 0;
        sum -= max_digit;
        ++pos;
    }
}

static i64 g_worker(int lo, int hi, const std::vector<int> &pow3, const std::vector<uint8_t> &is_pow3) {
    int d10[8], d3[15];
    int s10 = 0, s3 = 0;
    init_digits(lo, 10, d10, 8, s10);
    init_digits(lo, 3, d3, 15, s3);

    int pow_idx = 0;
    while (pow_idx + 1 < (int)pow3.size() && pow3[pow_idx + 1] <= lo) ++pow_idx;
    int max_target = pow3[pow_idx];

    i64 ans = 0;
    for (int i = lo; i <= hi; ++i) {
        while (pow_idx + 1 < (int)pow3.size() && pow3[pow_idx + 1] <= i) {
            ++pow_idx;
            max_target = pow3[pow_idx];
        }

        const int f10 = (i < 10) ? 0 : (s10 < 10 ? 1 : 2);
        int f3 = 0;
        if (i < 3) {
            f3 = 0;
        } else if (s3 < 3) {
            f3 = 1;
        } else if ((i & 1) == 0) {
            f3 = 2;
        } else if (f10 == 2) {
            int ms = 14;
            while (ms > 0 && d3[ms] == 0) --ms;
            int digits[15];
            int len = 0;
            for (int t = ms; t >= 0; --t) digits[len++] = d3[t];
            f3 = can_reach_pow3_one_step_digits(digits, len, max_target, is_pow3) ? 2 : 3;
        } else {
            f3 = 2;
        }

        if (f10 == f3) ans += i;
        if (i == hi) break;
        inc_digits(10, 9, d10, 8, s10);
        inc_digits(3, 2, d3, 15, s3);
    }
    return ans;
}

static i64 g(int n) {
    std::vector<int> pow3;
    for (int x = 3; x <= n; x *= 3) pow3.push_back(x);
    assert(!pow3.empty());

    const int max_pow = pow3.back();
    std::vector<uint8_t> is_pow3(max_pow + 1, 0);
    for (int x : pow3) is_pow3[x] = 1;

    int threads = (int)std::thread::hardware_concurrency();
    if (threads <= 0) threads = 4;
    if (threads > 12) threads = 12;
    if (n < 200'000) threads = 1;

    std::vector<std::thread> ts;
    std::vector<i64> partial(threads, 0);
    ts.reserve(threads);

    const int chunk = (n + threads - 1) / threads;
    for (int t = 0; t < threads; ++t) {
        const int lo = t * chunk + 1;
        const int hi = std::min(n, (t + 1) * chunk);
        if (lo > hi) continue;
        ts.emplace_back([&, t, lo, hi]() { partial[t] = g_worker(lo, hi, pow3, is_pow3); });
    }
    for (auto &th : ts) th.join();
    i64 ans = 0;
    for (i64 x : partial) ans += x;
    return ans;
}

static int f3_exact_small(int n, std::vector<int> &memo) {
    if (n < 3) return 0;
    int &res = memo[n];
    if (res != -1) return res;

    int digits[20];
    int len = 0;
    int x = n;
    while (x > 0) {
        digits[len++] = x % 3;
        x /= 3;
    }
    for (int i = 0; i < len / 2; ++i) std::swap(digits[i], digits[len - 1 - i]);

    int best = 100;
    const int masks = 1 << (len - 1);
    for (int mask = 1; mask < masks; ++mask) {
        int sum = 0;
        int cur = digits[0];
        for (int i = 1; i < len; ++i) {
            if (mask & (1 << (i - 1))) {
                sum += cur;
                cur = digits[i];
            } else {
                cur = cur * 3 + digits[i];
            }
        }
        sum += cur;
        best = std::min(best, 1 + f3_exact_small(sum, memo));
    }
    res = best;
    return res;
}

static int f3_fast_single(int n, const std::vector<int> &pow3, const std::vector<uint8_t> &is_pow3) {
    if (n < 3) return 0;
    int s3 = 0;
    int digits[20];
    int len = 0;
    int x = n;
    while (x > 0) {
        digits[len++] = x % 3;
        s3 += digits[len - 1];
        x /= 3;
    }
    for (int i = 0; i < len / 2; ++i) std::swap(digits[i], digits[len - 1 - i]);
    if (s3 < 3) return 1;
    if ((n & 1) == 0) return 2;

    int max_target = pow3[0];
    for (int t : pow3) {
        if (t > n) break;
        max_target = t;
    }
    return can_reach_pow3_one_step_digits(digits, len, max_target, is_pow3) ? 2 : 3;
}

int main() {
    assert(g(100) == 3302);
    {
        std::vector<int> pow3;
        for (int x = 3; x <= 10'000'000; x *= 3) pow3.push_back(x);
        std::vector<uint8_t> is_pow3(pow3.back() + 1, 0);
        for (int x : pow3) is_pow3[x] = 1;

        std::vector<int> memo(2001, -1);
        for (int i = 1; i <= 2000; ++i) assert(f3_fast_single(i, pow3, is_pow3) == f3_exact_small(i, memo));
    }
    std::cout << g(10'000'000) << "\n";
    return 0;
}
