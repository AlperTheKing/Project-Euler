#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

static constexpr long long MOD = 989898989LL;

static long long norm_mod(long long x) {
    x %= MOD;
    if (x < 0) x += MOD;
    return x;
}

static long long add_mod(long long a, long long b) {
    long long v = a + b;
    if (v >= MOD) v -= MOD;
    return v;
}

static long long sub_mod(long long a, long long b) {
    long long v = a - b;
    if (v < 0) v += MOD;
    return v;
}

static long long mul_mod(long long a, long long b) {
    return static_cast<long long>(
        (static_cast<unsigned long long>(a) * static_cast<unsigned long long>(b)) %
        static_cast<unsigned long long>(MOD));
}

static long long mod_pow(long long base, long long exp) {
    long long result = 1;
    long long cur = norm_mod(base);
    while (exp > 0) {
        if (exp & 1LL) result = mul_mod(result, cur);
        cur = mul_mod(cur, cur);
        exp >>= 1LL;
    }
    return result;
}

static long long mod_inv(long long a) {
    long long t = 0, new_t = 1;
    long long r = MOD, new_r = norm_mod(a);
    while (new_r != 0) {
        long long q = r / new_r;
        long long tmp_t = t - q * new_t;
        t = new_t;
        new_t = tmp_t;
        long long tmp_r = r - q * new_r;
        r = new_r;
        new_r = tmp_r;
    }
    if (r != 1) {
        throw runtime_error("mod inverse does not exist");
    }
    return norm_mod(t);
}

static long long g8(long long n, long long inv3) {
    if (n <= 0) return 0;
    if (n == 1) return 1;
    long long coeff = (n % 2 == 0) ? 14 : 20;
    long long term1 = mul_mod(coeff, mul_mod(mod_pow(2, n / 2), inv3));
    long long term2 = norm_mod(6 * n - 15);
    return add_mod(term1, term2);
}

static long long compute_G(int n) {
    if (n <= 1) return 0;

    long long inv3 = mod_inv(3);
    long long pow2 = mod_pow(2, n / 2);
    long long pref = (n % 2 == 0) ? 35 : 50;
    long long counter = mul_mod(pref, mul_mod(pow2, inv3));
    counter = sub_mod(counter, norm_mod(6LL * n));
    counter = sub_mod(counter, mul_mod(35, inv3));

    int lim = (n - 1) / 2;
    vector<long long> nums(lim + 1, 0);
    nums[0] = 1;
    if (lim >= 1) nums[1] = 10;

    for (int i = 2; i <= lim; ++i) {
        long long a = mul_mod(norm_mod(20LL * i - 10), nums[i - 1]);
        long long b = mul_mod(norm_mod(64LL * i - 64), nums[i - 2]);
        long long val = sub_mod(a, b);
        nums[i] = mul_mod(val, mod_inv(i));
    }

    for (int i = 0; i <= lim; ++i) {
        counter = add_mod(counter, mul_mod(nums[i], g8(n - 2LL * i, inv3)));
    }
    return counter;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    struct Check {
        int n;
        long long expected;
    };
    const Check checks[] = {
        {2, 6},
        {5, 348},
        {20, 125825982708LL},
    };

    for (const auto &check : checks) {
        long long got = compute_G(check.n);
        long long exp = check.expected % MOD;
        if (got != exp) {
            cerr << "Validation failed for G(" << check.n << "): got " << got
                 << ", expected " << exp << ".\n";
            return 1;
        }
    }
    cerr << "Validation checkpoints passed.\n";

    cout << compute_G(9898) << "\n";
    return 0;
}
