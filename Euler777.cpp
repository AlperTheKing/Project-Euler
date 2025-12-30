#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace std;

namespace {

using i128 = __int128_t;

inline long long sum_upto(long long n) {
    return static_cast<long long>((static_cast<i128>(n) * (n + 1)) / 2);
}

vector<int> mobius_sieve(int n) {
    vector<int> mu(n + 1, 1);
    vector<int> primes;
    vector<bool> composite(n + 1, false);
    mu[0] = 0;
    for (int i = 2; i <= n; ++i) {
        if (!composite[i]) {
            primes.push_back(i);
            mu[i] = -1;
        }
        for (int p : primes) {
            long long v = static_cast<long long>(i) * p;
            if (v > n) break;
            composite[static_cast<int>(v)] = true;
            if (i % p == 0) {
                mu[static_cast<int>(v)] = 0;
                break;
            }
            mu[static_cast<int>(v)] = -mu[i];
        }
    }
    return mu;
}

enum class SubsetType {
    S10,
    S5Odd,
    S2Not5,
    SCoprime10,
};

inline void count_sum_S10(int n, int d, long long& cnt, long long& sum) {
    int g = 1;
    if (d % 10 == 0) g = 10;
    else if (d % 5 == 0) g = 5;
    else if (d % 2 == 0) g = 2;
    long long L = static_cast<long long>(d / g) * 10LL;
    long long t = n / L;
    cnt = t;
    sum = static_cast<long long>((static_cast<i128>(L) * t * (t + 1)) / 2);
}

inline void count_sum_S5Odd(int n, int d, long long& cnt, long long& sum) {
    if (d % 2 == 0) {
        cnt = 0;
        sum = 0;
        return;
    }
    int g = (d % 5 == 0) ? 5 : 1;
    long long L = static_cast<long long>(d / g) * 5LL;
    long long t = n / L;
    long long k = (t + 1) / 2;
    cnt = k;
    sum = L * k * k;
}

inline void count_sum_S2Not5(int n, int d, long long& cnt, long long& sum) {
    if (d % 5 == 0) {
        cnt = 0;
        sum = 0;
        return;
    }
    long long L = (d % 2 == 0) ? d : 2LL * d;
    long long t = n / L;
    long long t5 = n / (5LL * L);
    cnt = t - t5;
    long long sumL = static_cast<long long>((static_cast<i128>(L) * t * (t + 1)) / 2);
    long long sum5L = static_cast<long long>((static_cast<i128>(5LL * L) * t5 * (t5 + 1)) / 2);
    sum = sumL - sum5L;
}

inline void count_sum_SCoprime10(int n, int d, long long& cnt, long long& sum) {
    if (d % 2 == 0 || d % 5 == 0) {
        cnt = 0;
        sum = 0;
        return;
    }
    long long x = n / d;
    long long t2 = x / 2;
    long long t5 = x / 5;
    long long t10 = x / 10;
    cnt = x - t2 - t5 + t10;
    long long sum_all = sum_upto(x);
    long long sum2 = static_cast<long long>(2LL * sum_upto(t2));
    long long sum5 = static_cast<long long>(5LL * sum_upto(t5));
    long long sum10 = static_cast<long long>(10LL * sum_upto(t10));
    long long sum_t = sum_all - sum2 - sum5 + sum10;
    sum = static_cast<long long>(static_cast<i128>(d) * sum_t);
}

void compute_subset(int n, const vector<int>& mu, SubsetType type,
                    vector<long long>& cnt, vector<long long>& sum) {
    for (int d = 1; d <= n; ++d) {
        int mud = mu[d];
        if (mud == 0) continue;
        long long c = 0;
        long long s = 0;
        switch (type) {
            case SubsetType::S10:
                count_sum_S10(n, d, c, s);
                break;
            case SubsetType::S5Odd:
                count_sum_S5Odd(n, d, c, s);
                break;
            case SubsetType::S2Not5:
                count_sum_S2Not5(n, d, c, s);
                break;
            case SubsetType::SCoprime10:
                count_sum_SCoprime10(n, d, c, s);
                break;
        }
        if (c == 0) continue;
        for (int a = d; a <= n; a += d) {
            cnt[a] += static_cast<long long>(mud) * c;
            sum[a] += static_cast<long long>(mud) * s;
        }
    }
}

inline int gcd10_bucket(int a) {
    if (a % 10 == 0) return 10;
    if (a % 5 == 0) return 5;
    if (a % 2 == 0) return 2;
    return 1;
}

i128 compute_s_num(int n) {
    vector<int> mu = mobius_sieve(n);

    i128 sum_ab = 0;
    i128 sum_a = 0;
    for (int d = 1; d <= n; ++d) {
        int mud = mu[d];
        if (mud == 0) continue;
        long long k = n / d;
        long long s = sum_upto(k);
        i128 dd = d;
        i128 ss = s;
        sum_ab += static_cast<i128>(mud) * dd * dd * ss * ss;
        sum_a += static_cast<i128>(mud) * dd * ss * k;
    }

    long long s_all = sum_upto(n);
    i128 sum_ab_2n = sum_ab - 2 * static_cast<i128>(s_all) + 1;
    i128 sum_a_2n = sum_a - n - static_cast<i128>(s_all) + 1;

    vector<long long> cnt10(n + 1, 0), sum10(n + 1, 0);
    vector<long long> cnt5o(n + 1, 0), sum5o(n + 1, 0);
    vector<long long> cnt2n(n + 1, 0), sum2n(n + 1, 0);
    vector<long long> cntg1(n + 1, 0), sumg1(n + 1, 0);

    unsigned hw = thread::hardware_concurrency();
    bool use_threads = hw >= 2 && n >= 2000;
    if (use_threads) {
        thread t1(compute_subset, n, cref(mu), SubsetType::S10, ref(cnt10), ref(sum10));
        thread t2(compute_subset, n, cref(mu), SubsetType::S5Odd, ref(cnt5o), ref(sum5o));
        thread t3(compute_subset, n, cref(mu), SubsetType::S2Not5, ref(cnt2n), ref(sum2n));
        compute_subset(n, mu, SubsetType::SCoprime10, cntg1, sumg1);
        t1.join();
        t2.join();
        t3.join();
    } else {
        compute_subset(n, mu, SubsetType::S10, cnt10, sum10);
        compute_subset(n, mu, SubsetType::S5Odd, cnt5o, sum5o);
        compute_subset(n, mu, SubsetType::S2Not5, cnt2n, sum2n);
        compute_subset(n, mu, SubsetType::SCoprime10, cntg1, sumg1);
    }

    i128 count10 = 0;
    i128 sum_a10 = 0;
    i128 sum_ab10 = 0;
    for (int a = 2; a <= n; ++a) {
        long long cnt = 0;
        long long s = 0;
        int g = gcd10_bucket(a);
        if (g == 1) {
            cnt = cnt10[a];
            s = sum10[a];
        } else if (g == 2) {
            cnt = cnt5o[a];
            s = sum5o[a];
        } else if (g == 5) {
            cnt = cnt2n[a];
            s = sum2n[a];
        } else {
            cnt = cntg1[a];
            s = sumg1[a];
            cnt -= 1;  // exclude b=1
            s -= 1;
        }
        count10 += cnt;
        sum_a10 += static_cast<i128>(a) * cnt;
        sum_ab10 += static_cast<i128>(a) * s;
    }

    i128 base = 2 * sum_ab_2n - 3 * sum_a_2n;
    i128 corr = 6 * sum_a10 + 4 * count10 - 6 * sum_ab10;
    i128 total_num = 4 * base + corr;
    return total_num;
}

string format_scientific(long double value) {
    ostringstream oss;
    oss.setf(ios::scientific);
    oss << setprecision(9) << value;
    string s = oss.str();
    size_t pos = s.find('e');
    if (pos == string::npos) return s;
    char sign = s[pos + 1];
    size_t i = pos + 2;
    while (i < s.size() && s[i] == '0') ++i;
    string exp = (i < s.size()) ? s.substr(i) : "0";
    string out = s.substr(0, pos + 1);
    if (sign == '-') out.push_back('-');
    out += exp;
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n = 1'000'000;
    if (argc >= 2) n = stoi(argv[1]);

    const i128 expect10 = 6410;      // 1602.5 * 4
    const i128 expect100 = 97026020; // 24256505 * 4

    i128 val10 = compute_s_num(10);
    if (val10 != expect10) {
        cerr << "Validation failed: s(10) expected 1602.5, got "
             << static_cast<long double>(val10) / 4.0L << "\n";
        return 1;
    }
    i128 val100 = compute_s_num(100);
    if (val100 != expect100) {
        cerr << "Validation failed: s(100) expected 24256505, got "
             << static_cast<long double>(val100) / 4.0L << "\n";
        return 1;
    }

    i128 total_num = compute_s_num(n);
    long double value = static_cast<long double>(total_num) / 4.0L;
    cout << format_scientific(value) << "\n";
    return 0;
}
