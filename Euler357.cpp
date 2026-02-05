#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

namespace {

constexpr int kLimit = 100000000;

vector<bool> sieve_primes(int n) {
    vector<bool> is_prime(n + 1, true);
    if (n >= 0) is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (int i = 4; i <= n; i += 2) {
        is_prime[i] = false;
    }
    for (int p = 3; 1LL * p * p <= n; p += 2) {
        if (!is_prime[p]) continue;
        for (long long j = 1LL * p * p; j <= n; j += 2LL * p) {
            is_prime[j] = false;
        }
    }
    return is_prime;
}

vector<uint32_t> build_spf_odd(int limit) {
    vector<uint32_t> spf(limit / 2 + 1, 0);
    int root = static_cast<int>(sqrt(limit));
    for (int p = 3; p <= root; p += 2) {
        if (spf[p / 2] != 0) continue;
        long long step = 2LL * p;
        long long start = 1LL * p * p;
        for (long long j = start; j <= limit; j += step) {
            if (spf[j / 2] == 0) spf[j / 2] = static_cast<uint32_t>(p);
        }
    }
    return spf;
}

bool check_n(int n,
             const vector<bool>& is_prime,
             const vector<uint32_t>& spf_odd,
             vector<int>& factors,
             vector<uint64_t>& divisors) {
    if (n == 1) return true;
    if (n & 1) return false;
    if ((n & 3) == 0) return false;

    int m = n / 2;
    factors.clear();
    factors.push_back(2);
    while (m > 1) {
        int p = spf_odd[m / 2];
        if (p == 0) p = m;
        factors.push_back(p);
        m /= p;
        if (m % p == 0) return false;
    }

    divisors.clear();
    divisors.push_back(1);
    for (int p : factors) {
        size_t size = divisors.size();
        for (size_t i = 0; i < size; ++i) {
            divisors.push_back(divisors[i] * static_cast<uint64_t>(p));
        }
    }

    uint64_t nn = static_cast<uint64_t>(n);
    for (uint64_t d : divisors) {
        uint64_t q = nn / d;
        if (d > q) continue;
        if (!is_prime[d + q]) return false;
    }
    return true;
}

uint64_t solve(int limit, int threads) {
    if (limit < 1) return 0;

    int max_prime = limit + 1;
    vector<bool> is_prime = sieve_primes(max_prime);

    vector<int> primes;
    primes.reserve(max_prime / 10);
    if (max_prime >= 2 && is_prime[2]) primes.push_back(2);
    for (int i = 3; i <= max_prime; i += 2) {
        if (is_prime[i]) primes.push_back(i);
    }

    int half = limit / 2;
    vector<uint32_t> spf_odd = build_spf_odd(half);

    if (threads < 1) threads = 1;
    if (threads > static_cast<int>(primes.size())) {
        threads = static_cast<int>(primes.size());
    }

    vector<uint64_t> partial(threads, 0);
    vector<thread> pool;
    pool.reserve(threads);
    size_t chunk = (primes.size() + threads - 1) / threads;

    for (int t = 0; t < threads; ++t) {
        size_t start = static_cast<size_t>(t) * chunk;
        size_t end = min(primes.size(), start + chunk);
        pool.emplace_back([&, start, end, t]() {
            vector<int> factors;
            factors.reserve(8);
            vector<uint64_t> divisors;
            divisors.reserve(256);
            uint64_t local = 0;
            for (size_t i = start; i < end; ++i) {
                int n = primes[i] - 1;
                if (n > limit) continue;
                if (check_n(n, is_prime, spf_odd, factors, divisors)) {
                    local += static_cast<uint64_t>(n);
                }
            }
            partial[t] = local;
        });
    }

    for (auto& th : pool) th.join();

    uint64_t total = 0;
    for (uint64_t s : partial) total += s;
    return total;
}

}  // namespace

int main(int argc, char** argv) {
    int limit = kLimit;
    int threads = static_cast<int>(thread::hardware_concurrency());
    if (threads <= 0) threads = 1;
    if (argc > 1) limit = max(1, atoi(argv[1]));
    if (argc > 2) threads = max(1, atoi(argv[2]));

    uint64_t answer = solve(limit, threads);
    cout << answer << endl;
    return 0;
}
