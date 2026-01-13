#include <iostream>
#include <vector>
#include <cmath>
#include <future>
#include <thread>
#include <numeric>

// Type definitions for handling large numbers
using int128 = __int128;
using ll = long long;

const ll MOD = 1000000007;

// Modular arithmetic helper to handle negative results
ll safe_mod(int128 val) {
    val %= MOD;
    if (val < 0) val += MOD;
    return (ll)val;
}

// Modular inverse of 2 for division
const ll INV2 = 500000004; // (MOD + 1) / 2

// Class to handle the Sieve of Eratosthenes for Mobius function
class Sieve {
public:
    std::vector<int> mu;
    
    Sieve(int n) {
        mu.resize(n + 1);
        std::vector<bool> is_prime(n + 1, true);
        std::vector<int> primes;
        
        mu[1] = 1;
        for (int i = 2; i <= n; ++i) {
            if (is_prime[i]) {
                primes.push_back(i);
                mu[i] = -1;
            }
            for (int p : primes) {
                if (i * p > n) break;
                is_prime[i * p] = false;
                if (i % p == 0) {
                    mu[i * p] = 0;
                    break;
                }
                mu[i * p] = -mu[i];
            }
        }
    }
};

// Calculate sum of sigma_1(k) for k=1 to x
// T(x) = sum_{i=1}^x i * floor(x/i)
ll calc_T(ll x) {
    ll total_sum = 0;
    for (ll l = 1, r; l <= x; l = r + 1) {
        r = x / (x / l);
        // Sum of arithmetic progression from l to r: (l+r)*(r-l+1)/2
        // We do calculation in int128 to prevent overflow before modulo
        int128 count = (r - l + 1);
        int128 sum_l_r = (l + r);
        
        // (count * sum_l_r) / 2 % MOD
        ll sum_seq = safe_mod((count * sum_l_r) % (2 * MOD) / 2);
        
        int128 term = (int128)sum_seq * (x / l);
        total_sum = safe_mod(total_sum + term);
    }
    return total_sum;
}

// Solver function
ll solve(ll N, const Sieve& sieve) {
    ll limit = sqrt(N);
    
    // We need to compute sum_{d=1}^{limit} d * mu[d] * T(N/d^2)
    // We will parallelize this loop.
    
    int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 4;
    
    std::vector<std::future<ll>> futures;
    
    auto task = [&](int thread_id) -> ll {
        ll local_sum = 0;
        // Cyclic distribution for load balancing (small d is expensive, large d is cheap)
        for (ll d = thread_id + 1; d <= limit; d += num_threads) {
            if (sieve.mu[d] == 0) continue;
            
            ll inner = calc_T(N / (d * d));
            int128 term = (int128)d * inner;
            term = safe_mod(term);
            
            if (sieve.mu[d] == 1) {
                local_sum = safe_mod(local_sum + term);
            } else {
                local_sum = safe_mod(local_sum - term);
            }
        }
        return local_sum;
    };
    
    for (int i = 0; i < num_threads; ++i) {
        futures.push_back(std::async(std::launch::async, task, i));
    }
    
    ll S_N = 0;
    for (auto& f : futures) {
        S_N = safe_mod(S_N + f.get());
    }
    
    // F(N) = S(N) - N(N+1)/2
    int128 n_sum = (int128)N * (N + 1);
    // Be careful with modulo division
    ll n_sum_mod = safe_mod(n_sum % (2 * MOD) / 2);
    
    return safe_mod(S_N - n_sum_mod);
}

int main() {
    // 1. Validation Check Point
    std::cout << "Running validation for N = 10^7..." << std::endl;
    
    // Precompute sieve up to sqrt(10^14) = 10^7
    // This covers the requirement for both 10^7 and 10^14 cases.
    Sieve sieve(10000000); 
    
    ll val_N = 10000000;
    ll val_result = solve(val_N, sieve);
    ll expected = 638042271;
    
    std::cout << "F(10^7) = " << val_result << std::endl;
    
    if (val_result == expected) {
        std::cout << "Validation PASSED. Proceeding to solve for 10^14..." << std::endl;
    } else {
        std::cout << "Validation FAILED. Expected " << expected << ". Aborting." << std::endl;
        return 1;
    }
    
    std::cout << "------------------------------------------------" << std::endl;
    
    // 2. Main Problem
    ll target_N = 100000000000000LL; // 10^14
    ll result = solve(target_N, sieve);
    
    std::cout << "F(10^14) modulo 10^9 + 7 is:" << std::endl;
    std::cout << result << std::endl;
    
    return 0;
}