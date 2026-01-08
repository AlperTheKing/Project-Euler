#include <iostream>
#include <vector>
#include <numeric>
#include <map>
#include <cmath>
#include <algorithm>
#include <future>
#include <mutex>
#include <iomanip>

using namespace std;

long long MOD = 1000000007;

long long power(long long base, long long exp) {
    long long res = 1;
    base %= MOD;
    while (exp > 0) {
        if (exp % 2 == 1) res = (res * base) % MOD;
        base = (base * base) % MOD;
        exp /= 2;
    }
    return res;
}

long long modInverse(long long n) {
    return power(n, MOD - 2);
}

long long factorial[101];
long long invFactorial[101];

void precomputeFactorials() {
    factorial[0] = 1;
    invFactorial[0] = 1;
    for (int i = 1; i <= 100; i++) {
        factorial[i] = (factorial[i - 1] * i) % MOD;
        invFactorial[i] = modInverse(factorial[i]);
    }
}

// Struct to represent an integer partition
struct Partition {
    vector<int> parts;
};

// Generate all partitions of k
void generatePartitions(int target, int minVal, vector<int>& current, vector<Partition>& result) {
    if (target == 0) {
        result.push_back({current});
        return;
    }
    for (int i = minVal; i <= target; i++) {
        current.push_back(i);
        generatePartitions(target - i, i, current, result);
        current.pop_back();
    }
}

// Get exponents of prime factorization of n!
map<int, int> getFactorialExponents(int n) {
    map<int, int> exponents;
    // Sieve
    vector<bool> is_prime(n + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int p = 2; p <= n; p++) {
        if (is_prime[p]) {
            for (int i = 2 * p; i <= n; i += p)
                is_prime[i] = false;
            
            // Legendre's Formula
            int count = 0;
            long long p_pow = p;
            while (p_pow <= n) {
                count += n / p_pow;
                p_pow *= p;
            }
            exponents[p] = count;
        }
    }
    return exponents;
}

// Get exponents for a specific number (for validation)
map<int, int> getNumberExponents(int n) {
    map<int, int> exponents;
    for (int i = 2; i * i <= n; ++i) {
        while (n % i == 0) {
            exponents[i]++;
            n /= i;
        }
    }
    if (n > 1) exponents[n]++;
    return exponents;
}

// Solver function
long long solveW(int n_val, int k_val, bool isFactorial) {
    // 1. Get Exponents
    map<int, int> prime_counts; // map exponent value -> frequency of primes with that exponent
    int max_exponent = 0;

    map<int, int> raw_exponents;
    if (isFactorial) {
        raw_exponents = getFactorialExponents(n_val);
    } else {
        raw_exponents = getNumberExponents(n_val);
    }

    for (auto const& [prime, exp] : raw_exponents) {
        prime_counts[exp]++;
        if (exp > max_exponent) max_exponent = exp;
    }

    // 2. Generate Partitions of k
    vector<Partition> partitions;
    vector<int> current;
    generatePartitions(k_val, 1, current, partitions);

    // 3. Multithreaded Processing
    long long total_W = 0;
    mutex total_mutex;

    // Split partitions among threads
    int num_threads = thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 4;
    vector<future<void>> futures;
    
    int part_count = partitions.size();
    int chunk_size = (part_count + num_threads - 1) / num_threads;

    auto worker = [&](int start, int end) {
        long long local_sum = 0;
        // DP array to be reused
        // dp[e] = ways to write e as sum of parts from partition
        vector<long long> dp(max_exponent + 1);

        for (int i = start; i < end; i++) {
            const auto& p = partitions[i];
            
            // Calculate Coefficient for this partition
            // Coeff = (1 / Product(count(v)!)) * Product( (-1)^(part-1) / part )
            map<int, int> part_counts;
            long long term_coeff = 1;
            
            for (int val : p.parts) {
                part_counts[val]++;
                
                long long inv_val = modInverse(val);
                // (-1)^(val-1)
                long long sign = ((val - 1) % 2 == 0) ? 1 : -1;
                long long term = (sign * inv_val) % MOD;
                if (term < 0) term += MOD;
                
                term_coeff = (term_coeff * term) % MOD;
            }

            for (auto const& [val, count] : part_counts) {
                term_coeff = (term_coeff * invFactorial[count]) % MOD;
            }

            // Calculate S(lambda)
            // Build DP table for this partition
            // Effectively coefficient of x^E in Product(1/(1-x^val))
            fill(dp.begin(), dp.end(), 0);
            dp[0] = 1;

            for (int val : p.parts) {
                for (int j = val; j <= max_exponent; j++) {
                    dp[j] = (dp[j] + dp[j - val]); 
                    if (dp[j] >= MOD) dp[j] -= MOD;
                }
            }

            long long S_lambda = 1;
            for (auto const& [exp, freq] : prime_counts) {
                long long ways = dp[exp];
                if (ways == 0) {
                    S_lambda = 0;
                    break;
                }
                // (ways ^ freq)
                S_lambda = (S_lambda * power(ways, freq)) % MOD;
            }

            long long total_term = (term_coeff * S_lambda) % MOD;
            local_sum = (local_sum + total_term) % MOD;
        }

        lock_guard<mutex> lock(total_mutex);
        total_W = (total_W + local_sum) % MOD;
    };

    for (int t = 0; t < num_threads; t++) {
        int start = t * chunk_size;
        int end = min(start + chunk_size, part_count);
        if (start < end) {
            futures.push_back(async(launch::async, worker, start, end));
        }
    }

    for (auto& f : futures) {
        f.get();
    }

    return total_W;
}

int main() {
    precomputeFactorials();

    cout << "--- Validation Checkpoints ---" << endl;

    // Check 1: W(144, 4) should be 7
    long long res1 = solveW(144, 4, false);
    cout << "W(144, 4) = " << res1 << (res1 == 7 ? " [PASS]" : " [FAIL]") << endl;

    // Check 2: W(100!, 10) should be 287549200
    long long res2 = solveW(100, 10, true);
    cout << "W(100!, 10) = " << res2 << (res2 == 287549200 ? " [PASS]" : " [FAIL]") << endl;

    cout << "\n--- Final Solution ---" << endl;
    
    // Final: W(10000!, 30)
    cout << "Calculating W(10000!, 30)..." << endl;
    long long final_res = solveW(10000, 30, true);
    
    cout << "W(10000!, 30) modulo 10^9+7 = " << final_res << endl;

    return 0;
}