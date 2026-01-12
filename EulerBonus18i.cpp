#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <future>
#include <thread>
#include <mutex>
#include <cstring>

const long long RANGE_START = 1000000000;
const long long RANGE_END   = 1100000000;
const int NUM_THREADS = std::thread::hardware_concurrency();

long long power(long long base, long long exp, long long mod) {
    long long res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (__int128)res * base % mod;
        base = (__int128)base * base % mod;
        exp /= 2;
    }
    return res;
}

struct Poly {
    long long a, b, c;

    static Poly multiply(const Poly& P, const Poly& Q, long long p) {
        __int128 P_a = P.a, P_b = P.b, P_c = P.c;
        __int128 Q_a = Q.a, Q_b = Q.b, Q_c = Q.c;

        __int128 c4 = P_a * Q_a;
        __int128 c3 = P_a * Q_b + P_b * Q_a;
        __int128 c2 = P_a * Q_c + P_b * Q_b + P_c * Q_a;
        __int128 c1 = P_b * Q_c + P_c * Q_b;
        __int128 c0 = P_c * Q_c;

        if (c4 != 0) {
            c4 %= p;
            c2 = (c2 + c4 * 3) % p;
            c1 = (c1 - c4 * 4) % p;
        }

        if (c3 != 0) {
            c3 %= p;
            c1 = (c1 + c3 * 3) % p;
            c0 = (c0 - c3 * 4) % p;
        }

        Poly res;
        res.a = (long long)((c2 % p + p) % p);
        res.b = (long long)((c1 % p + p) % p);
        res.c = (long long)((c0 % p + p) % p);
        return res;
    }

    static Poly pow(Poly base, long long exp, long long p) {
        Poly res = {0, 0, 1};
        while (exp > 0) {
            if (exp % 2 == 1) res = multiply(res, base, p);
            base = multiply(base, base, p);
            exp /= 2;
        }
        return res;
    }
};

long long resultant_mod(long long c2, long long c1, long long c0, long long p) {
    long long mat[5][5] = {
        {1, 0, -3, 4, 0},
        {0, 1, 0, -3, 4},
        {c2, c1, c0, 0, 0},
        {0, c2, c1, c0, 0},
        {0, 0, c2, c1, c0}
    };

    for(int i=2; i<5; ++i) 
        for(int j=0; j<5; ++j) 
            mat[i][j] = (mat[i][j] % p + p) % p;

    long long det = 1;

    for (int i = 0; i < 5; ++i) {
        int pivot = i;
        while (pivot < 5 && mat[pivot][i] == 0) pivot++;
        
        if (pivot == 5) return 0;
        
        if (pivot != i) {
            for (int j = 0; j < 5; ++j) std::swap(mat[i][j], mat[pivot][j]);
            det = -det;
        }
        
        det = (det * mat[i][i]) % p;
        long long inv = power(mat[i][i], p - 2, p);

        for (int k = i + 1; k < 5; ++k) {
            if (mat[k][i] != 0) {
                long long factor = (__int128)mat[k][i] * inv % p;
                for (int j = i; j < 5; ++j) {
                    long long sub = (__int128)factor * mat[i][j] % p;
                    mat[k][j] = (mat[k][j] - sub + p) % p;
                }
            }
        }
    }
    return (det + p) % p;
}

long long compute_Rp(long long p) {
    if (p % 4 == 3) return 0;

    Poly x = {0, 1, 0};
    Poly h = Poly::pow(x, p, p);

    long long c2 = (-h.a % p + p) % p;
    long long c1 = ((1 - h.b) % p + p) % p;
    long long c0 = (-h.c % p + p) % p;

    return resultant_mod(c2, c1, c0, p);
}

void validate() {
    std::cout << "Running validation checks...\n";
    
    long long r11 = compute_Rp(11);
    if (r11 != 0) std::cerr << "FAIL: R(11) should be 0, got " << r11 << "\n";
    else std::cout << "PASS: R(11) = 0\n";

    long long r29 = compute_Rp(29);
    if (r29 != 13) std::cerr << "FAIL: R(29) should be 13, got " << r29 << "\n";
    else std::cout << "PASS: R(29) = 13\n";

    long long r5 = compute_Rp(5);
    if (r5 != 1) std::cerr << "FAIL: R(5) should be 1, got " << r5 << "\n";
    else std::cout << "PASS: R(5) = 1\n";
    
    std::cout << "Validation complete.\n\n";
}

std::mutex total_sum_mutex;
unsigned __int128 total_sum = 0;

void process_chunk(long long start, long long end, const std::vector<int>& small_primes) {
    long long range_len = end - start;
    std::vector<bool> is_prime(range_len, true);

    for (int p : small_primes) {
        long long start_idx = (start + p - 1) / p;
        long long j = std::max(start_idx, 2LL) * p - start;
        for (; j < range_len; j += p) {
            is_prime[j] = false;
        }
    }

    unsigned __int128 local_sum = 0;
    for (int i = 0; i < range_len; ++i) {
        if (is_prime[i]) {
            long long p = start + i;
            if (p < 2) continue;
            local_sum += compute_Rp(p);
        }
    }

    std::lock_guard<std::mutex> lock(total_sum_mutex);
    total_sum += local_sum;
}

int main() {
    validate();

    std::cout << "Generating small primes for sieve...\n";
    long long limit = sqrt(RANGE_END) + 1;
    std::vector<int> small_primes;
    std::vector<bool> sieve(limit + 1, true);
    for (long long p = 2; p <= limit; ++p) {
        if (sieve[p]) {
            small_primes.push_back(p);
            for (long long i = p * p; i <= limit; i += p) sieve[i] = false;
        }
    }

    std::cout << "Processing range [" << RANGE_START << ", " << RANGE_END << "] with " << NUM_THREADS << " threads...\n";

    std::vector<std::thread> threads;
    long long chunk_size = (RANGE_END - RANGE_START) / NUM_THREADS;
    
    for (int i = 0; i < NUM_THREADS; ++i) {
        long long start = RANGE_START + i * chunk_size;
        long long end = (i == NUM_THREADS - 1) ? RANGE_END + 1 : start + chunk_size;
        threads.emplace_back(process_chunk, start, end, std::ref(small_primes));
    }

    for (auto& t : threads) t.join();

    long long upper = (long long)(total_sum / 1000000000000000000ULL);
    long long lower = (long long)(total_sum % 1000000000000000000ULL);
    
    std::cout << "Total Sum: ";
    if (upper > 0) std::cout << upper << lower;
    else std::cout << (long long)total_sum;
    std::cout << std::endl;

    return 0;
}