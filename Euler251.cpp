#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>
#include <tuple>
#include <string>
#include <sstream>

typedef __uint128_t ull128;
typedef unsigned long long ull;

std::mutex mtx;
ull total_triplet_count = 0; // Global count of total Cardano Triplets

// Custom function to handle input for __uint128_t
std::istream& operator>>(std::istream& in, ull128& value) {
    std::string str;
    in >> str;
    value = 0;
    for (char c : str) {
        if (c >= '0' && c <= '9') {
            value = value * 10 + (c - '0');
        } else {
            in.setstate(std::ios::failbit);
            break;
        }
    }
    return in;
}

// Output operator for 128-bit integers (__uint128_t)
std::ostream& operator<<(std::ostream& os, const ull128& value) {
    std::string str;
    ull128 tmp = value;
    do {
        str += '0' + tmp % 10;
        tmp /= 10;
    } while (tmp > 0);
    std::reverse(str.begin(), str.end());
    return os << str;
}

// Integer exponentiation for 128-bit numbers
ull128 int_pow(ull128 base, ull128 exp) {
    ull128 result = 1;
    while (exp > 0) {
        if (exp & 1)
            result *= base;
        base *= base;
        exp >>= 1;
    }
    return result;
}

// Function to factorize n using trial division
void factorize(ull128 n, std::vector<std::pair<ull128, ull128>>& factors) {
    for (ull128 i = 2; i * i <= n; ++i) {
        ull128 count = 0;
        while (n % i == 0) {
            n /= i;
            ++count;
        }
        if (count > 0) {
            factors.emplace_back(i, count);
        }
    }
    if (n > 1) {
        factors.emplace_back(n, 1); // n is prime
    }
}

// Recursive function to generate all possible (b, c) pairs and count them
void generate_bc(const std::vector<std::pair<ull128, ull128>>& factors, size_t idx, ull128 b, ull128 c, ull128 max_sum, ull128 a) {
    if (idx == factors.size()) {
        if (a + b + c <= max_sum && b > 0 && c > 0) {
            if (b > c) std::swap(b, c); // Ensure b <= c
            std::lock_guard<std::mutex> lock(mtx);
            total_triplet_count++;
        }
        return;
    }

    ull128 p = factors[idx].first;
    ull128 e = factors[idx].second;

    ull128 max_k = e / 2;
    for (ull128 k = 0; k <= max_k; ++k) {
        ull128 b_new = b * int_pow(p, k);
        ull128 c_new = c * int_pow(p, e - 2 * k);
        generate_bc(factors, idx + 1, b_new, c_new, max_sum, a);
    }
}

// Function to find Cardano Triplets in a given range of 'a'
void find_cardano_triplets(ull128 start_a, ull128 end_a, ull128 max_sum) {
    for (ull128 a = start_a; a <= end_a; a += 3) { // a ≡ 2 mod 3
        ull128 N = (1 + a) * (1 + a) * (8 * a - 1);
        if (N % 27 != 0)
            continue;
        ull128 N_div = N / 27;

        // Factorize N_div
        std::vector<std::pair<ull128, ull128>> factors;
        factorize(N_div, factors);

        // Generate all possible (b, c) pairs and count them
        generate_bc(factors, 0, 1, 1, max_sum, a);
    }
}

int main() {
    ull128 max_sum;
    std::cout << "Enter the maximum value for (a + b + c): ";
    std::cin >> max_sum;

    auto start_time = std::chrono::high_resolution_clock::now();

    // Determine the number of hardware threads available
    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0)
        num_threads = 4; // Default to 4 if unable to determine

    std::cout << "Number of threads: " << num_threads << std::endl;

    // Estimate MAX_A based on max_sum
    ull128 MAX_A = max_sum; // Conservative estimate

    ull128 range = MAX_A / num_threads + 1; // Ensure full coverage
    std::vector<std::thread> threads;

    for (unsigned int i = 0; i < num_threads; ++i) {
        ull128 start_a = 2 + i * range;
        if (start_a % 3 != 2)
            start_a += (3 - (start_a % 3) + 2) % 3; // Adjust to a ≡ 2 mod 3
        ull128 end_a = std::min(start_a + range - 1, MAX_A);

        if (start_a > end_a)
            continue;

        threads.emplace_back(find_cardano_triplets, start_a, end_a, max_sum);
    }

    for (auto& t : threads) {
        t.join();
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;

    std::cout << "Total Cardano Triplets: " << total_triplet_count << std::endl;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";

    return 0;
}