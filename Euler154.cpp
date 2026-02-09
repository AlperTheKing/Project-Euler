#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int n = 200000;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    if (options.threads <= 0) {
        options.threads = 1;
    }
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--n=", options.n) ||
            parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1 && options.threads >= 1;
}

void build_factorial_prime_exponents(const int n, std::vector<int>& v2, std::vector<int>& v5) {
    v2.assign(static_cast<std::size_t>(n + 1), 0);
    v5.assign(static_cast<std::size_t>(n + 1), 0);
    for (int i = 1; i <= n; ++i) {
        int x = i;
        int add2 = 0;
        while ((x & 1) == 0) {
            x >>= 1;
            ++add2;
        }

        x = i;
        int add5 = 0;
        while (x % 5 == 0) {
            x /= 5;
            ++add5;
        }

        v2[static_cast<std::size_t>(i)] = v2[static_cast<std::size_t>(i - 1)] + add2;
        v5[static_cast<std::size_t>(i)] = v5[static_cast<std::size_t>(i - 1)] + add5;
    }
}

u64 count_coefficients_fast(const int n, int threads) {
    std::vector<int> v2;
    std::vector<int> v5;
    build_factorial_prime_exponents(n, v2, v5);

    const int need2 = v2[static_cast<std::size_t>(n)] - 12;
    const int need5 = v5[static_cast<std::size_t>(n)] - 12;

    const int a_max = n / 3;
    threads = std::max(1, std::min(threads, a_max + 1));

    std::atomic<int> next_a{0};
    constexpr int kChunk = 64;

    std::vector<u64> partial(static_cast<std::size_t>(threads), 0);
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(threads));

    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            u64 local = 0;
            while (true) {
                const int start = next_a.fetch_add(kChunk, std::memory_order_relaxed);
                if (start > a_max) {
                    break;
                }
                const int end = std::min(a_max, start + kChunk - 1);
                for (int a = start; a <= end; ++a) {
                    const int a2 = v2[static_cast<std::size_t>(a)];
                    const int a5 = v5[static_cast<std::size_t>(a)];
                    const int b_max = (n - a) / 2;
                    for (int b = a; b <= b_max; ++b) {
                        const int c = n - a - b;
                        if (a2 + v2[static_cast<std::size_t>(b)] + v2[static_cast<std::size_t>(c)] > need2) {
                            continue;
                        }
                        if (a5 + v5[static_cast<std::size_t>(b)] + v5[static_cast<std::size_t>(c)] > need5) {
                            continue;
                        }

                        if (a == b && b == c) {
                            local += 1;
                        } else if (a == b || b == c) {
                            local += 3;
                        } else {
                            local += 6;
                        }
                    }
                }
            }
            partial[static_cast<std::size_t>(t)] = local;
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    u64 answer = 0;
    for (u64 v : partial) {
        answer += v;
    }
    return answer;
}

u64 count_coefficients_bruteforce(const int n) {
    std::vector<int> v2;
    std::vector<int> v5;
    build_factorial_prime_exponents(n, v2, v5);

    const int need2 = v2[static_cast<std::size_t>(n)] - 12;
    const int need5 = v5[static_cast<std::size_t>(n)] - 12;

    u64 answer = 0;
    for (int a = 0; a <= n; ++a) {
        const int a2 = v2[static_cast<std::size_t>(a)];
        const int a5 = v5[static_cast<std::size_t>(a)];
        for (int b = 0; b <= n - a; ++b) {
            const int c = n - a - b;
            if (a2 + v2[static_cast<std::size_t>(b)] + v2[static_cast<std::size_t>(c)] > need2) {
                continue;
            }
            if (a5 + v5[static_cast<std::size_t>(b)] + v5[static_cast<std::size_t>(c)] > need5) {
                continue;
            }
            ++answer;
        }
    }
    return answer;
}

bool run_checkpoints() {
    for (int n : {40, 60, 80, 120}) {
        const u64 slow = count_coefficients_bruteforce(n);
        const u64 fast = count_coefficients_fast(n, 1);
        if (slow != fast) {
            std::cerr << "Checkpoint failed at n=" << n << ": brute=" << slow
                      << ", fast=" << fast << '\n';
            return false;
        }
    }

    const int thread_check_n = 3000;
    const u64 single = count_coefficients_fast(thread_check_n, 1);
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) {
        hw = 2;
    }
    const int multi_threads = static_cast<int>(std::min<unsigned>(hw, 8));
    const u64 multi = count_coefficients_fast(thread_check_n, multi_threads);
    if (single != multi) {
        std::cerr << "Thread consistency failed at n=" << thread_check_n
                  << ": single=" << single << ", multi=" << multi << '\n';
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }
    std::cout << count_coefficients_fast(options.n, options.threads) << '\n';
    return 0;
}
