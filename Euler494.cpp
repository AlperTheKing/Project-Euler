#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

using boost::multiprecision::cpp_int;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

std::string to_string_u128(u128 v) {
    if (v == 0) return "0";
    std::string s;
    while (v > 0) {
        s.push_back(static_cast<char>('0' + static_cast<int>(v % 10)));
        v /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

u128 fib(const int n) {
    if (n <= 0) return 0;
    if (n <= 2) return 1;
    u128 a = 1;
    u128 b = 1;
    for (int i = 2; i < n; ++i) {
        u128 c = a + b;
        a = b;
        b = c;
    }
    return b;
}

int coll_length(u128 n) {
    int len = 0;
    while ((n & (n - 1)) != 0) {
        ++len;
        if ((n & 1) != 0) {
            n = 3 * n + 1;
        } else {
            n >>= 1;
        }
    }
    return len;
}

u128 count_back(const std::vector<std::pair<cpp_int, int>>& sources, int precompute_length = -1) {
    u128 ans = 0;
    std::vector<std::pair<int, cpp_int>> pending;
    pending.reserve(sources.size());
    for (const auto& [n, steps] : sources) {
        if (steps < 0) continue;
        const int adjusted = (n % 3 != 0) ? steps : 0;
        pending.push_back({adjusted, n});
    }
    if (pending.empty()) return 0;

    std::sort(pending.begin(), pending.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) return a.first > b.first;
        return a.second > b.second;
    });

    if (precompute_length < 0) {
        precompute_length = std::max(0, (pending.front().first - 2) / 3);
    }

    const int max_ones = (precompute_length + 1) / 2;
    std::vector<u64> pow3(max_ones + 3, 1);
    for (int i = 1; i < static_cast<int>(pow3.size()); ++i) {
        pow3[i] = pow3[i - 1] * 3ULL;
    }

    std::vector<std::vector<u64>> counts;
    counts.push_back({1, 1});  // length 0, mod 2.

    for (int length = 1; length <= precompute_length; ++length) {
        while (!pending.empty() && pending.back().first == length - 1) {
            const cpp_int n = pending.back().second;
            pending.pop_back();
            for (int ones = 0; ones < static_cast<int>(counts.size()); ++ones) {
                const u64 mod = 2ULL * pow3[ones];
                const u64 idx = static_cast<u64>(n % mod);
                ans += counts[ones][idx];
            }
        }

        if ((length & 1) != 0) {
            const int ones = (length + 1) / 2;
            std::vector<u64> layer(2ULL * pow3[ones], 0);
            layer[(pow3[ones] - 1ULL) << 1] = 1;
            counts.push_back(std::move(layer));
        }

        for (int ones = length / 2; ones >= 1; --ones) {
            const u64 p3 = pow3[ones];
            const u64 size = 2ULL * p3;
            std::vector<u64> next(size, 0);

            const auto& cur = counts[ones];
            for (u64 i = 0; i < size; i += 2) {
                const u64 c = cur[i];
                if (c == 0) continue;
                const u64 j = i >> 1;
                next[j] += c;
                next[j + p3] += c;
            }

            const auto& prev = counts[ones - 1];
            const u64 prev_size = 2ULL * pow3[ones - 1];
            for (u64 i = 1; i < prev_size; i += 2) {
                const u64 c = prev[i];
                if (c == 0) continue;
                next[3ULL * i + 1ULL] += c;
            }
            counts[ones] = std::move(next);
        }
    }

    static constexpr int count4[18] = {1, 2, 3, 1, 4, 3, 1, 2, 4, 1, 3, 3, 1, 2, 3, 1, 4, 4};
    const int mod4_target = precompute_length % 4;

    for (auto [steps, n] : pending) {
        std::vector<cpp_int> states;
        states.push_back(n);

        while (steps % 4 != mod4_target) {
            std::vector<cpp_int> next;
            next.reserve(states.size() * 2);
            for (const cpp_int& a : states) {
                next.push_back(a << 1);
                if (a % 6 == 4) next.push_back(a / 3);
            }
            states.swap(next);
            --steps;
        }

        while (steps != precompute_length) {
            std::vector<cpp_int> next;
            next.reserve(states.size() * 4);
            for (const cpp_int& a : states) {
                if (a % 3 == 0) {
                    next.push_back(0);
                    continue;
                }
                const int c = count4[static_cast<int>(a % 18)];
                if (c >= 1) next.push_back(a << 4);
                if (c >= 2) {
                    const cpp_int t = (a << 3) / 3;
                    next.push_back(t);
                    if (c >= 3) next.push_back(t + (a % 3) - 3);
                    if (c >= 4) next.push_back(((a - 1) << 2) / 9);
                }
            }
            states.swap(next);
            steps -= 4;
        }

        for (const cpp_int& a : states) {
            for (int ones = 0; ones < static_cast<int>(counts.size()); ++ones) {
                const u64 mod = 2ULL * pow3[ones];
                const u64 idx = static_cast<u64>(a % mod);
                ans += counts[ones][idx];
            }
        }
    }
    return ans;
}

u128 solve_f(const int steps) {
    static constexpr int seeds[6] = {9, 19, 37, 51, 155, 159};
    std::vector<std::pair<cpp_int, int>> sources;
    sources.reserve(6);
    for (int s : seeds) {
        sources.push_back({cpp_int(s), steps - coll_length(static_cast<u128>(s))});
    }
    return fib(steps) + count_back(sources);
}

bool run_checkpoints() {
    if (solve_f(5) != static_cast<u128>(5)) return false;
    if (solve_f(10) != static_cast<u128>(55)) return false;
    if (solve_f(20) != static_cast<u128>(6771)) return false;
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!skip_checkpoints && !run_checkpoints()) {
        std::cerr << "Checkpoint failed\n";
        return 2;
    }

    std::cout << to_string_u128(solve_f(90)) << '\n';
    return 0;
}
