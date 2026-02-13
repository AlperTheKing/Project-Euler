#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u128 = unsigned __int128;

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }
    std::string s;
    while (x > 0) {
        s.push_back(static_cast<char>('0' + static_cast<int>(x % 10)));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

u128 multiset_weighted_count(const std::vector<u128>& counts_by_size, int child_count, int size_sum) {
    std::vector<std::vector<u128>> dp(child_count + 1, std::vector<u128>(size_sum + 1, 0));
    dp[0][0] = 1;

    for (int s = 1; s <= size_sum; ++s) {
        const u128 c = counts_by_size[s];
        if (c == 0) {
            continue;
        }

        const int q_max = std::min(child_count, size_sum / s);
        std::vector<u128> choose(q_max + 1, 0);
        choose[0] = 1;
        for (int q = 1; q <= q_max; ++q) {
            choose[q] = choose[q - 1] * (c + static_cast<u128>(q - 1)) / static_cast<u128>(q);
        }

        std::vector<std::vector<u128>> next = dp;
        for (int used = 0; used <= child_count; ++used) {
            for (int w = 0; w <= size_sum; ++w) {
                const u128 base = dp[used][w];
                if (base == 0) {
                    continue;
                }
                for (int q = 1; q <= q_max; ++q) {
                    const int n_used = used + q;
                    const int n_w = w + q * s;
                    if (n_used > child_count || n_w > size_sum) {
                        break;
                    }
                    next[n_used][n_w] += base * choose[q];
                }
            }
        }
        dp.swap(next);
    }

    return dp[child_count][size_sum];
}

std::vector<u128> compute_peerless_unrooted(int N) {
    std::vector<std::vector<u128>> planted(N + 1, std::vector<u128>(N + 1, 0));
    std::vector<u128> planted_total_by_size(N + 1, 0);

    for (int n = 1; n <= N; ++n) {
        for (int d = 1; d <= n; ++d) {
            const int child_count = d - 1;
            const int child_size_sum = n - 1;
            if (child_count > child_size_sum) {
                continue;
            }

            std::vector<u128> allowed_by_size(n, 0);
            for (int s = 1; s < n; ++s) {
                allowed_by_size[s] = planted_total_by_size[s] - planted[d][s];
            }
            planted[d][n] = multiset_weighted_count(allowed_by_size, child_count, child_size_sum);
        }

        u128 total = 0;
        for (int d = 1; d <= n; ++d) {
            total += planted[d][n];
        }
        planted_total_by_size[n] = total;
    }

    std::vector<u128> vertex_rooted(N + 1, 0);
    for (int n = 1; n <= N; ++n) {
        u128 total = 0;
        for (int d = 1; d <= n; ++d) {
            const int child_count = d;
            const int child_size_sum = n - 1;
            if (child_count > child_size_sum) {
                continue;
            }

            std::vector<u128> allowed_by_size(n, 0);
            for (int s = 1; s < n; ++s) {
                allowed_by_size[s] = planted_total_by_size[s] - planted[d][s];
            }
            total += multiset_weighted_count(allowed_by_size, child_count, child_size_sum);
        }
        vertex_rooted[n] = total;
    }

    std::vector<u128> peerless_unrooted(N + 1, 0);
    for (int n = 1; n <= N; ++n) {
        u128 directed_edge_rooted = 0;

        for (int left = 1; left < n; ++left) {
            const int right = n - left;
            directed_edge_rooted += planted_total_by_size[left] * planted_total_by_size[right];
        }

        u128 same_degree = 0;
        for (int d = 1; d <= N; ++d) {
            for (int left = 1; left < n; ++left) {
                const int right = n - left;
                same_degree += planted[d][left] * planted[d][right];
            }
        }
        directed_edge_rooted -= same_degree;

        assert((directed_edge_rooted & 1U) == 0);
        peerless_unrooted[n] = vertex_rooted[n] - directed_edge_rooted / 2;
    }

    return peerless_unrooted;
}

u128 prefix_sum_from_3(const std::vector<u128>& values, int N) {
    u128 total = 0;
    for (int n = 3; n <= N; ++n) {
        total += values[n];
    }
    return total;
}

}  // namespace

int main() {
    const std::vector<u128> peerless = compute_peerless_unrooted(50);
    assert(peerless[7] == 6);
    assert(prefix_sum_from_3(peerless, 10) == 74);
    std::cout << to_string_u128(prefix_sum_from_3(peerless, 50)) << '\n';
    return 0;
}
