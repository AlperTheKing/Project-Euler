#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

std::vector<std::vector<int>> generate_partitions(const int n) {
    std::vector<std::vector<int>> parts;
    std::vector<int> cur;

    const auto dfs = [&](const auto& self, int rem, int max_part) -> void {
        if (rem == 0) {
            parts.push_back(cur);
            return;
        }
        for (int x = std::min(rem, max_part); x >= 1; --x) {
            cur.push_back(x);
            self(self, rem - x, x);
            cur.pop_back();
        }
    };

    dfs(dfs, n, n);
    return parts;
}

std::vector<int> build_representative_permutation(const int n, const std::vector<int>& partition) {
    std::vector<int> perm(static_cast<std::size_t>(n), 0);
    int cur = 0;
    for (const int len : partition) {
        if (len == 1) {
            perm[static_cast<std::size_t>(cur)] = cur;
            ++cur;
            continue;
        }
        for (int i = 0; i < len - 1; ++i) {
            perm[static_cast<std::size_t>(cur + i)] = cur + i + 1;
        }
        perm[static_cast<std::size_t>(cur + len - 1)] = cur;
        cur += len;
    }
    return perm;
}

std::vector<int> cycle_type_partition(const std::vector<int>& perm) {
    const int n = static_cast<int>(perm.size());
    std::vector<char> vis(static_cast<std::size_t>(n), 0);
    std::vector<int> lengths;
    lengths.reserve(static_cast<std::size_t>(n));

    for (int i = 0; i < n; ++i) {
        if (vis[static_cast<std::size_t>(i)] != 0) {
            continue;
        }
        int j = i;
        int len = 0;
        while (vis[static_cast<std::size_t>(j)] == 0) {
            vis[static_cast<std::size_t>(j)] = 1;
            j = perm[static_cast<std::size_t>(j)];
            ++len;
        }
        lengths.push_back(len);
    }

    std::sort(lengths.begin(), lengths.end(), std::greater<int>());
    return lengths;
}

u64 factorial(const int n) {
    u64 v = 1;
    for (int i = 2; i <= n; ++i) {
        v *= static_cast<u64>(i);
    }
    return v;
}

u64 class_size(const int n, const std::vector<int>& partition) {
    const u64 num = factorial(n);
    std::vector<int> freq(static_cast<std::size_t>(n + 1), 0);
    for (const int len : partition) {
        ++freq[static_cast<std::size_t>(len)];
    }

    u64 den = 1;
    for (int len = 1; len <= n; ++len) {
        const int cnt = freq[static_cast<std::size_t>(len)];
        if (cnt == 0) {
            continue;
        }
        u64 p = 1;
        for (int i = 0; i < cnt; ++i) {
            p *= static_cast<u64>(len);
        }
        den *= p;
        den *= factorial(cnt);
    }

    return num / den;
}

long double solve_average_expectation(const int n) {
    const std::vector<std::vector<int>> partitions = generate_partitions(n);
    const int ccount = static_cast<int>(partitions.size());

    std::map<std::vector<int>, int> class_index;
    for (int i = 0; i < ccount; ++i) {
        class_index[partitions[static_cast<std::size_t>(i)]] = i;
    }

    std::vector<std::vector<int>> reps;
    reps.reserve(static_cast<std::size_t>(ccount));
    for (const auto& part : partitions) {
        reps.push_back(build_representative_permutation(n, part));
    }

    const std::vector<int> identity_partition(static_cast<std::size_t>(n), 1);
    const int id_idx = class_index[identity_partition];

    std::vector<int> unknown_classes;
    unknown_classes.reserve(static_cast<std::size_t>(ccount - 1));
    for (int i = 0; i < ccount; ++i) {
        if (i != id_idx) {
            unknown_classes.push_back(i);
        }
    }
    const int m = static_cast<int>(unknown_classes.size());

    const int trans_total = n * (n - 1) / 2;
    const int cycle3_total = n * (n - 1) * (n - 2) / 3;  // oriented 3-cycles: 2*C(n,3)

    std::vector<std::vector<long double>> a(static_cast<std::size_t>(m),
                                            std::vector<long double>(static_cast<std::size_t>(m + 1), 0.0L));

    for (int ri = 0; ri < m; ++ri) {
        const int cls = unknown_classes[static_cast<std::size_t>(ri)];
        const std::vector<int>& perm = reps[static_cast<std::size_t>(cls)];

        std::vector<int> cnt2(static_cast<std::size_t>(ccount), 0);
        std::vector<int> cnt3(static_cast<std::size_t>(ccount), 0);

        // Multiply on the right by transpositions.
        for (int x = 0; x < n; ++x) {
            for (int y = x + 1; y < n; ++y) {
                std::vector<int> q = perm;
                std::swap(q[static_cast<std::size_t>(x)], q[static_cast<std::size_t>(y)]);
                const int to = class_index[cycle_type_partition(q)];
                ++cnt2[static_cast<std::size_t>(to)];
            }
        }

        // Multiply on the right by oriented 3-cycles.
        for (int x = 0; x < n; ++x) {
            for (int y = x + 1; y < n; ++y) {
                for (int z = y + 1; z < n; ++z) {
                    {
                        std::vector<int> q = perm;
                        const int px = q[static_cast<std::size_t>(x)];
                        const int py = q[static_cast<std::size_t>(y)];
                        const int pz = q[static_cast<std::size_t>(z)];
                        q[static_cast<std::size_t>(x)] = py;
                        q[static_cast<std::size_t>(y)] = pz;
                        q[static_cast<std::size_t>(z)] = px;
                        const int to = class_index[cycle_type_partition(q)];
                        ++cnt3[static_cast<std::size_t>(to)];
                    }
                    {
                        std::vector<int> q = perm;
                        const int px = q[static_cast<std::size_t>(x)];
                        const int py = q[static_cast<std::size_t>(y)];
                        const int pz = q[static_cast<std::size_t>(z)];
                        q[static_cast<std::size_t>(x)] = pz;
                        q[static_cast<std::size_t>(y)] = px;
                        q[static_cast<std::size_t>(z)] = py;
                        const int to = class_index[cycle_type_partition(q)];
                        ++cnt3[static_cast<std::size_t>(to)];
                    }
                }
            }
        }

        for (int cj = 0; cj < m; ++cj) {
            const int cls_to = unknown_classes[static_cast<std::size_t>(cj)];
            const long double p2 =
                static_cast<long double>(cnt2[static_cast<std::size_t>(cls_to)]) / static_cast<long double>(trans_total);
            const long double p3 =
                static_cast<long double>(cnt3[static_cast<std::size_t>(cls_to)]) / static_cast<long double>(cycle3_total);

            long double coef = -0.5L * p2 - (1.0L / 3.0L) * p3;
            if (cls == cls_to) {
                coef += 5.0L / 6.0L;
            }
            a[static_cast<std::size_t>(ri)][static_cast<std::size_t>(cj)] = coef;
        }

        a[static_cast<std::size_t>(ri)][static_cast<std::size_t>(m)] = 1.0L;
    }

    // Gaussian elimination with partial pivoting.
    for (int col = 0; col < m; ++col) {
        int pivot = col;
        for (int r = col + 1; r < m; ++r) {
            if (std::fabsl(a[static_cast<std::size_t>(r)][static_cast<std::size_t>(col)]) >
                std::fabsl(a[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(col)])) {
                pivot = r;
            }
        }
        std::swap(a[static_cast<std::size_t>(col)], a[static_cast<std::size_t>(pivot)]);

        const long double pv = a[static_cast<std::size_t>(col)][static_cast<std::size_t>(col)];
        for (int j = col; j <= m; ++j) {
            a[static_cast<std::size_t>(col)][static_cast<std::size_t>(j)] /= pv;
        }

        for (int r = 0; r < m; ++r) {
            if (r == col) {
                continue;
            }
            const long double factor = a[static_cast<std::size_t>(r)][static_cast<std::size_t>(col)];
            if (std::fabsl(factor) < 1e-30L) {
                continue;
            }
            for (int j = col; j <= m; ++j) {
                a[static_cast<std::size_t>(r)][static_cast<std::size_t>(j)] -=
                    factor * a[static_cast<std::size_t>(col)][static_cast<std::size_t>(j)];
            }
        }
    }

    std::vector<long double> h(static_cast<std::size_t>(ccount), 0.0L);
    for (int ri = 0; ri < m; ++ri) {
        const int cls = unknown_classes[static_cast<std::size_t>(ri)];
        h[static_cast<std::size_t>(cls)] = a[static_cast<std::size_t>(ri)][static_cast<std::size_t>(m)];
    }
    h[static_cast<std::size_t>(id_idx)] = 0.0L;

    const u64 total_perm = factorial(n);
    long double weighted_sum = 0.0L;
    for (int i = 0; i < ccount; ++i) {
        const u64 sz = class_size(n, partitions[static_cast<std::size_t>(i)]);
        weighted_sum += static_cast<long double>(sz) * h[static_cast<std::size_t>(i)];
    }

    return weighted_sum / static_cast<long double>(total_perm);
}

bool run_checkpoints() {
    const long double avg4 = solve_average_expectation(4);
    if (std::fabsl(avg4 - 27.5L) > 1e-10L) {
        std::cerr << "Checkpoint failed: n=4 average expectation\n";
        return false;
    }
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
        return 2;
    }

    const long double avg11 = solve_average_expectation(11);
    const long long answer = std::llround(avg11);
    std::cout << answer << '\n';
    return 0;
}
