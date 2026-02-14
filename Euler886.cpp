#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <functional>
#include <map>
#include <numeric>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

constexpr int MOD = 83'456'729;

int add_mod(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}

int sub_mod(int a, int b) {
    int s = a - b;
    if (s < 0) s += MOD;
    return s;
}

int mul_mod(i64 a, i64 b) {
    return static_cast<int>((static_cast<__int128>(a) * b) % MOD);
}

int mod_pow(int a, i64 e) {
    i64 r = 1;
    i64 x = a % MOD;
    if (x < 0) x += MOD;
    while (e > 0) {
        if (e & 1LL) r = (r * x) % MOD;
        x = (x * x) % MOD;
        e >>= 1LL;
    }
    return static_cast<int>(r);
}

int mod_inv(int a) {
    return mod_pow(a, MOD - 2);
}

std::vector<std::vector<int>> build_binom(int n) {
    std::vector<std::vector<int>> C(static_cast<std::size_t>(n + 1));
    C[0] = {1};
    for (int i = 1; i <= n; ++i) {
        C[static_cast<std::size_t>(i)].assign(static_cast<std::size_t>(i + 1), 1);
        for (int j = 1; j < i; ++j) {
            C[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                C[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j - 1)] +
                C[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j)];
        }
    }
    return C;
}

int determinant_mod(std::vector<std::vector<int>> A) {
    const int n = static_cast<int>(A.size());
    if (n == 0) return 1;
    if (n == 1) return A[0][0] % MOD;

    int result_inverse = 1;

    for (int i = 0; i < n - 2; ++i) {
        int row_to_use = -1;
        for (int r = i; r < n; ++r) {
            if (A[r][i] % MOD != 0) {
                row_to_use = r;
                break;
            }
        }
        if (row_to_use < 0) return 0;

        if (i != row_to_use) {
            std::swap(A[row_to_use], A[i]);
            result_inverse = (result_inverse == 0 ? 0 : MOD - result_inverse);
        }

        const int Aii = A[i][i] % MOD;
        result_inverse = mul_mod(result_inverse, mod_pow(Aii, n - i - 2));

        for (int r = i + 1; r < n; ++r) {
            const int x = A[r][i] % MOD;
            A[r][i] = 0;
            for (int c = i + 1; c < n; ++c) {
                int v = sub_mod(mul_mod(Aii, A[r][c]), mul_mod(x, A[i][c]));
                A[r][c] = v;
            }
        }
    }

    int last = sub_mod(mul_mod(A[n - 2][n - 2], A[n - 1][n - 1]),
                       mul_mod(A[n - 2][n - 1], A[n - 1][n - 2]));
    return mul_mod(last, mod_inv(result_inverse));
}

int permanent_mod_with_column_classes(const std::vector<std::vector<int>>& A,
                                      const std::vector<std::vector<int>>& binom) {
    const int n = static_cast<int>(A.size());
    if (n == 0) return 1;

    std::map<u64, std::vector<int>> column_classes;
    for (int col = 0; col < n; ++col) {
        u64 mask = 0;
        for (int r = 0; r < n; ++r) {
            if (A[r][col] != 0) mask |= (1ULL << r);
        }
        column_classes[mask].push_back(col);
    }

    std::vector<std::vector<int>> classes;
    for (const auto& [k, inds] : column_classes) {
        (void)k;
        if (static_cast<int>(inds.size()) >= 3) {
            classes.push_back(inds);
        } else {
            for (int idx : inds) {
                classes.push_back({idx});
            }
        }
    }
    std::sort(classes.begin(), classes.end(), [](const auto& a, const auto& b) {
        return a.size() < b.size();
    });

    std::vector<int> new_col_order;
    new_col_order.reserve(n);
    for (const auto& c : classes) {
        for (int idx : c) new_col_order.push_back(idx);
    }

    std::vector<std::vector<int>> B(n, std::vector<int>(n));
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            B[r][c] = A[r][new_col_order[c]];
        }
    }

    std::vector<int> class_sizes;
    class_sizes.reserve(classes.size());
    for (const auto& c : classes) class_sizes.push_back(static_cast<int>(c.size()));

    const int default_block_exponent = 10;
    int smallest_class_count = 0;
    int smallest_class_index_count = 0;
    while (smallest_class_count < static_cast<int>(class_sizes.size()) &&
           smallest_class_index_count + class_sizes[smallest_class_count] < default_block_exponent) {
        smallest_class_index_count += class_sizes[smallest_class_count];
        ++smallest_class_count;
    }

    int block_exponent = smallest_class_index_count;
    int block_size = 1 << block_exponent;

    std::vector<int> high_class_sizes;
    std::vector<int> high_class_offsets;
    int off = block_exponent;
    for (int i = smallest_class_count; i < static_cast<int>(class_sizes.size()); ++i) {
        high_class_sizes.push_back(class_sizes[i]);
        high_class_offsets.push_back(off);
        off += class_sizes[i];
    }

    std::vector<int> result_block(static_cast<std::size_t>(block_size), 0);
    std::vector<int> high_mask(n, 0);

    std::vector<int> popcnt(static_cast<std::size_t>(block_size), 0);
    for (int s = 1; s < block_size; ++s) popcnt[static_cast<std::size_t>(s)] = popcnt[s >> 1] + (s & 1);

    std::function<void(int, int, int)> rec = [&](int cls_idx, int total_ways, int high_ones) {
        if (cls_idx == static_cast<int>(high_class_sizes.size())) {
            for (int state = 0; state < block_size; ++state) {
                int term = total_ways;
                for (int r = 0; r < n; ++r) {
                    int sum = 0;
                    for (int c = 0; c < block_exponent; ++c) {
                        if ((state >> c) & 1) sum += B[r][c];
                    }
                    for (int c = block_exponent; c < n; ++c) {
                        if (high_mask[c]) sum += B[r][c];
                    }
                    term = mul_mod(term, sum);
                    if (term == 0) break;
                }

                if (((popcnt[static_cast<std::size_t>(state)] + high_ones) & 1) != 0 && term != 0) {
                    term = MOD - term;
                }
                result_block[static_cast<std::size_t>(state)] =
                    add_mod(result_block[static_cast<std::size_t>(state)], term);
            }
            return;
        }

        const int sz = high_class_sizes[static_cast<std::size_t>(cls_idx)];
        const int base = high_class_offsets[static_cast<std::size_t>(cls_idx)];

        for (int ones = 0; ones <= sz; ++ones) {
            for (int j = 0; j < sz; ++j) {
                high_mask[base + j] = (j >= sz - ones) ? 1 : 0;
            }
            int ways2 = mul_mod(total_ways, binom[static_cast<std::size_t>(sz)][static_cast<std::size_t>(ones)]);
            rec(cls_idx + 1, ways2, high_ones + ones);
        }
    };

    rec(0, 1, 0);

    int result = 0;
    for (int v : result_block) {
        result = add_mod(result, v);
    }
    if (n & 1) {
        result = (result == 0) ? 0 : (MOD - result);
    }
    return result;
}

std::vector<std::pair<int, int>> row_classes(const std::vector<std::vector<int>>& A) {
    std::map<u64, std::vector<int>> groups;
    for (int i = 0; i < static_cast<int>(A.size()); ++i) {
        u64 mask = 0;
        for (int j = 0; j < static_cast<int>(A[i].size()); ++j) {
            if (A[i][j] != 0) mask |= (1ULL << j);
        }
        groups[mask].push_back(i);
    }

    std::vector<std::pair<int, int>> out;
    out.reserve(groups.size());
    for (const auto& [k, idxs] : groups) {
        (void)k;
        out.push_back({idxs[0], static_cast<int>(idxs.size())});
    }
    return out;
}

template <class F>
void for_each_combination(int n, int k, F&& fn) {
    if (k < 0 || k > n) return;
    if (k == 0) {
        std::vector<int> empty;
        fn(empty);
        return;
    }
    std::vector<int> comb(static_cast<std::size_t>(k));
    for (int i = 0; i < k; ++i) comb[static_cast<std::size_t>(i)] = i;

    while (true) {
        fn(comb);
        int i = k - 1;
        while (i >= 0 && comb[static_cast<std::size_t>(i)] == n - k + i) --i;
        if (i < 0) break;
        ++comb[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < k; ++j) {
            comb[static_cast<std::size_t>(j)] = comb[static_cast<std::size_t>(j - 1)] + 1;
        }
    }
}

int P_mod(int n) {
    assert(n >= 2 && (n % 2 == 0));
    const int half_n = n / 2;

    std::vector<int> even_indices, odd_indices;
    even_indices.reserve(half_n);
    odd_indices.reserve(half_n);

    for (int x = 1; x < n; x += 2) even_indices.push_back(x);
    for (int x = 2; x <= n; x += 2) odd_indices.push_back(x);

    std::vector<std::vector<int>> B(half_n, std::vector<int>(half_n, 0));
    for (int c = 0; c < half_n; ++c) {
        for (int r = 0; r < half_n; ++r) {
            B[c][r] = (std::gcd(even_indices[r], odd_indices[c]) == 1) ? 1 : 0;
        }
    }

    std::vector<std::vector<int>> binom = build_binom(half_n);

    std::vector<std::pair<int, int>> rows_and_ways = row_classes(B);
    if (!rows_and_ways.empty()) {
        int row_to_remove = -1;
        for (int i = 0; i < static_cast<int>(rows_and_ways.size()); ++i) {
            if (rows_and_ways[static_cast<std::size_t>(i)].second == 1) {
                row_to_remove = i;
                break;
            }
        }
        if (row_to_remove >= 0) {
            rows_and_ways.erase(rows_and_ways.begin() + row_to_remove);
        } else {
            rows_and_ways[0].second -= 1;
        }
    }

    std::vector<std::vector<int>> BT(half_n, std::vector<int>(half_n, 0));
    for (int r = 0; r < half_n; ++r) {
        for (int c = 0; c < half_n; ++c) {
            BT[r][c] = B[c][r];
        }
    }
    std::vector<std::pair<int, int>> cols_and_ways = row_classes(BT);

    int result = 0;

    for (int included = 0; included < half_n; ++included) {
        for_each_combination(static_cast<int>(rows_and_ways.size()), included,
                             [&](const std::vector<int>& rows_pick) {
                                 i64 rows_ways = 1;
                                 std::vector<int> included_rows;
                                 included_rows.reserve(rows_pick.size());
                                 for (int idx : rows_pick) {
                                     const auto [row_idx, ways] = rows_and_ways[static_cast<std::size_t>(idx)];
                                     rows_ways *= ways;
                                     included_rows.push_back(row_idx);
                                 }

                                 std::vector<char> in_row(half_n, 0);
                                 for (int r : included_rows) in_row[r] = 1;
                                 std::vector<int> not_rows;
                                 not_rows.reserve(half_n - included);
                                 for (int r = 0; r < half_n; ++r) {
                                     if (!in_row[r]) not_rows.push_back(r);
                                 }

                                 for_each_combination(static_cast<int>(cols_and_ways.size()), included,
                                                      [&](const std::vector<int>& cols_pick) {
                                                          i64 ways = rows_ways;
                                                          std::vector<int> included_cols;
                                                          included_cols.reserve(cols_pick.size());
                                                          for (int idx : cols_pick) {
                                                              const auto [col_idx, cw] =
                                                                  cols_and_ways[static_cast<std::size_t>(idx)];
                                                              ways *= cw;
                                                              included_cols.push_back(col_idx);
                                                          }

                                                          std::vector<std::vector<int>> det_mat(
                                                              static_cast<std::size_t>(included),
                                                              std::vector<int>(static_cast<std::size_t>(included),
                                                                               0));
                                                          for (int i = 0; i < included; ++i) {
                                                              for (int j = 0; j < included; ++j) {
                                                                  det_mat[static_cast<std::size_t>(i)]
                                                                         [static_cast<std::size_t>(j)] =
                                                                      B[included_rows[static_cast<std::size_t>(i)]]
                                                                       [included_cols[static_cast<std::size_t>(j)]];
                                                              }
                                                          }

                                                          int det_val = determinant_mod(det_mat);
                                                          if (det_val == 0) return;

                                                          det_val = mul_mod(det_val, det_val);
                                                          if (included & 1) {
                                                              det_val = (det_val == 0) ? 0 : (MOD - det_val);
                                                          }

                                                          std::vector<char> in_col(half_n, 0);
                                                          for (int c : included_cols) in_col[c] = 1;
                                                          std::vector<int> not_cols;
                                                          not_cols.reserve(half_n - included);
                                                          for (int c = 0; c < half_n; ++c) {
                                                              if (!in_col[c]) not_cols.push_back(c);
                                                          }

                                                          const int rest = half_n - included;
                                                          std::vector<std::vector<int>> per_mat(
                                                              static_cast<std::size_t>(rest),
                                                              std::vector<int>(static_cast<std::size_t>(rest), 0));
                                                          for (int i = 0; i < rest; ++i) {
                                                              for (int j = 0; j < rest; ++j) {
                                                                  per_mat[static_cast<std::size_t>(i)]
                                                                         [static_cast<std::size_t>(j)] =
                                                                      B[not_rows[static_cast<std::size_t>(i)]]
                                                                       [not_cols[static_cast<std::size_t>(j)]];
                                                              }
                                                          }

                                                          int per_val = permanent_mod_with_column_classes(per_mat, binom);
                                                          per_val = mul_mod(per_val, per_val);

                                                          int term = mul_mod(static_cast<int>(ways % MOD),
                                                                             mul_mod(det_val, per_val));
                                                          result = add_mod(result, term);
                                                      });
                             });
    }

    return result;
}

int brute_count(int n) {
    std::vector<int> v;
    for (int x = 2; x <= n; ++x) v.push_back(x);

    int cnt = 0;
    std::sort(v.begin(), v.end());
    do {
        bool ok = true;
        for (int i = 1; i < static_cast<int>(v.size()); ++i) {
            if (std::gcd(v[i - 1], v[i]) != 1) {
                ok = false;
                break;
            }
        }
        if (ok) ++cnt;
    } while (std::next_permutation(v.begin(), v.end()));

    return cnt;
}

}  // namespace

int main() {
    assert(P_mod(4) == 2);
    assert(P_mod(10) == 576);
    assert(P_mod(6) == brute_count(6));
    assert(P_mod(8) == brute_count(8));

    std::cout << P_mod(34) << '\n';
    return 0;
}
