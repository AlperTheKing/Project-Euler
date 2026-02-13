#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr std::array<int, 5> kCardChoose = {1, 4, 6, 4, 1};
constexpr std::array<std::array<int, 5>, 5> kSubChoose = {
    std::array<int, 5>{1, 0, 0, 0, 0},
    std::array<int, 5>{1, 1, 0, 0, 0},
    std::array<int, 5>{1, 2, 1, 0, 0},
    std::array<int, 5>{1, 3, 3, 1, 0},
    std::array<int, 5>{1, 4, 6, 4, 1},
};

struct ScoreParts {
    int hand_score = 0;
    int pair_score = 0;
    int run_score = 0;
    int fifteen_count = 0;
};

ScoreParts evaluate_hand(const std::vector<int>& counts, const std::vector<int>& values) {
    const int n = static_cast<int>(counts.size());

    ScoreParts out;

    for (int i = 0; i < n; ++i) {
        const int c = counts[i];
        out.hand_score += values[i] * c;
        out.pair_score += c * (c - 1);
    }

    int len = 0;
    int prod = 1;
    for (int i = 0; i < n; ++i) {
        const int c = counts[i];
        if (c == 0) {
            if (len >= 3) {
                out.run_score += len * prod;
            }
            len = 0;
            prod = 1;
        } else {
            ++len;
            prod *= c;
        }
    }
    if (len >= 3) {
        out.run_score += len * prod;
    }

    std::array<int, 16> ways{};
    ways[0] = 1;
    for (int i = 0; i < n; ++i) {
        const int c = counts[i];
        const int v = values[i];

        std::array<int, 16> next{};
        for (int s = 0; s <= 15; ++s) {
            if (ways[s] == 0) {
                continue;
            }
            for (int d = 0; d <= c; ++d) {
                const int t = s + d * v;
                if (t > 15) {
                    break;
                }
                next[t] += ways[s] * kSubChoose[c][d];
            }
        }
        ways = next;
    }
    out.fifteen_count = ways[15];

    return out;
}

u64 brute_count_equal(const std::vector<int>& values) {
    const int n = static_cast<int>(values.size());
    std::vector<int> counts(n, 0);

    u64 total = 0;

    auto rec = [&](auto&& self, int idx, int cards, u64 weight) -> void {
        if (idx == n) {
            if (cards == 0) {
                return;
            }
            const ScoreParts sc = evaluate_hand(counts, values);
            const int crib = sc.pair_score + sc.run_score + 2 * sc.fifteen_count;
            if (sc.hand_score == crib) {
                total += weight;
            }
            return;
        }

        for (int c = 0; c <= 4; ++c) {
            counts[idx] = c;
            self(self, idx + 1, cards + c, weight * static_cast<u64>(kCardChoose[c]));
        }
    };

    rec(rec, 0, 0, 1);
    return total;
}

struct HalfState {
    int hand_score = 0;
    int pair_score = 0;
    int run_score = 0;
    int base = 0;

    int prefix_len = 0;
    int prefix_prod = 1;
    int prefix_score = 0;

    int suffix_len = 0;
    int suffix_prod = 1;
    int suffix_score = 0;

    std::array<std::uint32_t, 16> poly{};
    std::array<std::uint32_t, 16> poly_rev{};

    u64 hand_multiplicity = 0;
};

HalfState make_half_state(const std::vector<int>& counts, const std::vector<int>& values) {
    const int n = static_cast<int>(counts.size());
    HalfState st;

    st.hand_multiplicity = 1;

    for (int i = 0; i < n; ++i) {
        const int c = counts[i];
        st.hand_score += values[i] * c;
        st.pair_score += c * (c - 1);
        st.hand_multiplicity *= static_cast<u64>(kCardChoose[c]);
    }

    int len = 0;
    int prod = 1;
    for (int i = 0; i < n; ++i) {
        const int c = counts[i];
        if (c == 0) {
            if (len >= 3) {
                st.run_score += len * prod;
            }
            len = 0;
            prod = 1;
        } else {
            ++len;
            prod *= c;
        }
    }
    if (len >= 3) {
        st.run_score += len * prod;
    }

    st.prefix_len = 0;
    st.prefix_prod = 1;
    while (st.prefix_len < n && counts[st.prefix_len] > 0) {
        st.prefix_prod *= counts[st.prefix_len];
        ++st.prefix_len;
    }
    st.prefix_score = (st.prefix_len >= 3) ? (st.prefix_len * st.prefix_prod) : 0;

    st.suffix_len = 0;
    st.suffix_prod = 1;
    int i = n - 1;
    while (i >= 0 && counts[i] > 0) {
        st.suffix_prod *= counts[i];
        ++st.suffix_len;
        --i;
    }
    st.suffix_score = (st.suffix_len >= 3) ? (st.suffix_len * st.suffix_prod) : 0;

    std::array<int, 16> ways{};
    ways[0] = 1;

    for (int j = 0; j < n; ++j) {
        const int c = counts[j];
        const int v = values[j];

        std::array<int, 16> next{};
        for (int s = 0; s <= 15; ++s) {
            if (ways[s] == 0) {
                continue;
            }
            for (int d = 0; d <= c; ++d) {
                const int t = s + d * v;
                if (t > 15) {
                    break;
                }
                next[t] += ways[s] * kSubChoose[c][d];
            }
        }
        ways = next;
    }

    for (int s = 0; s <= 15; ++s) {
        st.poly[s] = static_cast<std::uint32_t>(ways[s]);
        st.poly_rev[s] = static_cast<std::uint32_t>(ways[15 - s]);
    }

    st.base = st.hand_score - st.pair_score - st.run_score;

    return st;
}

std::vector<HalfState> enumerate_half_states(const std::vector<int>& values) {
    const int n = static_cast<int>(values.size());
    std::size_t total = 1;
    for (int i = 0; i < n; ++i) {
        total *= 5;
    }

    std::vector<HalfState> states;
    states.reserve(total);

    std::vector<int> counts(n, 0);

    auto rec = [&](auto&& self, int idx) -> void {
        if (idx == n) {
            states.push_back(make_half_state(counts, values));
            return;
        }
        for (int c = 0; c <= 4; ++c) {
            counts[idx] = c;
            self(self, idx + 1);
        }
    };

    rec(rec, 0);
    return states;
}

u64 count_equal_mitm(const std::vector<int>& values) {
    const int n = static_cast<int>(values.size());
    const int split = n / 2;

    const std::vector<int> left_values(values.begin(), values.begin() + split);
    const std::vector<int> right_values(values.begin() + split, values.end());

    const std::vector<HalfState> left = enumerate_half_states(left_values);
    const std::vector<HalfState> right = enumerate_half_states(right_values);

    const unsigned threads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<u64> local(threads, 0);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (unsigned tid = 0; tid < threads; ++tid) {
        const std::size_t begin = (left.size() * tid) / threads;
        const std::size_t end = (left.size() * (tid + 1)) / threads;

        workers.emplace_back([&, tid, begin, end]() {
            u64 subtotal = 0;

            for (std::size_t i = begin; i < end; ++i) {
                const HalfState& L = left[i];

                for (const HalfState& R : right) {
                    int adjust = 0;
                    if (L.suffix_len > 0 && R.prefix_len > 0) {
                        const int merged_len = L.suffix_len + R.prefix_len;
                        const int merged_prod = L.suffix_prod * R.prefix_prod;
                        const int merged_score = (merged_len >= 3) ? (merged_len * merged_prod) : 0;
                        adjust = merged_score - L.suffix_score - R.prefix_score;
                    }

                    const int target_twice = L.base + R.base - adjust;
                    if (target_twice < 0 || target_twice > 340 || (target_twice & 1) != 0) {
                        continue;
                    }

                    const int target_f = target_twice / 2;
                    i64 f = 0;

                    for (int k = 0; k <= 15; ++k) {
                        f += static_cast<i64>(L.poly[k]) * static_cast<i64>(R.poly_rev[k]);
                        if (f > target_f) {
                            break;
                        }
                    }

                    if (f == target_f) {
                        subtotal += L.hand_multiplicity * R.hand_multiplicity;
                    }
                }
            }

            local[tid] = subtotal;
        });
    }

    for (std::thread& t : workers) {
        t.join();
    }

    u64 total = 0;
    for (u64 v : local) {
        total += v;
    }

    return total - 1;
}

void run_validations() {
    {
        std::vector<int> values = {1, 2, 3, 4, 5, 6, 7, 8};
        assert(count_equal_mitm(values) == brute_count_equal(values));
    }

    {
        std::vector<int> values(13);
        for (int i = 0; i < 9; ++i) {
            values[i] = i + 1;
        }
        for (int i = 9; i < 13; ++i) {
            values[i] = 10;
        }

        std::vector<int> sample1(13, 0);
        sample1[4] = 3;
        sample1[12] = 1;
        const ScoreParts s1 = evaluate_hand(sample1, values);
        assert(s1.pair_score + s1.run_score + 2 * s1.fifteen_count == 14);

        std::vector<int> sample2(13, 0);
        sample2[0] = 2;
        sample2[1] = 1;
        sample2[2] = 1;
        sample2[3] = 1;
        sample2[4] = 1;
        const ScoreParts s2 = evaluate_hand(sample2, values);
        assert(s2.hand_score == 16);
        assert(s2.pair_score + s2.run_score + 2 * s2.fifteen_count == 16);
    }
}

}  // namespace

int main() {
    run_validations();

    std::vector<int> values(13);
    for (int i = 0; i < 9; ++i) {
        values[i] = i + 1;
    }
    for (int i = 9; i < 13; ++i) {
        values[i] = 10;
    }

    std::cout << count_equal_mitm(values) << '\n';
    return 0;
}
