#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <vector>

using int64 = long long;

namespace {

constexpr int kTargetSum = 10;
constexpr int kDigits = 10;
constexpr int64 kMod = 1000000007LL;

struct State {
    uint16_t history = 0;  // bits for reachable prefix sums in [S_c-10, S_c]
    uint16_t segment = 0;  // bits for segment sums from c to current position
};

struct Automaton {
    int start = 0;
    std::vector<std::array<int, kDigits>> trans;
    std::vector<char> accept;
};

constexpr uint32_t kHistoryMask = (1u << (kTargetSum + 1)) - 1u;

uint32_t encode_state(uint16_t history, uint16_t segment) {
    return static_cast<uint32_t>(history) | (static_cast<uint32_t>(segment) << 11);
}

int highest_bit(uint16_t mask) {
    return 31 - __builtin_clz(static_cast<unsigned>(mask));
}

bool next_state(uint16_t history, uint16_t segment, int digit, State& out) {
    int delta = highest_bit(segment);
    int new_delta = delta + digit;
    if (new_delta > kTargetSum) return false;
    uint16_t new_segment = static_cast<uint16_t>(segment | (1u << new_delta));
    if ((history >> (kTargetSum - new_delta)) & 1u) {
        uint32_t shifted = (static_cast<uint32_t>(history) << new_delta) & kHistoryMask;
        uint16_t new_history = static_cast<uint16_t>(shifted);
        for (int k = 0; k <= new_delta; ++k) {
            if ((new_segment >> k) & 1u) {
                new_history |= static_cast<uint16_t>(1u << (new_delta - k));
            }
        }
        out.history = new_history;
        out.segment = 1u;
    } else {
        out.history = history;
        out.segment = new_segment;
    }
    return true;
}

struct RawAutomaton {
    std::vector<State> states;
    std::vector<std::array<int, kDigits>> trans;
    std::vector<char> accept;
    int start = 0;
};

RawAutomaton build_raw() {
    RawAutomaton raw;
    raw.states.reserve(20000);
    std::unordered_map<uint32_t, int> index;
    index.reserve(20000);

    State start{1u, 1u};
    raw.states.push_back(start);
    index[encode_state(start.history, start.segment)] = 0;

    for (size_t i = 0; i < raw.states.size(); ++i) {
        const State& cur = raw.states[i];
        for (int d = 0; d < kDigits; ++d) {
            State nxt;
            if (!next_state(cur.history, cur.segment, d, nxt)) continue;
            uint32_t key = encode_state(nxt.history, nxt.segment);
            auto it = index.find(key);
            if (it == index.end()) {
                int id = static_cast<int>(raw.states.size());
                index[key] = id;
                raw.states.push_back(nxt);
            }
        }
    }

    int n = static_cast<int>(raw.states.size());
    raw.trans.assign(n, {});
    raw.accept.assign(n, 0);
    for (int i = 0; i < n; ++i) {
        const State& cur = raw.states[i];
        raw.accept[i] = (cur.segment == 1u) ? 1 : 0;
        for (int d = 0; d < kDigits; ++d) {
            State nxt;
            if (!next_state(cur.history, cur.segment, d, nxt)) {
                raw.trans[i][d] = -1;
            } else {
                uint32_t key = encode_state(nxt.history, nxt.segment);
                raw.trans[i][d] = index[key];
            }
        }
    }
    raw.start = 0;
    return raw;
}

Automaton minimize_automaton(const RawAutomaton& raw) {
    int n = static_cast<int>(raw.states.size());
    int dead = n;
    int total = n + 1;

    std::vector<std::array<int, kDigits>> trans(total);
    std::vector<char> accept(total, 0);
    for (int i = 0; i < n; ++i) {
        accept[i] = raw.accept[i];
        for (int d = 0; d < kDigits; ++d) {
            int t = raw.trans[i][d];
            trans[i][d] = (t == -1) ? dead : t;
        }
    }
    accept[dead] = 0;
    trans[dead].fill(dead);

    std::vector<std::vector<int>> inv[kDigits];
    for (int d = 0; d < kDigits; ++d) {
        inv[d].assign(total, {});
    }
    for (int s = 0; s < total; ++s) {
        for (int d = 0; d < kDigits; ++d) {
            inv[d][trans[s][d]].push_back(s);
        }
    }

    std::vector<std::vector<int>> classes;
    classes.reserve(total);
    classes.push_back({});
    classes.push_back({});
    for (int s = 0; s < total; ++s) {
        classes[accept[s] ? 0 : 1].push_back(s);
    }
    if (classes[0].empty()) classes.erase(classes.begin());
    if (classes.size() == 2 && classes[1].empty()) classes.pop_back();

    std::vector<int> class_of(total, -1);
    for (int i = 0; i < static_cast<int>(classes.size()); ++i) {
        for (int s : classes[i]) class_of[s] = i;
    }

    std::deque<int> work;
    std::vector<char> in_work(classes.size(), 0);
    for (int i = 0; i < static_cast<int>(classes.size()); ++i) {
        work.push_back(i);
        in_work[i] = 1;
    }

    std::vector<char> marked(total, 0);
    std::vector<std::vector<int>> bucket(classes.size());

    while (!work.empty()) {
        int A = work.front();
        work.pop_front();
        in_work[A] = 0;
        for (int d = 0; d < kDigits; ++d) {
            std::vector<int> X;
            X.reserve(total / 4);
            for (int t : classes[A]) {
                for (int s : inv[d][t]) {
                    if (!marked[s]) {
                        marked[s] = 1;
                        X.push_back(s);
                    }
                }
            }
            if (X.empty()) continue;
            std::vector<int> touched;
            for (int s : X) {
                int c = class_of[s];
                if (bucket[c].empty()) touched.push_back(c);
                bucket[c].push_back(s);
            }
            for (int c : touched) {
                if (bucket[c].size() < classes[c].size()) {
                    std::vector<int> inX;
                    inX.swap(bucket[c]);
                    std::vector<int> outX;
                    outX.reserve(classes[c].size() - inX.size());
                    for (int s : classes[c]) {
                        if (!marked[s]) outX.push_back(s);
                    }
                    classes[c].swap(inX);
                    int new_id = static_cast<int>(classes.size());
                    classes.push_back(std::move(outX));
                    if (bucket.size() < classes.size()) bucket.resize(classes.size());
                    if (in_work.size() < classes.size()) in_work.resize(classes.size(), 0);
                    for (int s : classes[new_id]) class_of[s] = new_id;
                    if (in_work[c]) {
                        if (!in_work[new_id]) {
                            work.push_back(new_id);
                            in_work[new_id] = 1;
                        }
                    } else {
                        int pick = (classes[c].size() <= classes[new_id].size()) ? c : new_id;
                        if (!in_work[pick]) {
                            work.push_back(pick);
                            in_work[pick] = 1;
                        }
                    }
                } else {
                    bucket[c].clear();
                }
            }
            for (int s : X) marked[s] = 0;
        }
    }

    Automaton min;
    min.trans.assign(classes.size(), {});
    min.accept.assign(classes.size(), 0);
    min.start = class_of[raw.start];
    for (int i = 0; i < static_cast<int>(classes.size()); ++i) {
        int rep = classes[i][0];
        min.accept[i] = accept[rep];
        for (int d = 0; d < kDigits; ++d) {
            min.trans[i][d] = class_of[trans[rep][d]];
        }
    }
    return min;
}

struct Sequences {
    std::vector<int64> counts;
    std::vector<int64> prefix;
};

Sequences build_sequences(const Automaton& aut, int terms) {
    int n = static_cast<int>(aut.trans.size());
    std::vector<int> accept_list;
    accept_list.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (aut.accept[i]) accept_list.push_back(i);
    }

    std::vector<int64> counts(terms + 1, 0);
    std::vector<int64> prefix(terms + 1, 0);
    std::vector<int64> cur(n, 0), next(n, 0);
    cur[aut.start] = 1;

    auto apply_digits = [&](int d_start, int d_end) {
        std::fill(next.begin(), next.end(), 0);
        for (int i = 0; i < n; ++i) {
            int64 val = cur[i];
            if (!val) continue;
            for (int d = d_start; d <= d_end; ++d) {
                int j = aut.trans[i][d];
                int64& slot = next[j];
                slot += val;
                if (slot >= kMod) slot -= kMod;
            }
        }
        cur.swap(next);
    };

    apply_digits(1, 9);
    int64 sum = 0;
    for (int idx : accept_list) {
        sum += cur[idx];
        if (sum >= kMod) sum -= kMod;
    }
    counts[1] = sum;
    prefix[1] = sum;

    for (int len = 2; len <= terms; ++len) {
        apply_digits(0, 9);
        sum = 0;
        for (int idx : accept_list) {
            sum += cur[idx];
            if (sum >= kMod) sum -= kMod;
        }
        counts[len] = sum;
        prefix[len] = prefix[len - 1] + sum;
        if (prefix[len] >= kMod) prefix[len] -= kMod;
    }
    return {counts, prefix};
}

int64 mod_pow(int64 base, uint64_t exp) {
    int64 result = 1 % kMod;
    int64 cur = base % kMod;
    while (exp > 0) {
        if (exp & 1) result = (result * cur) % kMod;
        cur = (cur * cur) % kMod;
        exp >>= 1;
    }
    return result;
}

int64 mod_inv(int64 x) {
    return mod_pow(x, kMod - 2);
}

std::vector<int64> berlekamp_massey(const std::vector<int64>& seq) {
    std::vector<int64> C(1, 1), B(1, 1);
    int L = 0;
    int m = 1;
    int64 b = 1;
    for (int n = 0; n < static_cast<int>(seq.size()); ++n) {
        int64 d = 0;
        for (int i = 0; i <= L; ++i) {
            d = (d + C[i] * seq[n - i]) % kMod;
        }
        if (d == 0) {
            ++m;
            continue;
        }
        std::vector<int64> T = C;
        int64 coef = d * mod_inv(b) % kMod;
        if (C.size() < B.size() + m) C.resize(B.size() + m, 0);
        for (size_t i = 0; i < B.size(); ++i) {
            int64 val = (coef * B[i]) % kMod;
            C[i + m] -= val;
            if (C[i + m] < 0) C[i + m] += kMod;
        }
        if (2 * L <= n) {
            L = n + 1 - L;
            B.swap(T);
            b = d;
            m = 1;
        } else {
            ++m;
        }
    }
    C.erase(C.begin());
    for (auto& x : C) x = (kMod - x) % kMod;
    return C;
}

std::vector<int64> combine_poly(const std::vector<int64>& a, const std::vector<int64>& b,
                                const std::vector<int64>& coef, int threads) {
    int k = static_cast<int>(coef.size());
    std::vector<int64> res(2 * k, 0);
    if (threads <= 1 || k < 512) {
        for (int i = 0; i < k; ++i) {
            if (a[i] == 0) continue;
            for (int j = 0; j < k; ++j) {
                if (b[j] == 0) continue;
                int64 add = (a[i] * b[j]) % kMod;
                int64& slot = res[i + j];
                slot += add;
                if (slot >= kMod) slot -= kMod;
            }
        }
    } else {
        int tcount = std::min(threads, k);
        std::vector<std::vector<int64>> partial(tcount, std::vector<int64>(2 * k, 0));
        std::vector<std::thread> pool;
        pool.reserve(tcount);
        for (int t = 0; t < tcount; ++t) {
            int start = t * k / tcount;
            int end = (t + 1) * k / tcount;
            pool.emplace_back([&, t, start, end]() {
                auto& loc = partial[t];
                for (int i = start; i < end; ++i) {
                    if (a[i] == 0) continue;
                    for (int j = 0; j < k; ++j) {
                        if (b[j] == 0) continue;
                        int64 add = (a[i] * b[j]) % kMod;
                        int64& slot = loc[i + j];
                        slot += add;
                        if (slot >= kMod) slot -= kMod;
                    }
                }
            });
        }
        for (auto& th : pool) th.join();
        for (int i = 0; i < 2 * k - 1; ++i) {
            int64 sum = 0;
            for (int t = 0; t < tcount; ++t) {
                sum += partial[t][i];
            }
            res[i] = sum % kMod;
        }
    }

    for (int i = 2 * k - 2; i >= k; --i) {
        int64 val = res[i];
        if (val == 0) continue;
        for (int j = 0; j < k; ++j) {
            int64 add = (val * coef[j]) % kMod;
            int64& slot = res[i - 1 - j];
            slot += add;
            if (slot >= kMod) slot -= kMod;
        }
    }
    res.resize(k);
    return res;
}

int64 linear_recurrence(const std::vector<int64>& coef, const std::vector<int64>& init,
                        uint64_t n, int threads) {
    int k = static_cast<int>(coef.size());
    if (k == 0) return 0;
    if (n < init.size()) return init[n];
    if (k == 1) {
        return init[0] * mod_pow(coef[0], n) % kMod;
    }
    std::vector<int64> pol(k, 0), e(k, 0);
    pol[0] = 1;
    e[1] = 1;
    while (n > 0) {
        if (n & 1) pol = combine_poly(pol, e, coef, threads);
        e = combine_poly(e, e, coef, threads);
        n >>= 1;
    }
    int64 res = 0;
    for (int i = 0; i < k; ++i) {
        res = (res + pol[i] * init[i]) % kMod;
    }
    return res;
}

bool run_validations(const Sequences& seq) {
    struct Check {
        int n;
        int64 expected;
    };
    const std::vector<Check> checks = {
        {1, 0},
        {2, 9},
        {5, 3492},
    };

    for (const auto& check : checks) {
        if (check.n >= static_cast<int>(seq.prefix.size())) {
            std::cerr << "Validation skipped for T(" << check.n << "): term not computed.\n";
            continue;
        }
        int64 got = seq.prefix[check.n];
        if (got != check.expected) {
            std::cerr << "Validation failed for T(" << check.n << "): got=" << got
                      << " expected=" << check.expected << "\n";
            return false;
        }
    }
    if (seq.counts.size() > 2 && seq.counts[2] != 9) {
        std::cerr << "Validation failed for length 2: got=" << seq.counts[2] << " expected=9\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    uint64_t n = 1000000000000000000ULL;
    if (argc > 1) n = static_cast<uint64_t>(std::strtoull(argv[1], nullptr, 10));

    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads == 0) threads = 4;
    if (argc > 2) threads = std::max(1, std::atoi(argv[2]));

    RawAutomaton raw = build_raw();
    Automaton aut = minimize_automaton(raw);

    int terms = 2 * static_cast<int>(aut.trans.size()) + 5;
    Sequences seq = build_sequences(aut, terms);
    if (!run_validations(seq)) return 1;

    std::vector<int64> prefix(seq.prefix.begin(), seq.prefix.end());
    std::vector<int64> coef = berlekamp_massey(prefix);
    std::vector<int64> init(prefix.begin(), prefix.begin() + static_cast<int>(coef.size()));
    int64 answer = linear_recurrence(coef, init, n, threads);
    std::cout << answer % kMod << "\n";
    return 0;
}
