#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Guess {
    std::string digits;
    int matches = 0;
};

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

struct State {
    std::vector<int> assign;  // -1 or [0..9]
    std::vector<std::array<bool, 10>> banned;
};

bool propagate_singletons(State& st) {
    bool changed = true;
    while (changed) {
        changed = false;
        for (std::size_t pos = 0; pos < st.assign.size(); ++pos) {
            if (st.assign[pos] != -1) {
                continue;
            }
            int last_digit = -1;
            int allowed = 0;
            for (int d = 0; d <= 9; ++d) {
                if (!st.banned[pos][static_cast<std::size_t>(d)]) {
                    ++allowed;
                    last_digit = d;
                }
            }
            if (allowed == 0) {
                return false;
            }
            if (allowed == 1) {
                st.assign[pos] = last_digit;
                for (int d = 0; d <= 9; ++d) {
                    st.banned[pos][static_cast<std::size_t>(d)] = (d != last_digit);
                }
                changed = true;
            }
        }
    }
    return true;
}

bool apply_equal(State& st, const int pos, const int digit) {
    if (st.assign[static_cast<std::size_t>(pos)] == digit) {
        return true;
    }
    if (st.assign[static_cast<std::size_t>(pos)] != -1) {
        return false;
    }
    if (st.banned[static_cast<std::size_t>(pos)][static_cast<std::size_t>(digit)]) {
        return false;
    }
    st.assign[static_cast<std::size_t>(pos)] = digit;
    for (int d = 0; d <= 9; ++d) {
        st.banned[static_cast<std::size_t>(pos)][static_cast<std::size_t>(d)] = (d != digit);
    }
    return true;
}

bool apply_not_equal(State& st, const int pos, const int digit) {
    if (st.assign[static_cast<std::size_t>(pos)] == digit) {
        return false;
    }
    if (st.assign[static_cast<std::size_t>(pos)] != -1) {
        return true;
    }
    st.banned[static_cast<std::size_t>(pos)][static_cast<std::size_t>(digit)] = true;
    return true;
}

bool feasible(const State& st, const std::vector<Guess>& guesses) {
    for (const auto& g : guesses) {
        int min_match = 0;
        int max_match = 0;
        for (std::size_t i = 0; i < st.assign.size(); ++i) {
            const int d = g.digits[i] - '0';
            if (st.assign[i] == d) {
                ++min_match;
                ++max_match;
            } else if (st.assign[i] == -1 && !st.banned[i][static_cast<std::size_t>(d)]) {
                ++max_match;
            }
        }
        if (min_match > g.matches || max_match < g.matches) {
            return false;
        }
    }

    for (std::size_t i = 0; i < st.assign.size(); ++i) {
        if (st.assign[i] != -1) {
            continue;
        }
        bool any = false;
        for (int d = 0; d <= 9; ++d) {
            if (!st.banned[i][static_cast<std::size_t>(d)]) {
                any = true;
                break;
            }
        }
        if (!any) {
            return false;
        }
    }
    return true;
}

std::string assignment_to_string(const State& st) {
    std::string s(st.assign.size(), '0');
    for (std::size_t i = 0; i < st.assign.size(); ++i) {
        if (st.assign[i] < 0) {
            return "";
        }
        s[i] = static_cast<char>('0' + st.assign[i]);
    }
    return s;
}

bool satisfies_all(const std::string& candidate, const std::vector<Guess>& guesses) {
    for (const auto& g : guesses) {
        int match = 0;
        for (std::size_t i = 0; i < candidate.size(); ++i) {
            if (candidate[i] == g.digits[i]) {
                ++match;
            }
        }
        if (match != g.matches) {
            return false;
        }
    }
    return true;
}

long long combination_count(const int n, const int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    if (k == 0 || k == n) {
        return 1;
    }
    const int kk = std::min(k, n - k);
    long long value = 1;
    for (int i = 1; i <= kk; ++i) {
        value = value * (n - kk + i) / i;
    }
    return value;
}

void branch_guess(const std::vector<Guess>& guesses,
                  const int guess_idx,
                  const std::vector<int>& candidates,
                  const int need,
                  int idx,
                  std::vector<int>& chosen,
                  const State& base,
                  std::string& first_solution,
                  int& solutions);

void dfs(const std::vector<Guess>& guesses, const State& st, std::string& first_solution, int& solutions) {
    if (solutions >= 2) {
        return;
    }

    State cur = st;
    if (!propagate_singletons(cur) || !feasible(cur, guesses)) {
        return;
    }

    bool done = true;
    for (int d : cur.assign) {
        if (d < 0) {
            done = false;
            break;
        }
    }
    if (done) {
        const std::string candidate = assignment_to_string(cur);
        if (!candidate.empty() && satisfies_all(candidate, guesses)) {
            ++solutions;
            if (solutions == 1) {
                first_solution = candidate;
            }
        }
        return;
    }

    int best_guess = -1;
    long long best_branching = (1LL << 60);
    std::vector<int> best_candidates;
    int best_need = 0;

    for (std::size_t gi = 0; gi < guesses.size(); ++gi) {
        const Guess& g = guesses[gi];
        int fixed_match = 0;
        std::vector<int> candidates;
        for (std::size_t pos = 0; pos < cur.assign.size(); ++pos) {
            const int digit = g.digits[pos] - '0';
            if (cur.assign[pos] == digit) {
                ++fixed_match;
            } else if (cur.assign[pos] == -1 && !cur.banned[pos][static_cast<std::size_t>(digit)]) {
                candidates.push_back(static_cast<int>(pos));
            }
        }
        const int need = g.matches - fixed_match;
        if (need < 0 || need > static_cast<int>(candidates.size())) {
            return;
        }
        if (need == 0 && candidates.empty()) {
            continue;
        }
        const long long branches = combination_count(static_cast<int>(candidates.size()), need);
        if (branches < best_branching) {
            best_branching = branches;
            best_guess = static_cast<int>(gi);
            best_candidates = std::move(candidates);
            best_need = need;
        }
    }

    if (best_guess < 0) {
        return;
    }

    std::vector<int> chosen;
    branch_guess(guesses, best_guess, best_candidates, best_need, 0, chosen,
                 cur, first_solution, solutions);
}

void branch_guess(const std::vector<Guess>& guesses,
                  const int guess_idx,
                  const std::vector<int>& candidates,
                  const int need,
                  int idx,
                  std::vector<int>& chosen,
                  const State& base,
                  std::string& first_solution,
                  int& solutions) {
    if (solutions >= 2) {
        return;
    }
    const int remaining = static_cast<int>(candidates.size()) - idx;
    if (static_cast<int>(chosen.size()) > need || static_cast<int>(chosen.size()) + remaining < need) {
        return;
    }
    if (idx == static_cast<int>(candidates.size())) {
        State next = base;
        const Guess& g = guesses[static_cast<std::size_t>(guess_idx)];
        std::vector<bool> selected(candidates.size(), false);
        for (int p : chosen) {
            for (std::size_t i = 0; i < candidates.size(); ++i) {
                if (candidates[i] == p) {
                    selected[i] = true;
                    break;
                }
            }
        }

        for (std::size_t i = 0; i < candidates.size(); ++i) {
            const int pos = candidates[i];
            const int digit = g.digits[static_cast<std::size_t>(pos)] - '0';
            const bool must_equal = selected[i];
            const bool ok = must_equal ? apply_equal(next, pos, digit) : apply_not_equal(next, pos, digit);
            if (!ok) {
                return;
            }
        }
        dfs(guesses, next, first_solution, solutions);
        return;
    }

    chosen.push_back(candidates[static_cast<std::size_t>(idx)]);
    branch_guess(guesses, guess_idx, candidates, need, idx + 1, chosen,
                 base, first_solution, solutions);
    chosen.pop_back();

    branch_guess(guesses, guess_idx, candidates, need, idx + 1, chosen,
                 base, first_solution, solutions);
}

std::string solve_puzzle(std::vector<Guess> guesses, bool* unique) {
    if (guesses.empty()) {
        if (unique != nullptr) {
            *unique = false;
        }
        return "";
    }
    const int length = static_cast<int>(guesses.front().digits.size());

    std::sort(guesses.begin(), guesses.end(), [](const Guess& a, const Guess& b) {
        return a.matches < b.matches;
    });

    State initial;
    initial.assign.assign(static_cast<std::size_t>(length), -1);
    initial.banned.assign(static_cast<std::size_t>(length), {});
    for (auto& row : initial.banned) {
        row.fill(false);
    }

    for (const auto& g : guesses) {
        if (static_cast<int>(g.digits.size()) != length) {
            if (unique != nullptr) {
                *unique = false;
            }
            return "";
        }
        if (g.matches == 0) {
            for (int i = 0; i < length; ++i) {
                const int d = g.digits[static_cast<std::size_t>(i)] - '0';
                initial.banned[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)] = true;
            }
        }
    }

    std::string first_solution;
    int solutions = 0;
    dfs(guesses, initial, first_solution, solutions);
    if (unique != nullptr) {
        *unique = (solutions == 1);
    }
    return first_solution;
}

std::vector<Guess> sample_puzzle() {
    return {
        {"90342", 2}, {"70794", 0}, {"39458", 2},
        {"34109", 1}, {"51545", 2}, {"12531", 1},
    };
}

std::vector<Guess> full_puzzle() {
    return {
        {"5616185650518293", 2}, {"3847439647293047", 1},
        {"5855462940810587", 3}, {"9742855507068353", 3},
        {"4296849643607543", 3}, {"3174248439465858", 1},
        {"4513559094146117", 2}, {"7890971548908067", 3},
        {"8157356344118483", 1}, {"2615250744386899", 2},
        {"8690095851526254", 3}, {"6375711915077050", 1},
        {"6913859173121360", 1}, {"6442889055042768", 2},
        {"2321386104303845", 0}, {"2326509471271448", 2},
        {"5251583379644322", 2}, {"1748270476758276", 3},
        {"4895722652190306", 1}, {"3041631117224635", 3},
        {"1841236454324589", 3}, {"2659862637316867", 2},
    };
}

bool run_checkpoints() {
    bool unique = false;
    const std::string sample = solve_puzzle(sample_puzzle(), &unique);
    if (sample != "39542" || !unique) {
        std::cerr << "Checkpoint failed for sample puzzle" << '\n';
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

    bool unique = false;
    const std::string answer = solve_puzzle(full_puzzle(), &unique);
    if (!unique || answer.empty()) {
        std::cerr << "Failed to find a unique solution" << '\n';
        return 3;
    }
    std::cout << answer << '\n';
    return 0;
}
