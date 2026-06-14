#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    std::string file = "resources/documents/1002_input.txt";
    bool run_checkpoints = true;
};

struct Interval {
    int left;
    int right;
};

struct SolveResult {
    bool bipartite;
    int answer;
    i64 crossings;
    int components;
};

bool parse_string_after_prefix(const std::string& arg,
                               const std::string& prefix,
                               std::string& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    value = arg.substr(prefix.size());
    return !value.empty();
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_string_after_prefix(arg, "--file=", options.file)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

class ParityDsu {
  public:
    explicit ParityDsu(const int n) : parent_(n), parity_(n, 0), size_(n, 1), counts_(n, {1, 0}) {
        for (int i = 0; i < n; ++i) {
            parent_[i] = i;
        }
    }

    std::pair<int, int> find(const int x) {
        if (parent_[x] == x) {
            return {x, 0};
        }
        const auto root = find(parent_[x]);
        parity_[x] ^= root.second;
        parent_[x] = root.first;
        return {parent_[x], parity_[x]};
    }

    bool unite(const int a, const int b) {
        const auto [ra, pa] = find(a);
        const auto [rb, pb] = find(b);
        if (ra == rb) {
            return (pa ^ pb) == 1;
        }

        const int link = pa ^ pb ^ 1;
        if (size_[ra] < size_[rb]) {
            parent_[ra] = rb;
            parity_[ra] = link;
            counts_[rb][0] += counts_[ra][link];
            counts_[rb][1] += counts_[ra][link ^ 1];
            size_[rb] += size_[ra];
        } else {
            parent_[rb] = ra;
            parity_[rb] = link;
            counts_[ra][0] += counts_[rb][link];
            counts_[ra][1] += counts_[rb][link ^ 1];
            size_[ra] += size_[rb];
        }
        return true;
    }

    int best_sum() const {
        int total = 0;
        for (std::size_t i = 0; i < parent_.size(); ++i) {
            if (parent_[i] == static_cast<int>(i)) {
                total += std::max(counts_[i][0], counts_[i][1]);
            }
        }
        return total;
    }

    int component_count() const {
        int total = 0;
        for (std::size_t i = 0; i < parent_.size(); ++i) {
            if (parent_[i] == static_cast<int>(i)) {
                ++total;
            }
        }
        return total;
    }

  private:
    std::vector<int> parent_;
    std::vector<int> parity_;
    std::vector<int> size_;
    std::vector<std::array<int, 2>> counts_;
};

std::vector<int> parse_csv_array(const std::string& text) {
    std::vector<int> values;
    std::istringstream input(text);
    std::string token;
    while (std::getline(input, token, ',')) {
        if (token.find_first_not_of(" \t\r\n") == std::string::npos) {
            continue;
        }
        const long long value = std::stoll(token);
        if (value < 0 || value > std::numeric_limits<int>::max()) {
            throw std::runtime_error("Array value is out of range");
        }
        values.push_back(static_cast<int>(value));
    }
    return values;
}

std::vector<int> read_csv_array(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open input file: " + file_path);
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    return parse_csv_array(buffer.str());
}

std::vector<Interval> build_intervals(const std::vector<int>& a) {
    if (a.size() % 2 != 0) {
        throw std::runtime_error("Array length is odd");
    }

    const int n = static_cast<int>(a.size() / 2);
    std::vector<int> first(static_cast<std::size_t>(n), -1);
    std::vector<int> occurrences(static_cast<std::size_t>(n), 0);
    std::vector<Interval> intervals;
    intervals.reserve(static_cast<std::size_t>(n));

    for (int pos = 0; pos < static_cast<int>(a.size()); ++pos) {
        const int value = a[static_cast<std::size_t>(pos)];
        if (value < 0 || value >= n) {
            throw std::runtime_error("Array values must be in [0,n)");
        }
        ++occurrences[static_cast<std::size_t>(value)];
        if (occurrences[static_cast<std::size_t>(value)] == 1) {
            first[static_cast<std::size_t>(value)] = pos;
        } else if (occurrences[static_cast<std::size_t>(value)] == 2) {
            intervals.push_back({first[static_cast<std::size_t>(value)], pos});
        } else {
            throw std::runtime_error("A value occurs more than twice");
        }
    }

    for (const int count : occurrences) {
        if (count != 2) {
            throw std::runtime_error("A value does not occur exactly twice");
        }
    }

    std::sort(intervals.begin(), intervals.end(), [](const Interval& lhs, const Interval& rhs) {
        return lhs.left < rhs.left;
    });
    return intervals;
}

bool crosses(const Interval& lhs, const Interval& rhs) {
    if (lhs.left < rhs.left) {
        return lhs.left < rhs.left && rhs.left < lhs.right && lhs.right < rhs.right;
    }
    return rhs.left < lhs.left && lhs.left < rhs.right && rhs.right < lhs.right;
}

SolveResult solve_array(const std::vector<int>& a) {
    const std::vector<Interval> intervals = build_intervals(a);
    const int n = static_cast<int>(intervals.size());
    ParityDsu dsu(n);
    std::set<std::pair<int, int>> active_by_right;
    i64 crossings = 0;

    for (int i = 0; i < n; ++i) {
        const Interval current = intervals[static_cast<std::size_t>(i)];
        auto it = active_by_right.lower_bound({current.left + 1, -1});
        while (it != active_by_right.end() && it->first < current.right) {
            if (!dsu.unite(i, it->second)) {
                return {false, -1, crossings + 1, 0};
            }
            ++crossings;
            ++it;
        }
        active_by_right.insert({current.right, i});
    }

    return {true, dsu.best_sum(), crossings, dsu.component_count()};
}

SolveResult brute_force_array(const std::vector<int>& a) {
    const std::vector<Interval> intervals = build_intervals(a);
    const int n = static_cast<int>(intervals.size());
    std::vector<std::pair<int, int>> edges;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (crosses(intervals[static_cast<std::size_t>(i)], intervals[static_cast<std::size_t>(j)])) {
                edges.push_back({i, j});
            }
        }
    }

    int best = -1;
    for (int mask = 0; mask < (1 << n); ++mask) {
        bool ok = true;
        for (const auto [u, v] : edges) {
            if (((mask >> u) & 1) == ((mask >> v) & 1)) {
                ok = false;
                break;
            }
        }
        if (ok) {
            best = std::max(best, __builtin_popcount(static_cast<unsigned>(mask)));
        }
    }

    return {best >= 0, best, static_cast<i64>(edges.size()), 0};
}

void validate_all_words(const int n,
                        std::vector<int>& word,
                        std::vector<int>& counts,
                        const int pos) {
    if (pos == 2 * n) {
        const SolveResult fast = solve_array(word);
        const SolveResult brute = brute_force_array(word);
        assert(fast.bipartite == brute.bipartite);
        if (fast.bipartite) {
            assert(fast.answer == brute.answer);
            assert(fast.crossings == brute.crossings);
        }
        return;
    }

    for (int value = 0; value < n; ++value) {
        if (counts[static_cast<std::size_t>(value)] == 2) {
            continue;
        }
        ++counts[static_cast<std::size_t>(value)];
        word[static_cast<std::size_t>(pos)] = value;
        validate_all_words(n, word, counts, pos + 1);
        --counts[static_cast<std::size_t>(value)];
    }
}

void run_checkpoints() {
    assert(solve_array({0, 1, 2, 1, 0, 2}).answer == 2);
    assert(solve_array({0, 0, 1, 1, 2, 2}).answer == 3);
    assert(!solve_array({0, 1, 2, 0, 1, 2}).bipartite);

    for (int n = 1; n <= 4; ++n) {
        std::vector<int> word(static_cast<std::size_t>(2 * n), 0);
        std::vector<int> counts(static_cast<std::size_t>(n), 0);
        validate_all_words(n, word, counts, 0);
    }
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    if (options.run_checkpoints) {
        run_checkpoints();
    }

    try {
        const std::vector<int> input = read_csv_array(options.file);
        const SolveResult result = solve_array(input);
        if (!result.bipartite) {
            throw std::runtime_error("The interval crossing graph is not bipartite");
        }
        std::cout << result.answer << '\n';
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 2;
    }

    return 0;
}
