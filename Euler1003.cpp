#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using i64 = long long;
using i128 = __int128_t;

constexpr int kTarget = 80;
constexpr int kTaskPrefix = 24;

struct Options {
    int target_k = kTarget;
    bool run_checkpoints = true;
    unsigned threads = 0;
};

struct Power {
    i64 constant = 0;
    i64 coeff = 0;
};

struct Entry {
    i64 coeff = 0;
    i64 constant = 0;
};

struct RightTable {
    std::vector<Entry> entries;
    std::vector<i128> prefix_constant;
};

struct LeftTask {
    int pos = 0;
    int last = -1000000;
    i64 constant = 0;
    i64 coeff = 0;
};

bool entry_less(const Entry& a, const Entry& b) {
    if (a.coeff != b.coeff) {
        return a.coeff < b.coeff;
    }
    return a.constant < b.constant;
}

std::string to_string_i128(i128 value) {
    if (value == 0) {
        return "0";
    }

    bool negative = false;
    if (value < 0) {
        negative = true;
        value = -value;
    }

    std::string result;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        result.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    if (negative) {
        result.push_back('-');
    }
    std::reverse(result.begin(), result.end());
    return result;
}

unsigned resolve_threads(const Options& options) {
    if (options.threads != 0) {
        return std::max(1U, options.threads);
    }
    const unsigned detected = std::thread::hardware_concurrency();
    return detected == 0 ? 1U : detected;
}

std::vector<Power> build_powers(const int k) {
    std::vector<Power> powers(std::max(2, k));
    powers[0] = {1, 0};
    powers[1] = {0, 1};
    for (int i = 2; i < k; ++i) {
        const i64 a = powers[i - 1].constant;
        const i64 b = powers[i - 1].coeff;
        powers[i] = {-2 * b, a - b};
    }
    powers.resize(k);
    return powers;
}

std::pair<i64, i64> evaluate_positions(const std::vector<int>& positions, const std::vector<Power>& powers) {
    i64 constant = 0;
    i64 coeff = 0;
    for (const int pos : positions) {
        constant += powers[pos].constant;
        coeff += powers[pos].coeff;
    }
    return {constant, coeff};
}

std::vector<int> singleton_positions(const i64 n, const int limit) {
    std::vector<i64> stones(static_cast<std::size_t>(limit + 4), 0);
    std::vector<int> positions;
    stones[0] = n;
    for (int i = 0; i < limit; ++i) {
        if ((stones[i] & 1LL) != 0) {
            positions.push_back(i);
        }
        const i64 moved = stones[i] / 2;
        stones[i + 1] += moved;
        stones[i + 3] += moved;
    }
    return positions;
}

int left_category(const int last, const int mid) {
    if (last < 0 || last <= mid - 3) {
        return 0;
    }
    if (last == mid - 2) {
        return 1;
    }
    return 2;
}

int right_category(const int first, const int mid) {
    if (first < 0 || first >= mid + 2) {
        return 0;
    }
    if (first == mid + 1) {
        return 1;
    }
    return 2;
}

bool compatible_boundary(const int left_cat, const int right_cat) {
    if (left_cat == 0) {
        return true;
    }
    if (left_cat == 1) {
        return right_cat != 2;
    }
    return right_cat == 0;
}

void enumerate_right(const int pos,
                     const int end,
                     const int mid,
                     const int last,
                     const int first,
                     const i64 constant,
                     const i64 coeff,
                     const std::vector<Power>& powers,
                     std::array<std::vector<Entry>, 3>& buckets) {
    if (pos == end) {
        buckets[static_cast<std::size_t>(right_category(first, mid))].push_back({coeff, constant});
        return;
    }

    enumerate_right(pos + 1, end, mid, last, first, constant, coeff, powers, buckets);

    if (last < 0 || pos - last >= 3) {
        const int next_first = first < 0 ? pos : first;
        enumerate_right(pos + 1,
                        end,
                        mid,
                        pos,
                        next_first,
                        constant + powers[pos].constant,
                        coeff + powers[pos].coeff,
                        powers,
                        buckets);
    }
}

void build_left_tasks(const int pos,
                      const int stop,
                      const int mid,
                      const int last,
                      const i64 constant,
                      const i64 coeff,
                      const std::vector<Power>& powers,
                      std::vector<LeftTask>& tasks) {
    if (pos == stop) {
        tasks.push_back({pos, last, constant, coeff});
        return;
    }

    build_left_tasks(pos + 1, stop, mid, last, constant, coeff, powers, tasks);

    if (last < 0 || pos - last >= 3) {
        build_left_tasks(pos + 1,
                         stop,
                         mid,
                         pos,
                         constant + powers[pos].constant,
                         coeff + powers[pos].coeff,
                         powers,
                         tasks);
    }
}

void prepare_table(std::vector<Entry>& entries, RightTable& table) {
    std::sort(entries.begin(), entries.end(), entry_less);
    table.entries.swap(entries);
    table.prefix_constant.assign(table.entries.size() + 1, 0);
    for (std::size_t i = 0; i < table.entries.size(); ++i) {
        table.prefix_constant[i + 1] = table.prefix_constant[i] + table.entries[i].constant;
    }
}

i128 query_table(const RightTable& table, const i64 target_coeff, const i64 left_constant) {
    const Entry low{target_coeff, std::numeric_limits<i64>::min()};
    const Entry high{target_coeff, std::numeric_limits<i64>::max()};

    const auto lo = std::lower_bound(table.entries.begin(), table.entries.end(), low, entry_less);
    const auto hi = std::upper_bound(table.entries.begin(), table.entries.end(), high, entry_less);
    if (lo == hi) {
        return 0;
    }

    const Entry threshold{target_coeff, -left_constant};
    const auto first_positive = std::upper_bound(lo, hi, threshold, entry_less);
    const std::size_t first_index = static_cast<std::size_t>(first_positive - table.entries.begin());
    const std::size_t end_index = static_cast<std::size_t>(hi - table.entries.begin());
    const i128 count = static_cast<i128>(end_index - first_index);
    const i128 right_sum = table.prefix_constant[end_index] - table.prefix_constant[first_index];

    return count * static_cast<i128>(left_constant) + right_sum;
}

i128 add_matches(const int left_cat,
                 const i64 left_coeff,
                 const i64 left_constant,
                 const std::array<RightTable, 3>& tables) {
    i128 total = 0;
    const i64 target_coeff = -left_coeff;
    for (int right_cat = 0; right_cat < 3; ++right_cat) {
        if (compatible_boundary(left_cat, right_cat)) {
            total += query_table(tables[static_cast<std::size_t>(right_cat)], target_coeff, left_constant);
        }
    }
    return total;
}

i128 combine_left_dfs(const int pos,
                      const int mid,
                      const int last,
                      const i64 constant,
                      const i64 coeff,
                      const std::vector<Power>& powers,
                      const std::array<RightTable, 3>& tables) {
    if (pos == mid) {
        return add_matches(left_category(last, mid), coeff, constant, tables);
    }

    i128 total = combine_left_dfs(pos + 1, mid, last, constant, coeff, powers, tables);
    if (last < 0 || pos - last >= 3) {
        total += combine_left_dfs(pos + 1,
                                  mid,
                                  pos,
                                  constant + powers[pos].constant,
                                  coeff + powers[pos].coeff,
                                  powers,
                                  tables);
    }
    return total;
}

i128 solve_s(const int k, const unsigned requested_threads) {
    if (k == 0) {
        return 0;
    }

    const int mid = k / 2;
    const std::vector<Power> powers = build_powers(k);

    std::array<std::vector<Entry>, 3> right_buckets;
    enumerate_right(mid, k, mid, -1000000, -1, 0, 0, powers, right_buckets);

    std::array<RightTable, 3> tables;
    for (int i = 0; i < 3; ++i) {
        prepare_table(right_buckets[static_cast<std::size_t>(i)], tables[static_cast<std::size_t>(i)]);
    }

    const int task_stop = std::min(mid, kTaskPrefix);
    std::vector<LeftTask> tasks;
    build_left_tasks(0, task_stop, mid, -1000000, 0, 0, powers, tasks);

    const unsigned workers =
        std::max(1U, std::min<unsigned>(requested_threads, static_cast<unsigned>(std::max<std::size_t>(1, tasks.size()))));

    if (workers == 1) {
        i128 total = 0;
        for (const LeftTask& task : tasks) {
            total += combine_left_dfs(task.pos, mid, task.last, task.constant, task.coeff, powers, tables);
        }
        return total;
    }

    std::vector<i128> partial(workers, 0);
    std::vector<std::thread> pool;
    pool.reserve(workers);
    for (unsigned worker = 0; worker < workers; ++worker) {
        pool.emplace_back([&, worker]() {
            i128 local = 0;
            for (std::size_t i = worker; i < tasks.size(); i += workers) {
                const LeftTask& task = tasks[i];
                local += combine_left_dfs(task.pos, mid, task.last, task.constant, task.coeff, powers, tables);
            }
            partial[worker] = local;
        });
    }
    for (std::thread& thread : pool) {
        thread.join();
    }

    i128 total = 0;
    for (const i128 value : partial) {
        total += value;
    }
    return total;
}

void require_checkpoint(const bool ok, const std::string& message) {
    if (!ok) {
        std::cerr << "Checkpoint failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

void run_checkpoints(const unsigned threads) {
    const std::vector<Power> powers = build_powers(kTarget);

    const std::vector<int> p68{2, 5, 8, 13};
    const std::vector<int> p90{1, 13};
    require_checkpoint(evaluate_positions({0}, powers) == std::make_pair<i64, i64>(1, 0), "n=1 polynomial");
    require_checkpoint(evaluate_positions(p68, powers) == std::make_pair<i64, i64>(68, 0), "n=68 polynomial");
    require_checkpoint(evaluate_positions(p90, powers) == std::make_pair<i64, i64>(90, 0), "n=90 polynomial");

    require_checkpoint(singleton_positions(1, 40) == std::vector<int>{0}, "n=1 trace");
    require_checkpoint(singleton_positions(68, 50) == p68, "n=68 trace");
    require_checkpoint(singleton_positions(90, 50) == p90, "n=90 trace");

    require_checkpoint(solve_s(14, 1) == 159, "S(14)");
    require_checkpoint(solve_s(30, 1) == 33438, "S(30)");
    require_checkpoint(solve_s(30, threads) == 33438, "threaded S(30)");

    std::cerr << "Validation checkpoints passed.\n";
}

void usage() {
    std::cerr << "Usage:\n"
              << "  ./Euler1003 [k] [--skip-checkpoints] [--single-thread] [--threads=N]\n";
}

Options parse_options(const int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
        } else if (arg == "--single-thread") {
            options.threads = 1;
        } else if (arg.rfind("--threads=", 0) == 0) {
            options.threads = static_cast<unsigned>(std::stoul(arg.substr(10)));
            if (options.threads == 0) {
                options.threads = 1;
            }
        } else if (!arg.empty() && arg[0] == '-') {
            usage();
            std::exit(EXIT_FAILURE);
        } else {
            options.target_k = std::stoi(arg);
        }
    }

    if (options.target_k < 0 || options.target_k > kTarget) {
        std::cerr << "k must satisfy 0 <= k <= " << kTarget << ".\n";
        std::exit(EXIT_FAILURE);
    }

    return options;
}

}  // namespace

int main(int argc, char** argv) {
    const Options options = parse_options(argc, argv);
    const unsigned threads = resolve_threads(options);

    if (options.run_checkpoints) {
        run_checkpoints(threads);
    }

    std::cout << to_string_i128(solve_s(options.target_k, threads)) << '\n';
    return 0;
}
