#include <atomic>
#include <cstdint>
#include <deque>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kThresholdNumerator = 1;     // x < 1/100
constexpr u64 kThresholdDenominator = 100; // x < 1/100
constexpr u64 kTargetQ = 100'000'000ULL;

struct Node {
    u64 a;
    u64 b;
    u64 c;
    u64 d;
};

inline bool exceeds_product_bound(const Node& node, u64 max_product) {
    return node.b > max_product / node.d;
}

inline bool left_endpoint_not_below_threshold(const Node& node) {
    return static_cast<u128>(node.a) * kThresholdDenominator >=
           static_cast<u128>(kThresholdNumerator) * node.b;
}

inline bool midpoint_below_threshold(const Node& node) {
    // (a/b + c/d) / 2 < T  <=>  (a*d + b*c) * T_den < 2*T_num*b*d
    const u128 lhs = (static_cast<u128>(node.a) * node.d +
                      static_cast<u128>(node.b) * node.c) *
                     kThresholdDenominator;
    const u128 rhs = static_cast<u128>(2) * kThresholdNumerator * node.b * node.d;
    return lhs < rhs;
}

u64 count_subtree_iterative(const Node& seed, u64 max_product) {
    std::vector<Node> stack;
    stack.reserve(4096);
    stack.push_back(seed);

    u64 count = 0;

    while (!stack.empty()) {
        const Node node = stack.back();
        stack.pop_back();

        if (exceeds_product_bound(node, max_product)) {
            continue;
        }

        if (left_endpoint_not_below_threshold(node)) {
            continue;
        }

        if (midpoint_below_threshold(node)) {
            ++count;
        }

        const u64 mediant_n = node.a + node.c;
        const u64 mediant_d = node.b + node.d;

        // Push left first, then right; right is processed first with LIFO.
        // This keeps memory low on long left spines like (0/1, 1/n).
        stack.push_back(Node{node.a, node.b, mediant_n, mediant_d});
        stack.push_back(Node{mediant_n, mediant_d, node.c, node.d});
    }

    return count;
}

struct SeedFrontier {
    std::vector<Node> seeds;
    u64 expanded_node_count = 0;
};

SeedFrontier build_seed_frontier_with_prefix_count(u64 max_product,
                                                   std::size_t desired_seeds) {
    std::deque<Node> queue;
    queue.push_back(Node{0, 1, 1, 1});
    u64 expanded_node_count = 0;

    while (!queue.empty() && queue.size() < desired_seeds) {
        const Node node = queue.front();
        queue.pop_front();

        if (exceeds_product_bound(node, max_product) ||
            left_endpoint_not_below_threshold(node)) {
            continue;
        }

        if (midpoint_below_threshold(node)) {
            ++expanded_node_count;
        }

        const u64 mediant_n = node.a + node.c;
        const u64 mediant_d = node.b + node.d;

        queue.push_back(Node{node.a, node.b, mediant_n, mediant_d});
        queue.push_back(Node{mediant_n, mediant_d, node.c, node.d});
    }

    return SeedFrontier{
        std::vector<Node>(queue.begin(), queue.end()),
        expanded_node_count,
    };
}

u64 count_ambiguous_fast(u64 q_limit, bool allow_multithreading, unsigned requested_threads = 0) {
    if (q_limit < 2) {
        return 0;
    }

    const u64 max_product = q_limit / 2;

    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 1;
        }
    }

    if (!allow_multithreading || threads <= 1 || q_limit < 1'000'000ULL) {
        return count_subtree_iterative(Node{0, 1, 1, 1}, max_product);
    }

    const std::size_t desired_seeds = static_cast<std::size_t>(threads) * 64;
    const SeedFrontier frontier =
        build_seed_frontier_with_prefix_count(max_product, desired_seeds);
    const std::vector<Node>& seeds = frontier.seeds;

    if (seeds.empty()) {
        return frontier.expanded_node_count;
    }

    if (seeds.size() == 1) {
        return frontier.expanded_node_count +
               count_subtree_iterative(seeds[0], max_product);
    }

    threads = static_cast<unsigned>(std::min<std::size_t>(threads, seeds.size()));

    std::atomic<std::size_t> next_index{0};
    std::vector<u64> partial(threads, 0);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (unsigned tid = 0; tid < threads; ++tid) {
        workers.emplace_back([&, tid]() {
            u64 local = 0;
            while (true) {
                const std::size_t index = next_index.fetch_add(1, std::memory_order_relaxed);
                if (index >= seeds.size()) {
                    break;
                }
                local += count_subtree_iterative(seeds[index], max_product);
            }
            partial[tid] = local;
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    u64 total = 0;
    for (u64 value : partial) {
        total += value;
    }
    return frontier.expanded_node_count + total;
}

u64 count_ambiguous_bruteforce(u64 q_limit) {
    struct Frac {
        u64 p;
        u64 q;
    };

    std::vector<Frac> targets;
    for (u64 q = 1; q <= q_limit; ++q) {
        for (u64 p = 1; p < q; ++p) {
            if (std::gcd(p, q) != 1) {
                continue;
            }
            if (p * kThresholdDenominator >= kThresholdNumerator * q) {
                continue;
            }
            targets.push_back(Frac{p, q});
        }
    }

    std::vector<Frac> approximants;
    std::vector<std::vector<Frac>> approximants_by_bound(q_limit + 1);

    for (u64 den = 1; den <= q_limit; ++den) {
        for (u64 num = 0; num <= den; ++num) {
            if (std::gcd(num, den) == 1) {
                approximants.push_back(Frac{num, den});
            }
        }
        approximants_by_bound[den] = approximants;
    }

    u64 ambiguous_count = 0;

    for (const Frac target : targets) {
        bool ambiguous = false;

        for (u64 bound = 1; bound <= q_limit && !ambiguous; ++bound) {
            u64 best_num = 0;
            u64 best_den = 1;
            bool has_best = false;
            int best_count = 0;

            for (const Frac cand : approximants_by_bound[bound]) {
                const std::int64_t diff = static_cast<std::int64_t>(target.p * cand.q) -
                                          static_cast<std::int64_t>(cand.p * target.q);
                const u64 num = static_cast<u64>(diff < 0 ? -diff : diff);
                const u64 den = target.q * cand.q;

                if (!has_best) {
                    has_best = true;
                    best_num = num;
                    best_den = den;
                    best_count = 1;
                    continue;
                }

                const u128 lhs = static_cast<u128>(num) * best_den;
                const u128 rhs = static_cast<u128>(best_num) * den;

                if (lhs < rhs) {
                    best_num = num;
                    best_den = den;
                    best_count = 1;
                } else if (lhs == rhs) {
                    ++best_count;
                }
            }

            if (best_count >= 2) {
                ambiguous = true;
            }
        }

        if (ambiguous) {
            ++ambiguous_count;
        }
    }

    return ambiguous_count;
}

bool run_validation_checkpoints() {
    struct Checkpoint {
        u64 q_limit;
        u64 expected;
    };

    const std::vector<Checkpoint> checkpoints = {
        {120, 10},
        {200'000, 100'967},
        {1'000'000, 509'763},
    };

    for (const Checkpoint cp : checkpoints) {
        const u64 got = count_ambiguous_fast(cp.q_limit, false, 1);
        if (got != cp.expected) {
            std::cerr << "Checkpoint failed for q <= " << cp.q_limit << ": got " << got
                      << ", expected " << cp.expected << "\n";
            return false;
        }
    }

    const u64 brute_q_limit = 120;
    const u64 brute_expected = checkpoints.front().expected;
    const u64 brute_got = count_ambiguous_bruteforce(brute_q_limit);

    if (brute_got != brute_expected) {
        std::cerr << "Brute-force checkpoint failed for q <= " << brute_q_limit << ": got "
                  << brute_got << ", expected " << brute_expected << "\n";
        return false;
    }

    return true;
}

} // namespace

int main() {
    if (!run_validation_checkpoints()) {
        return 1;
    }

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 4;
    }

    const u64 answer = count_ambiguous_fast(kTargetQ, true, threads);
    std::cout << answer << '\n';
    return 0;
}
