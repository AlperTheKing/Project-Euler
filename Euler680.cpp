#include <algorithm>
#include <cstdint>
#include <exception>
#include <future>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    u64 n = 1000000000000000000ULL;
    u64 k = 1000000ULL;
    u64 mod = 1000000000ULL;
    bool run_checkpoints = true;
    unsigned requested_threads = 0U;
};

bool parse_u64_after_prefix(const std::string& arg,
                            const std::string& prefix,
                            u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }

        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            throw std::overflow_error("u64 argument overflow");
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const std::string& prefix,
                                 unsigned& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }

        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<unsigned>::max() - digit) / 10ULL) {
            throw std::overflow_error("unsigned argument overflow");
        }
        parsed = parsed * 10ULL + digit;
    }

    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--N=", options.n)) {
            continue;
        }
        if (parse_u64_after_prefix(arg, "--K=", options.k)) {
            continue;
        }
        if (parse_u64_after_prefix(arg, "--mod=", options.mod)) {
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", options.requested_threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

unsigned pick_thread_count(const unsigned requested) {
    if (requested > 0U) {
        return requested;
    }
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0U) {
        hw = 4U;
    }
    return hw;
}

u64 add_mod_u64(const u64 a, const u64 b, const u64 mod) {
    if (mod == 1ULL) {
        return 0ULL;
    }
    const u64 x = a % mod;
    const u64 y = b % mod;
    const u64 threshold = mod - y;
    if (x >= threshold) {
        return x - threshold;
    }
    return x + y;
}

u64 sub_mod_u64(const u64 a, const u64 b, const u64 mod) {
    if (mod == 1ULL) {
        return 0ULL;
    }
    const u64 x = a % mod;
    const u64 y = b % mod;
    if (x >= y) {
        return x - y;
    }
    return mod - (y - x);
}

u64 mul_mod_u64(const u64 a, const u64 b, const u64 mod) {
    if (mod == 1ULL) {
        return 0ULL;
    }
    return static_cast<u64>((static_cast<u128>(a % mod) * static_cast<u128>(b % mod)) %
                            static_cast<u128>(mod));
}

u64 fibonacci_add_mod_n(const u64 a, const u64 b, const u64 n) {
    if (n == 0ULL) {
        return 0ULL;
    }
    return static_cast<u64>((static_cast<u128>(a) + static_cast<u128>(b)) %
                            static_cast<u128>(n));
}

class ReversalTreapSolver {
public:
    ReversalTreapSolver(const u64 n, const u64 k, const u64 mod)
        : n_(n), k_(k), mod_(mod), rng_state_(0x9e3779b97f4a7c15ULL), root_(nullptr) {
        if (mod_ == 0ULL) {
            throw std::invalid_argument("mod must be positive");
        }
        if (n_ == 0ULL) {
            throw std::invalid_argument("N must be positive");
        }

        const u128 required_nodes = static_cast<u128>(2) * static_cast<u128>(k_) + 8U;
        if (required_nodes > static_cast<u128>(std::numeric_limits<std::size_t>::max())) {
            throw std::overflow_error("node pool size exceeds size_t");
        }
        nodes_.reserve(static_cast<std::size_t>(required_nodes));

        root_ = make_node(n_, 0ULL, +1);
    }

    u64 solve() {
        u64 s = 1ULL % n_;
        u64 t = 1ULL % n_;

        for (u64 step = 0; step < k_; ++step) {
            const u64 left = std::min(s, t);
            const u64 right = std::max(s, t);
            reverse_interval(left, right);

            const u64 next_s = fibonacci_add_mod_n(s, t, n_);
            const u64 next_t = fibonacci_add_mod_n(t, next_s, n_);
            s = next_s;
            t = next_t;
        }

        return root_->weighted_sum % mod_;
    }

private:
    struct Node {
        Node* left = nullptr;
        Node* right = nullptr;

        u64 len_segment = 0ULL;
        u64 first_value = 0ULL;
        int diff = +1;

        u64 total_len = 0ULL;
        u64 segment_sum = 0ULL;
        u64 segment_weighted_sum = 0ULL;
        u64 sum = 0ULL;
        u64 weighted_sum = 0ULL;

        std::uint32_t priority = 0U;
        bool rev = false;
    };

    const u64 n_;
    const u64 k_;
    const u64 mod_;
    u64 rng_state_;
    std::vector<Node> nodes_;
    Node* root_;

    u64 next_random() {
        u64 z = (rng_state_ += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30U)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27U)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31U);
    }

    u64 triangular_mod(const u64 n) const {
        if (mod_ == 1ULL) {
            return 0ULL;
        }
        u128 a = static_cast<u128>(n);
        u128 b = static_cast<u128>(n - 1ULL);
        if ((a & 1U) == 0U) {
            a /= 2U;
        } else {
            b /= 2U;
        }

        const u64 am = static_cast<u64>(a % mod_);
        const u64 bm = static_cast<u64>(b % mod_);
        return mul_mod_u64(am, bm, mod_);
    }

    u64 square_sum_mod(const u64 n) const {
        if (mod_ == 1ULL) {
            return 0ULL;
        }
        u128 a = static_cast<u128>(n);
        u128 b = static_cast<u128>(n - 1ULL);
        u128 c = static_cast<u128>(2ULL) * static_cast<u128>(n) - 1U;

        auto divide_factor = [&](const u64 factor) {
            if (a % factor == 0U) {
                a /= factor;
                return;
            }
            if (b % factor == 0U) {
                b /= factor;
                return;
            }
            c /= factor;
        };

        divide_factor(2ULL);
        divide_factor(3ULL);

        const u64 am = static_cast<u64>(a % mod_);
        const u64 bm = static_cast<u64>(b % mod_);
        const u64 cm = static_cast<u64>(c % mod_);
        return mul_mod_u64(mul_mod_u64(am, bm, mod_), cm, mod_);
    }

    void recompute_segment(Node* node) {
        const u64 len = node->len_segment;
        const u64 a = node->first_value % mod_;
        const u64 tri = triangular_mod(len);
        const u64 sq = square_sum_mod(len);

        const u64 linear_sum = mul_mod_u64(len % mod_, a, mod_);
        if (node->diff > 0) {
            node->segment_sum = add_mod_u64(linear_sum, tri, mod_);
        } else {
            node->segment_sum = sub_mod_u64(linear_sum, tri, mod_);
        }

        const u64 weighted_linear = mul_mod_u64(a, tri, mod_);
        if (node->diff > 0) {
            node->segment_weighted_sum = add_mod_u64(weighted_linear, sq, mod_);
        } else {
            node->segment_weighted_sum = sub_mod_u64(weighted_linear, sq, mod_);
        }
    }

    Node* make_node(const u64 len_segment, const u64 first_value, const int diff) {
        if (len_segment == 0ULL) {
            return nullptr;
        }
        if (nodes_.size() >= nodes_.capacity()) {
            throw std::runtime_error("Treap node pool exhausted");
        }

        nodes_.push_back(Node{});
        Node* node = &nodes_.back();
        node->left = nullptr;
        node->right = nullptr;
        node->len_segment = len_segment;
        node->first_value = first_value;
        node->diff = diff;
        node->priority = static_cast<std::uint32_t>(next_random());
        node->rev = false;
        recompute_segment(node);
        node->total_len = len_segment;
        node->sum = node->segment_sum;
        node->weighted_sum = node->segment_weighted_sum;
        return node;
    }

    static u64 node_len(const Node* node) {
        return node == nullptr ? 0ULL : node->total_len;
    }

    static u64 node_sum(const Node* node) {
        return node == nullptr ? 0ULL : node->sum;
    }

    static u64 node_weighted_sum(const Node* node) {
        return node == nullptr ? 0ULL : node->weighted_sum;
    }

    void pull(Node* node) {
        if (node == nullptr) {
            return;
        }

        const u64 left_len = node_len(node->left);
        const u64 right_len = node_len(node->right);
        node->total_len = left_len + node->len_segment + right_len;

        const u64 left_sum = node_sum(node->left);
        const u64 right_sum = node_sum(node->right);
        node->sum =
            add_mod_u64(add_mod_u64(left_sum, node->segment_sum, mod_), right_sum, mod_);

        const u64 left_weighted = node_weighted_sum(node->left);
        const u64 right_weighted = node_weighted_sum(node->right);

        u64 total_weighted = left_weighted;

        const u64 left_shift = left_len % mod_;
        const u64 segment_shifted = add_mod_u64(
            node->segment_weighted_sum,
            mul_mod_u64(left_shift, node->segment_sum, mod_),
            mod_);
        total_weighted = add_mod_u64(total_weighted, segment_shifted, mod_);

        const u64 right_shift = add_mod_u64(left_shift, node->len_segment % mod_, mod_);
        const u64 right_shifted =
            add_mod_u64(right_weighted, mul_mod_u64(right_shift, right_sum, mod_), mod_);
        total_weighted = add_mod_u64(total_weighted, right_shifted, mod_);

        node->weighted_sum = total_weighted;
    }

    void flip_segment(Node* node) {
        if (node->diff > 0) {
            node->first_value += node->len_segment - 1ULL;
        } else {
            node->first_value -= node->len_segment - 1ULL;
        }
        node->diff = -node->diff;
        recompute_segment(node);
    }

    void mark_reversed(Node* node) {
        if (node == nullptr) {
            return;
        }

        node->rev = !node->rev;
        std::swap(node->left, node->right);
        flip_segment(node);

        const u64 factor = (node->total_len - 1ULL) % mod_;
        const u64 reflected = mul_mod_u64(factor, node->sum, mod_);
        node->weighted_sum = sub_mod_u64(reflected, node->weighted_sum, mod_);
    }

    void push(Node* node) {
        if (node == nullptr || !node->rev) {
            return;
        }
        mark_reversed(node->left);
        mark_reversed(node->right);
        node->rev = false;
    }

    Node* merge(Node* left, Node* right) {
        if (left == nullptr) {
            return right;
        }
        if (right == nullptr) {
            return left;
        }

        if (left->priority > right->priority) {
            push(left);
            left->right = merge(left->right, right);
            pull(left);
            return left;
        }

        push(right);
        right->left = merge(left, right->left);
        pull(right);
        return right;
    }

    void split(Node* node, const u64 take_left, Node*& left, Node*& right) {
        if (node == nullptr) {
            left = nullptr;
            right = nullptr;
            return;
        }

        push(node);
        const u64 left_len = node_len(node->left);

        if (take_left < left_len) {
            split(node->left, take_left, left, node->left);
            pull(node);
            right = node;
            return;
        }

        const u64 left_and_segment = left_len + node->len_segment;
        if (take_left > left_and_segment) {
            split(node->right, take_left - left_and_segment, node->right, right);
            pull(node);
            left = node;
            return;
        }

        if (take_left == left_len) {
            left = node->left;
            node->left = nullptr;
            pull(node);
            right = node;
            return;
        }

        if (take_left == left_and_segment) {
            right = node->right;
            node->right = nullptr;
            pull(node);
            left = node;
            return;
        }

        const u64 inside = take_left - left_len;
        Node* left_sub = node->left;
        Node* right_sub = node->right;

        const u64 old_len = node->len_segment;
        const u64 old_first = node->first_value;
        const int old_diff = node->diff;

        node->left = nullptr;
        node->right = nullptr;
        node->len_segment = inside;
        node->first_value = old_first;
        node->diff = old_diff;
        node->rev = false;
        recompute_segment(node);
        pull(node);

        const u64 right_first = old_diff > 0 ? old_first + inside : old_first - inside;
        Node* right_piece = make_node(old_len - inside, right_first, old_diff);

        left = merge(left_sub, node);
        right = merge(right_piece, right_sub);
    }

    void reverse_interval(const u64 left_idx, const u64 right_idx) {
        Node* left = nullptr;
        Node* middle = nullptr;
        Node* right = nullptr;

        split(root_, left_idx, left, middle);
        split(middle, right_idx - left_idx + 1ULL, middle, right);
        mark_reversed(middle);
        root_ = merge(left, merge(middle, right));
    }
};

u128 brute_force_r_exact(const u64 n, const u64 k) {
    std::vector<u64> array(static_cast<std::size_t>(n));
    std::iota(array.begin(), array.end(), 0ULL);

    u64 s = 1ULL % n;
    u64 t = 1ULL % n;

    for (u64 step = 0; step < k; ++step) {
        const u64 left = std::min(s, t);
        const u64 right = std::max(s, t);
        std::reverse(array.begin() + static_cast<std::ptrdiff_t>(left),
                     array.begin() + static_cast<std::ptrdiff_t>(right + 1ULL));

        const u64 next_s = fibonacci_add_mod_n(s, t, n);
        const u64 next_t = fibonacci_add_mod_n(t, next_s, n);
        s = next_s;
        t = next_t;
    }

    u128 result = 0;
    for (u64 i = 0; i < n; ++i) {
        result += static_cast<u128>(i) * static_cast<u128>(array[static_cast<std::size_t>(i)]);
    }
    return result;
}

void run_random_consistency_checks(const unsigned threads, const u64 mod) {
    const unsigned worker_count = std::max(1U, std::min(threads, 8U));
    const int total_cases = 32;

    std::vector<std::future<void>> futures;
    futures.reserve(worker_count);

    for (unsigned worker = 0U; worker < worker_count; ++worker) {
        futures.emplace_back(std::async(std::launch::async, [worker, worker_count, mod]() {
            std::mt19937_64 rng(0x123456789abcdef0ULL + worker * 0x9e3779b97f4a7c15ULL);
            std::uniform_int_distribution<int> n_dist(1, 120);
            std::uniform_int_distribution<int> k_dist(1, 140);

            for (int case_idx = static_cast<int>(worker); case_idx < total_cases;
                 case_idx += static_cast<int>(worker_count)) {
                const u64 n = static_cast<u64>(n_dist(rng));
                const u64 k = static_cast<u64>(k_dist(rng));

                const u128 brute = brute_force_r_exact(n, k);
                ReversalTreapSolver fast_solver(n, k, mod);
                const u64 fast = fast_solver.solve();
                const u64 expected = static_cast<u64>(brute % static_cast<u128>(mod));

                if (fast != expected) {
                    throw std::runtime_error("Random consistency check failed for N=" +
                                             std::to_string(n) + ", K=" + std::to_string(k));
                }
            }
        }));
    }

    for (auto& future : futures) {
        future.get();
    }
}

void run_checkpoints(const unsigned threads, const u64 mod) {
    auto cp1 = std::async(std::launch::async, []() { return brute_force_r_exact(5ULL, 4ULL); });
    auto cp2 =
        std::async(std::launch::async, []() { return brute_force_r_exact(100ULL, 100ULL); });
    auto cp3 = std::async(std::launch::async,
                          []() { return brute_force_r_exact(10000ULL, 10000ULL); });

    const u128 v1 = cp1.get();
    const u128 v2 = cp2.get();
    const u128 v3 = cp3.get();

    if (v1 != static_cast<u128>(27ULL)) {
        throw std::runtime_error("Checkpoint failed: R(5,4)");
    }
    if (v2 != static_cast<u128>(246597ULL)) {
        throw std::runtime_error("Checkpoint failed: R(100,100)");
    }
    if (v3 != static_cast<u128>(249275481640ULL)) {
        throw std::runtime_error("Checkpoint failed: R(10000,10000)");
    }

    auto fast1 = std::async(std::launch::async, [mod]() {
        ReversalTreapSolver solver(5ULL, 4ULL, mod);
        return solver.solve();
    });
    auto fast2 = std::async(std::launch::async, [mod]() {
        ReversalTreapSolver solver(100ULL, 100ULL, mod);
        return solver.solve();
    });
    auto fast3 = std::async(std::launch::async, [mod]() {
        ReversalTreapSolver solver(10000ULL, 10000ULL, mod);
        return solver.solve();
    });

    if (fast1.get() != static_cast<u64>(v1 % static_cast<u128>(mod))) {
        throw std::runtime_error("Fast solver mismatch: R(5,4)");
    }
    if (fast2.get() != static_cast<u64>(v2 % static_cast<u128>(mod))) {
        throw std::runtime_error("Fast solver mismatch: R(100,100)");
    }
    if (fast3.get() != static_cast<u64>(v3 % static_cast<u128>(mod))) {
        throw std::runtime_error("Fast solver mismatch: R(10000,10000)");
    }

    run_random_consistency_checks(threads, mod);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        Options options;
        if (!parse_arguments(argc, argv, options)) {
            return 1;
        }
        if (options.n == 0ULL) {
            throw std::invalid_argument("N must be positive");
        }
        if (options.mod == 0ULL) {
            throw std::invalid_argument("mod must be positive");
        }

        const unsigned threads = pick_thread_count(options.requested_threads);
        if (options.run_checkpoints) {
            run_checkpoints(threads, options.mod);
        }

        ReversalTreapSolver solver(options.n, options.k, options.mod);
        const u64 answer = solver.solve();
        std::cout << answer << '\n';
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
}
