#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace std;

static constexpr uint64_t MASK = (1ULL << 60) - 1;

struct State {
    uint64_t x;
    uint64_t parent;
};

struct Block {
    uint64_t start;
    int depth;
    bool right;
    uint64_t value;
};

static inline uint64_t combine(uint64_t a, uint64_t b, uint64_t c, uint64_t d) {
    unsigned __int128 v = static_cast<unsigned __int128>(a) * b
                          + static_cast<unsigned __int128>(c) * d;
    return static_cast<uint64_t>(v) & MASK;
}

static State compute_state(uint64_t k) {
    if (k == 0) return {0, 0};
    uint64_t x = 1;
    uint64_t p = 0;
    int msb = 63 - __builtin_clzll(k);
    for (int i = msb - 1; i >= 0; --i) {
        if ((k >> i) & 1ULL) {
            uint64_t nx = combine(2, x, 3, p);
            p = x;
            x = nx;
        } else {
            uint64_t nx = combine(3, x, 2, p);
            p = x;
            x = nx;
        }
    }
    return {x, p};
}

static uint64_t minimax(uint64_t x, uint64_t parent, int depth,
                        uint64_t alpha, uint64_t beta) {
    if (depth == 0) return x;

    uint64_t left_x = combine(3, x, 2, parent);
    uint64_t right_x = combine(2, x, 3, parent);
    bool is_max = (depth & 1) != 0;

    if (is_max) {
        if (left_x < right_x) std::swap(left_x, right_x);
        if (depth == 1) return left_x;
        uint64_t v = minimax(left_x, x, depth - 1, alpha, beta);
        if (v > alpha) alpha = v;
        if (alpha >= beta) return alpha;
        uint64_t v2 = minimax(right_x, x, depth - 1, alpha, beta);
        if (v2 > alpha) alpha = v2;
        return alpha;
    }

    if (left_x > right_x) std::swap(left_x, right_x);
    if (depth == 1) return left_x;
    uint64_t v = minimax(left_x, x, depth - 1, alpha, beta);
    if (v < beta) beta = v;
    if (alpha >= beta) return beta;
    uint64_t v2 = minimax(right_x, x, depth - 1, alpha, beta);
    if (v2 < beta) beta = v2;
    return beta;
}

static uint64_t eval_subtree(uint64_t k, int depth) {
    State s = compute_state(k);
    return minimax(s.x, s.parent, depth, 0, MASK);
}

static void collect_blocks(uint64_t start, int depth, uint64_t left_len,
                           vector<Block>& blocks) {
    uint64_t size = 1ULL << depth;
    if (start + size <= left_len) {
        blocks.push_back({start, depth, false, 0});
        return;
    }
    if (start >= left_len) {
        blocks.push_back({start, depth, true, 0});
        return;
    }
    collect_blocks(start, depth - 1, left_len, blocks);
    collect_blocks(start + (size >> 1), depth - 1, left_len, blocks);
}

static uint64_t compute_block_value(const Block& block, uint64_t base) {
    uint64_t j0 = base + block.start;
    if (!block.right) {
        uint64_t k = j0 >> block.depth;
        return eval_subtree(k, block.depth);
    }
    if (block.depth == 0) {
        uint64_t k = j0 >> 1;
        return MASK - compute_state(k).x;
    }
    // Right blocks are duplicated complements, so one max level collapses.
    uint64_t k = j0 >> block.depth;
    return MASK - eval_subtree(k, block.depth - 1);
}

static uint64_t compute_A(uint64_t n, int threads) {
    if (n == 1) return 1;
    uint64_t total = 2 * n - 1;
    int h = 63 - __builtin_clzll(total);
    uint64_t base = 1ULL << h;
    uint64_t boundary = 2 * n;
    uint64_t left_len = boundary - base;

    // Split the leaf level into maximal power-of-two blocks on each side.
    vector<Block> blocks;
    collect_blocks(0, h, left_len, blocks);

    if (threads < 1) threads = 1;
    if (threads > static_cast<int>(blocks.size())) {
        threads = static_cast<int>(blocks.size());
        if (threads < 1) threads = 1;
    }

    atomic<size_t> index(0);
    vector<thread> pool;
    pool.reserve(threads);
    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&]() {
            size_t i;
            while ((i = index.fetch_add(1)) < blocks.size()) {
                blocks[i].value = compute_block_value(blocks[i], base);
            }
        });
    }
    for (auto& th : pool) th.join();

    unordered_map<uint64_t, uint64_t> values;
    values.reserve(blocks.size() * 2);
    for (const auto& block : blocks) {
        uint64_t key = (block.start << 6) | static_cast<uint64_t>(block.depth);
        values.emplace(key, block.value);
    }

    auto seg = [&](auto&& self, uint64_t start, int depth) -> uint64_t {
        uint64_t key = (start << 6) | static_cast<uint64_t>(depth);
        auto it = values.find(key);
        if (it != values.end()) return it->second;
        uint64_t half = 1ULL << (depth - 1);
        uint64_t left = self(self, start, depth - 1);
        uint64_t right = self(self, start + half, depth - 1);
        return (depth & 1) ? max(left, right) : min(left, right);
    };

    uint64_t z = seg(seg, 0, h);
    return (h & 1) ? (MASK - z) : z;
}

static bool run_validation(int threads) {
    struct Check {
        uint64_t n;
        uint64_t expected;
    };

    const Check checks[] = {
        {4, 8},
        {10, MASK - 33},
        {1000, 101881},
    };

    bool ok = true;
    for (const auto& c : checks) {
        uint64_t got = compute_A(c.n, threads);
        if (got != c.expected) {
            cerr << "Validation failed for n=" << c.n
                 << ": got " << got << ", expected " << c.expected << "\n";
            ok = false;
        }
    }

    if (ok) {
        cerr << "Validation checkpoints passed.\n";
    }
    return ok;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    uint64_t n = 1000000000000ULL;
    unsigned hw = thread::hardware_concurrency();
    int threads = hw ? static_cast<int>(hw) : 1;
    threads = max(1, min(threads, 8));
    bool validate = true;

    // Optional CLI: ./a.out [n] [threads] [validate(0/1)]
    if (argc >= 2) n = stoull(argv[1]);
    if (argc >= 3) threads = max(1, stoi(argv[2]));
    if (argc >= 4) validate = (stoi(argv[3]) != 0);

    if (validate && !run_validation(min(threads, 2))) {
        return 1;
    }

    uint64_t answer = compute_A(n, threads);
    cout << answer << "\n";
    return 0;
}
