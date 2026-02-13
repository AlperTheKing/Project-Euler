#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }
    std::string s;
    while (x > 0) {
        const int d = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + d));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

struct GameKey {
    std::vector<int> left;
    std::vector<int> right;

    bool operator==(const GameKey& other) const {
        return left == other.left && right == other.right;
    }
};

struct GameKeyHash {
    std::size_t operator()(const GameKey& key) const {
        std::size_t h = 1469598103934665603ULL;
        auto mix = [&](std::size_t x) {
            h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        };

        mix(key.left.size());
        for (int v : key.left) {
            mix(static_cast<std::size_t>(static_cast<u64>(v) + 1ULL));
        }
        mix(0xabcdef01ULL);
        mix(key.right.size());
        for (int v : key.right) {
            mix(static_cast<std::size_t>(static_cast<u64>(v) + 1ULL));
        }
        return h;
    }
};

class GameEngine {
public:
    GameEngine() {
        left_options_.push_back({});
        right_options_.push_back({});
        canonical_.emplace(GameKey{{}, {}}, 0);
    }

    int canonicalize(std::vector<int> left, std::vector<int> right) {
        sort_unique(left);
        sort_unique(right);

        while (true) {
            bool changed = false;

            {
                std::vector<int> kept;
                kept.reserve(left.size());
                for (std::size_t i = 0; i < left.size(); ++i) {
                    bool dominated = false;
                    for (std::size_t j = 0; j < left.size(); ++j) {
                        if (i == j) {
                            continue;
                        }
                        if (leq(left[i], left[j])) {
                            dominated = true;
                            break;
                        }
                    }
                    if (!dominated) {
                        kept.push_back(left[i]);
                    }
                }
                if (kept.size() != left.size()) {
                    left.swap(kept);
                    changed = true;
                }
            }
            if (changed) {
                sort_unique(left);
                continue;
            }

            {
                std::vector<int> kept;
                kept.reserve(right.size());
                for (std::size_t i = 0; i < right.size(); ++i) {
                    bool dominated = false;
                    for (std::size_t j = 0; j < right.size(); ++j) {
                        if (i == j) {
                            continue;
                        }
                        if (leq(right[j], right[i])) {
                            dominated = true;
                            break;
                        }
                    }
                    if (!dominated) {
                        kept.push_back(right[i]);
                    }
                }
                if (kept.size() != right.size()) {
                    right.swap(kept);
                    changed = true;
                }
            }
            if (changed) {
                sort_unique(right);
                continue;
            }

            {
                RawComparator cmp(*this, left, right);
                for (std::size_t i = 0; i < left.size(); ++i) {
                    const int l = left[i];
                    int reversible_via = -1;
                    for (int r : right_options_[static_cast<std::size_t>(l)]) {
                        if (cmp.leq_id_raw(r)) {
                            reversible_via = r;
                            break;
                        }
                    }
                    if (reversible_via >= 0) {
                        std::vector<int> next_left;
                        next_left.reserve(left.size() - 1 +
                                          left_options_[static_cast<std::size_t>(reversible_via)].size());
                        for (std::size_t j = 0; j < left.size(); ++j) {
                            if (j != i) {
                                next_left.push_back(left[j]);
                            }
                        }
                        for (int x : left_options_[static_cast<std::size_t>(reversible_via)]) {
                            next_left.push_back(x);
                        }
                        left.swap(next_left);
                        sort_unique(left);
                        changed = true;
                        break;
                    }
                }
            }
            if (changed) {
                continue;
            }

            {
                RawComparator cmp(*this, left, right);
                for (std::size_t i = 0; i < right.size(); ++i) {
                    const int r = right[i];
                    int reversible_via = -1;
                    for (int l : left_options_[static_cast<std::size_t>(r)]) {
                        if (cmp.leq_raw_id(l)) {
                            reversible_via = l;
                            break;
                        }
                    }
                    if (reversible_via >= 0) {
                        std::vector<int> next_right;
                        next_right.reserve(right.size() - 1 +
                                           right_options_[static_cast<std::size_t>(reversible_via)].size());
                        for (std::size_t j = 0; j < right.size(); ++j) {
                            if (j != i) {
                                next_right.push_back(right[j]);
                            }
                        }
                        for (int x : right_options_[static_cast<std::size_t>(reversible_via)]) {
                            next_right.push_back(x);
                        }
                        right.swap(next_right);
                        sort_unique(right);
                        changed = true;
                        break;
                    }
                }
            }
            if (changed) {
                continue;
            }

            break;
        }

        sort_unique(left);
        sort_unique(right);

        GameKey key{left, right};
        auto it = canonical_.find(key);
        if (it != canonical_.end()) {
            return it->second;
        }

        const int id = static_cast<int>(left_options_.size());
        left_options_.push_back(std::move(left));
        right_options_.push_back(std::move(right));
        canonical_.emplace(GameKey{left_options_.back(), right_options_.back()}, id);
        return id;
    }

    bool leq(int a, int b) {
        const u64 key = pack_ordered(a, b);
        auto it = leq_cache_.find(key);
        if (it != leq_cache_.end()) {
            return it->second != 0;
        }

        bool ok = true;
        for (int l : left_options_[static_cast<std::size_t>(a)]) {
            if (leq(b, l)) {
                ok = false;
                break;
            }
        }
        if (ok) {
            for (int r : right_options_[static_cast<std::size_t>(b)]) {
                if (leq(r, a)) {
                    ok = false;
                    break;
                }
            }
        }

        leq_cache_[key] = static_cast<std::uint8_t>(ok ? 1 : 0);
        return ok;
    }

    int add_games(int a, int b) {
        if (a > b) {
            std::swap(a, b);
        }
        const u64 key = pack_ordered(a, b);
        auto it = add_cache_.find(key);
        if (it != add_cache_.end()) {
            return it->second;
        }

        std::vector<int> left;
        std::vector<int> right;

        left.reserve(left_options_[static_cast<std::size_t>(a)].size() +
                     left_options_[static_cast<std::size_t>(b)].size());
        right.reserve(right_options_[static_cast<std::size_t>(a)].size() +
                      right_options_[static_cast<std::size_t>(b)].size());

        for (int l : left_options_[static_cast<std::size_t>(a)]) {
            left.push_back(add_games(l, b));
        }
        for (int l : left_options_[static_cast<std::size_t>(b)]) {
            left.push_back(add_games(a, l));
        }
        for (int r : right_options_[static_cast<std::size_t>(a)]) {
            right.push_back(add_games(r, b));
        }
        for (int r : right_options_[static_cast<std::size_t>(b)]) {
            right.push_back(add_games(a, r));
        }

        const int sum_id = canonicalize(std::move(left), std::move(right));
        add_cache_[key] = sum_id;
        return sum_id;
    }

    int game_count() const {
        return static_cast<int>(left_options_.size());
    }

private:
    struct RawComparator {
        RawComparator(GameEngine& engine, const std::vector<int>& left, const std::vector<int>& right)
            : engine_(engine), left_(left), right_(right) {}

        bool leq_id_raw(int a) {
            auto it = memo_id_raw_.find(a);
            if (it != memo_id_raw_.end()) {
                return it->second != 0;
            }

            bool ok = true;
            for (int l : engine_.left_options_[static_cast<std::size_t>(a)]) {
                if (leq_raw_id(l)) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                for (int r : right_) {
                    if (engine_.leq(r, a)) {
                        ok = false;
                        break;
                    }
                }
            }

            memo_id_raw_[a] = static_cast<std::uint8_t>(ok ? 1 : 0);
            return ok;
        }

        bool leq_raw_id(int b) {
            auto it = memo_raw_id_.find(b);
            if (it != memo_raw_id_.end()) {
                return it->second != 0;
            }

            bool ok = true;
            for (int l : left_) {
                if (engine_.leq(b, l)) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                for (int r : engine_.right_options_[static_cast<std::size_t>(b)]) {
                    if (leq_id_raw(r)) {
                        ok = false;
                        break;
                    }
                }
            }

            memo_raw_id_[b] = static_cast<std::uint8_t>(ok ? 1 : 0);
            return ok;
        }

        GameEngine& engine_;
        const std::vector<int>& left_;
        const std::vector<int>& right_;
        std::unordered_map<int, std::uint8_t> memo_id_raw_;
        std::unordered_map<int, std::uint8_t> memo_raw_id_;
    };

    static void sort_unique(std::vector<int>& v) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    static u64 pack_ordered(int a, int b) {
        return (static_cast<u64>(static_cast<std::uint32_t>(a)) << 32) |
               static_cast<u64>(static_cast<std::uint32_t>(b));
    }

    std::vector<std::vector<int>> left_options_;
    std::vector<std::vector<int>> right_options_;
    std::unordered_map<GameKey, int, GameKeyHash> canonical_;
    std::unordered_map<u64, std::uint8_t> leq_cache_;
    std::unordered_map<u64, int> add_cache_;
};

struct Children {
    std::array<std::vector<int>, 3> by_digit;
};

std::vector<Children> build_children(int limit) {
    std::vector<Children> children(static_cast<std::size_t>(limit + 1));

    for (int n = 1; n <= limit; ++n) {
        int x = n;
        std::string s;
        while (x > 0) {
            s.push_back(static_cast<char>('0' + (x % 3)));
            x /= 3;
        }
        std::reverse(s.begin(), s.end());

        for (std::size_t i = 0; i < s.size(); ++i) {
            std::string t = s.substr(0, i) + s.substr(i + 1);
            const int m = t.empty() ? 0 : std::stoi(t, nullptr, 3);
            const int d = s[i] - '0';
            children[static_cast<std::size_t>(n)].by_digit[static_cast<std::size_t>(d)].push_back(m);
        }
    }

    return children;
}

size_t triangular_index(size_t i, size_t j, size_t m) {
    return i * m - (i * (i - 1)) / 2 + (j - i);
}

u128 solve_fast(int limit, unsigned thread_count) {
    const auto children = build_children(limit);

    GameEngine engine;
    std::vector<int> gid(static_cast<std::size_t>(limit + 1), 0);

    for (int n = 1; n <= limit; ++n) {
        std::vector<int> left;
        std::vector<int> right;

        for (int m : children[static_cast<std::size_t>(n)].by_digit[0]) {
            left.push_back(gid[static_cast<std::size_t>(m)]);
        }
        for (int m : children[static_cast<std::size_t>(n)].by_digit[2]) {
            left.push_back(gid[static_cast<std::size_t>(m)]);
        }

        for (int m : children[static_cast<std::size_t>(n)].by_digit[1]) {
            right.push_back(gid[static_cast<std::size_t>(m)]);
        }
        for (int m : children[static_cast<std::size_t>(n)].by_digit[2]) {
            right.push_back(gid[static_cast<std::size_t>(m)]);
        }

        gid[static_cast<std::size_t>(n)] = engine.canonicalize(std::move(left), std::move(right));
    }

    std::unordered_map<int, u64> freq_map;
    freq_map.reserve(static_cast<std::size_t>(limit));
    for (int n = 1; n <= limit; ++n) {
        ++freq_map[gid[static_cast<std::size_t>(n)]];
    }

    std::vector<int> ids;
    ids.reserve(freq_map.size());
    for (const auto& [id, _] : freq_map) {
        ids.push_back(id);
    }
    std::sort(ids.begin(), ids.end());

    std::vector<u64> freq;
    freq.reserve(ids.size());
    for (int id : ids) {
        freq.push_back(freq_map[id]);
    }

    const size_t m = ids.size();
    const size_t tri_size = m * (m + 1) / 2;
    std::vector<int> sum_ids(tri_size, 0);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = i; j < m; ++j) {
            sum_ids[triangular_index(i, j, m)] =
                engine.add_games(ids[i], ids[j]);
        }
    }

    const size_t game_count = static_cast<size_t>(engine.game_count());

    if (thread_count == 0) {
        thread_count = 1;
    }
    thread_count = std::min<unsigned>(thread_count, static_cast<unsigned>(m == 0 ? 1 : m));

    std::vector<std::vector<u64>> locals(static_cast<size_t>(thread_count),
                                         std::vector<u64>(game_count, 0ULL));
    std::vector<std::thread> threads;
    threads.reserve(thread_count);

    const size_t chunk = (m + static_cast<size_t>(thread_count) - 1) / static_cast<size_t>(thread_count);

    for (unsigned t = 0; t < thread_count; ++t) {
        const size_t begin = static_cast<size_t>(t) * chunk;
        const size_t end = std::min(m, begin + chunk);

        threads.emplace_back([&, t, begin, end]() {
            auto& local = locals[static_cast<size_t>(t)];
            for (size_t i = begin; i < end; ++i) {
                const u64 fi = freq[i];
                local[static_cast<size_t>(sum_ids[triangular_index(i, i, m)])] += fi * (fi + 1ULL) / 2ULL;
                for (size_t j = i + 1; j < m; ++j) {
                    local[static_cast<size_t>(sum_ids[triangular_index(i, j, m)])] += fi * freq[j];
                }
            }
        });
    }

    for (auto& th : threads) {
        th.join();
    }

    std::vector<u64> pair_count(game_count, 0ULL);
    for (const auto& local : locals) {
        for (size_t i = 0; i < game_count; ++i) {
            pair_count[i] += local[i];
        }
    }

    u128 ans = 0;
    for (u64 c : pair_count) {
        ans += static_cast<u128>(c) * static_cast<u128>(c);
    }

    return ans;
}

class BruteSolver {
public:
    explicit BruteSolver(int limit) : limit_(limit), children_(build_children(limit)) {}

    bool is_fair_state(int a, int b, int c, int d) {
        const auto [sa, sb] = ordered(a, b);
        const auto [sc, sd] = ordered(c, d);
        return !can_win(sa, sb, sc, sd, 0) && !can_win(sa, sb, sc, sd, 1);
    }

    u64 F() {
        u64 count = 0;
        for (int a = 1; a <= limit_; ++a) {
            for (int b = a; b <= limit_; ++b) {
                for (int c = 1; c <= limit_; ++c) {
                    for (int d = c; d <= limit_; ++d) {
                        if (is_fair_state(a, b, c, d)) {
                            ++count;
                        }
                    }
                }
            }
        }
        return count;
    }

private:
    static std::pair<int, int> ordered(int x, int y) {
        if (x <= y) {
            return {x, y};
        }
        return {y, x};
    }

    static u64 pack_state(int a, int b, int c, int d, int turn) {
        u64 key = static_cast<u64>(a);
        key = (key << 17) | static_cast<u64>(b);
        key = (key << 17) | static_cast<u64>(c);
        key = (key << 17) | static_cast<u64>(d);
        key = (key << 1) | static_cast<u64>(turn);
        return key;
    }

    bool can_win(int a, int b, int c, int d, int turn) {
        const u64 key = pack_state(a, b, c, d, turn);
        auto it = memo_.find(key);
        if (it != memo_.end()) {
            return it->second;
        }

        bool win = false;

        auto try_next = [&](int na, int nb, int nc, int nd, int nt) {
            const auto [sa, sb] = ordered(na, nb);
            const auto [sc, sd] = ordered(nc, nd);
            if (!can_win(sa, sb, sc, sd, nt)) {
                win = true;
            }
        };

        if (turn == 0) {
            for (int y : children_[static_cast<size_t>(a)].by_digit[0]) {
                try_next(y, b, c, d, 1);
                if (win) break;
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(b)].by_digit[0]) {
                    try_next(a, y, c, d, 1);
                    if (win) break;
                }
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(c)].by_digit[1]) {
                    try_next(a, b, y, d, 1);
                    if (win) break;
                }
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(d)].by_digit[1]) {
                    try_next(a, b, c, y, 1);
                    if (win) break;
                }
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(a)].by_digit[2]) {
                    try_next(y, b, c, d, 1);
                    if (win) break;
                }
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(b)].by_digit[2]) {
                    try_next(a, y, c, d, 1);
                    if (win) break;
                }
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(c)].by_digit[2]) {
                    try_next(a, b, y, d, 1);
                    if (win) break;
                }
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(d)].by_digit[2]) {
                    try_next(a, b, c, y, 1);
                    if (win) break;
                }
            }
        } else {
            for (int y : children_[static_cast<size_t>(c)].by_digit[0]) {
                try_next(a, b, y, d, 0);
                if (win) break;
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(d)].by_digit[0]) {
                    try_next(a, b, c, y, 0);
                    if (win) break;
                }
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(a)].by_digit[1]) {
                    try_next(y, b, c, d, 0);
                    if (win) break;
                }
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(b)].by_digit[1]) {
                    try_next(a, y, c, d, 0);
                    if (win) break;
                }
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(a)].by_digit[2]) {
                    try_next(y, b, c, d, 0);
                    if (win) break;
                }
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(b)].by_digit[2]) {
                    try_next(a, y, c, d, 0);
                    if (win) break;
                }
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(c)].by_digit[2]) {
                    try_next(a, b, y, d, 0);
                    if (win) break;
                }
            }
            if (!win) {
                for (int y : children_[static_cast<size_t>(d)].by_digit[2]) {
                    try_next(a, b, c, y, 0);
                    if (win) break;
                }
            }
        }

        memo_[key] = win;
        return win;
    }

    int limit_;
    std::vector<Children> children_;
    std::unordered_map<u64, bool> memo_;
};

void run_checkpoints() {
    BruteSolver brute(5);
    assert(brute.is_fair_state(1, 5, 2, 4));
    assert(brute.F() == 21ULL);

    const u128 fast_small = solve_fast(5, 1);
    assert(fast_small == 21ULL);
}

}  // namespace

int main() {
    run_checkpoints();
    const unsigned threads = std::max(1U, std::thread::hardware_concurrency());
    const u128 answer = solve_fast(100'000, threads);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
