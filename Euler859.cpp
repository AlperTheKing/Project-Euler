#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using u64 = std::uint64_t;

struct Game {
    std::vector<int> L;
    std::vector<int> R;
};

class GameDB {
public:
    GameDB() {
        add_game({}, {});
    }

    int zero() const { return 0; }

    int canonical(std::vector<int> L, std::vector<int> R) {
        dedup(L);
        dedup(R);

        bool changed = true;
        while (changed) {
            changed = false;

            std::vector<int> keepL;
            keepL.reserve(L.size());
            for (std::size_t i = 0; i < L.size(); ++i) {
                bool dominated = false;
                for (std::size_t j = 0; j < L.size(); ++j) {
                    if (i == j) continue;
                    if (leq(L[i], L[j])) {
                        dominated = true;
                        break;
                    }
                }
                if (!dominated) keepL.push_back(L[i]);
            }
            if (keepL.size() != L.size()) changed = true;
            L.swap(keepL);

            std::vector<int> keepR;
            keepR.reserve(R.size());
            for (std::size_t i = 0; i < R.size(); ++i) {
                bool dominated = false;
                for (std::size_t j = 0; j < R.size(); ++j) {
                    if (i == j) continue;
                    if (leq(R[j], R[i])) {
                        dominated = true;
                        break;
                    }
                }
                if (!dominated) keepR.push_back(R[i]);
            }
            if (keepR.size() != R.size()) changed = true;
            R.swap(keepR);

            int gtemp = add_game(L, R);

            std::vector<int> newL;
            bool revL = false;
            for (int a : L) {
                const auto& AR = games_[a].R;
                bool reversible = false;
                for (int ar : AR) {
                    if (leq(ar, gtemp)) {
                        const auto& ARL = games_[ar].L;
                        newL.insert(newL.end(), ARL.begin(), ARL.end());
                        reversible = true;
                        revL = true;
                        break;
                    }
                }
                if (!reversible) newL.push_back(a);
            }
            if (revL) {
                dedup(newL);
                L.swap(newL);
                changed = true;
            }

            gtemp = add_game(L, R);
            std::vector<int> newR;
            bool revR = false;
            for (int a : R) {
                const auto& AL = games_[a].L;
                bool reversible = false;
                for (int al : AL) {
                    if (leq(gtemp, al)) {
                        const auto& ALR = games_[al].R;
                        newR.insert(newR.end(), ALR.begin(), ALR.end());
                        reversible = true;
                        revR = true;
                        break;
                    }
                }
                if (!reversible) newR.push_back(a);
            }
            if (revR) {
                dedup(newR);
                R.swap(newR);
                changed = true;
            }
        }

        return add_game(L, R);
    }

    int add(int g, int h) {
        if (g > h) std::swap(g, h);
        u64 key = pair_key(g, h);
        auto it = add_cache_.find(key);
        if (it != add_cache_.end()) return it->second;

        std::vector<int> L;
        std::vector<int> R;

        for (int gl : games_[g].L) L.push_back(add(gl, h));
        for (int hl : games_[h].L) L.push_back(add(g, hl));
        for (int gr : games_[g].R) R.push_back(add(gr, h));
        for (int hr : games_[h].R) R.push_back(add(g, hr));

        int res = canonical(std::move(L), std::move(R));
        add_cache_[key] = res;
        return res;
    }

    char cls(int g) {
        bool o = win(g, 0);
        bool e = win(g, 1);
        if (o && e) return 'L';
        if (!o && !e) return 'R';
        if (o && !e) return 'N';
        return 'P';
    }

private:
    std::vector<Game> games_;
    std::unordered_map<std::string, int> intern_;
    std::unordered_map<u64, bool> leq_cache_;
    std::unordered_map<u64, int> add_cache_;
    std::unordered_map<u64, bool> win_cache_;

    static void dedup(std::vector<int>& v) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    static u64 pair_key(int a, int b) {
        return (static_cast<u64>(static_cast<std::uint32_t>(a)) << 32) |
               static_cast<u64>(static_cast<std::uint32_t>(b));
    }

    static std::string encode(const std::vector<int>& L, const std::vector<int>& R) {
        std::string s;
        for (int x : L) {
            s += std::to_string(x);
            s.push_back(',');
        }
        s.push_back('|');
        for (int x : R) {
            s += std::to_string(x);
            s.push_back(',');
        }
        return s;
    }

    int add_game(const std::vector<int>& L, const std::vector<int>& R) {
        std::string key = encode(L, R);
        auto it = intern_.find(key);
        if (it != intern_.end()) return it->second;

        int id = static_cast<int>(games_.size());
        games_.push_back({L, R});
        intern_[std::move(key)] = id;
        return id;
    }

    bool leq(int g, int h) {
        u64 key = pair_key(g, h);
        auto it = leq_cache_.find(key);
        if (it != leq_cache_.end()) return it->second;

        leq_cache_[key] = true;

        for (int gl : games_[g].L) {
            if (leq(h, gl)) {
                leq_cache_[key] = false;
                return false;
            }
        }
        for (int hr : games_[h].R) {
            if (leq(hr, g)) {
                leq_cache_[key] = false;
                return false;
            }
        }
        return true;
    }

    bool win(int g, int turn) {
        u64 key = pair_key(g, turn);
        auto it = win_cache_.find(key);
        if (it != win_cache_.end()) return it->second;

        const auto& moves = (turn == 0) ? games_[g].L : games_[g].R;
        if (moves.empty()) {
            win_cache_[key] = false;
            return false;
        }

        for (int nxt : moves) {
            if (!win(nxt, 1 - turn)) {
                win_cache_[key] = true;
                return true;
            }
        }

        win_cache_[key] = false;
        return false;
    }
};

static u64 solve_C(int N) {
    GameDB db;

    std::vector<int> heap(N + 1, 0);
    heap[0] = db.zero();

    for (int n = 1; n <= N; ++n) {
        if (n & 1) {
            int m = (n - 1) / 2;
            int x = db.add(heap[m], heap[m]);
            heap[n] = db.canonical({x}, {});
        } else {
            int m = (n - 2) / 2;
            int x = db.add(heap[m], heap[m]);
            heap[n] = db.canonical({}, {x});
        }
    }

    std::vector<std::unordered_map<int, u64>> dp(N + 1);
    dp[0][db.zero()] = 1;

    for (int size = 1; size <= N; ++size) {
        int g = heap[size];
        for (int s = size; s <= N; ++s) {
            if (dp[s - size].empty()) continue;
            for (const auto& [gid, cnt] : dp[s - size]) {
                int nxt = db.add(gid, g);
                dp[s][nxt] += cnt;
            }
        }
    }

    u64 ans = 0;
    for (const auto& [gid, cnt] : dp[N]) {
        char c = db.cls(gid);
        if (c == 'P' || c == 'R') ans += cnt;
    }
    return ans;
}

int main() {
    assert(solve_C(5) == 2);
    assert(solve_C(16) == 64);
    std::cout << solve_C(300) << '\n';
    return 0;
}
