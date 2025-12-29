#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr int kMod = 1000000007;

struct Dyadic {
    long long num;
    int exp;
    Dyadic(long long n = 0, int e = 0) : num(n), exp(e) { normalize(); }

    void normalize() {
        if (num == 0) {
            exp = 0;
            return;
        }
        while (exp > 0 && (num & 1LL) == 0) {
            num >>= 1;
            --exp;
        }
    }
};

static int cmp_dyadic(const Dyadic& a, const Dyadic& b) {
    if (a.num == b.num && a.exp == b.exp) {
        return 0;
    }
    int e = std::max(a.exp, b.exp);
    __int128 na = static_cast<__int128>(a.num) << (e - a.exp);
    __int128 nb = static_cast<__int128>(b.num) << (e - b.exp);
    if (na < nb) {
        return -1;
    }
    if (na > nb) {
        return 1;
    }
    return 0;
}

static Dyadic add_dyadic(const Dyadic& a, const Dyadic& b) {
    int e = std::max(a.exp, b.exp);
    __int128 na = static_cast<__int128>(a.num) << (e - a.exp);
    __int128 nb = static_cast<__int128>(b.num) << (e - b.exp);
    __int128 s = na + nb;
    Dyadic r(static_cast<long long>(s), e);
    r.normalize();
    return r;
}

static Dyadic neg_dyadic(const Dyadic& a) {
    return Dyadic(-a.num, a.exp);
}

static long long floor_dyadic(const Dyadic& a) {
    if (a.exp == 0) {
        return a.num;
    }
    long long den = 1LL << a.exp;
    if (a.num >= 0) {
        return a.num / den;
    }
    return -(((-a.num) + den - 1) / den);
}

static long long ceil_dyadic(const Dyadic& a) {
    Dyadic na = neg_dyadic(a);
    long long f = floor_dyadic(na);
    return -f;
}

static Dyadic simplest_between(const Dyadic& l, const Dyadic& r) {
    for (int n = 0; n <= 60; ++n) {
        auto floor_scaled = [&](const Dyadic& x) -> long long {
            if (n >= x.exp) {
                __int128 val = static_cast<__int128>(x.num) << (n - x.exp);
                return static_cast<long long>(val);
            }
            int sh = x.exp - n;
            if (x.num >= 0) {
                return static_cast<long long>(static_cast<__int128>(x.num) >> sh);
            }
            __int128 pos = -static_cast<__int128>(x.num);
            __int128 q = (pos + (static_cast<__int128>(1) << sh) - 1) >> sh;
            return static_cast<long long>(-q);
        };
        auto ceil_scaled = [&](const Dyadic& x) -> long long {
            Dyadic nx = neg_dyadic(x);
            return -floor_scaled(nx);
        };

        long long low = floor_scaled(l) + 1;
        long long high = ceil_scaled(r) - 1;
        if (low > high) {
            continue;
        }

        long long t;
        if (low <= 0 && 0 <= high) {
            t = 0;
        } else if (high < 0) {
            t = high;
        } else {
            t = low;
        }
        Dyadic ans(t, n);
        ans.normalize();
        return ans;
    }
    return Dyadic(0, 0);
}

template <typename V>
class FlatHash64 {
public:
    FlatHash64() = default;

    V* find(uint64_t key) {
        if (keys_.empty()) {
            return nullptr;
        }
        std::size_t idx = mix(key) & mask_;
        while (true) {
            uint64_t cur = keys_[idx];
            if (cur == kEmpty) {
                return nullptr;
            }
            if (cur == key) {
                return &values_[idx];
            }
            idx = (idx + 1) & mask_;
        }
    }

    void insert(uint64_t key, V value) {
        if (keys_.empty() || (size_ + 1) * 10 >= keys_.size() * 7) {
            rehash(keys_.empty() ? 1024 : keys_.size() * 2);
        }
        insertNoResize(key, value);
    }

private:
    static constexpr uint64_t kEmpty = std::numeric_limits<uint64_t>::max();

    std::vector<uint64_t> keys_;
    std::vector<V> values_;
    std::size_t mask_ = 0;
    std::size_t size_ = 0;

    static uint64_t mix(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    void insertNoResize(uint64_t key, V value) {
        std::size_t idx = mix(key) & mask_;
        while (true) {
            uint64_t cur = keys_[idx];
            if (cur == kEmpty) {
                keys_[idx] = key;
                values_[idx] = value;
                ++size_;
                return;
            }
            if (cur == key) {
                values_[idx] = value;
                return;
            }
            idx = (idx + 1) & mask_;
        }
    }

    void rehash(std::size_t cap) {
        std::size_t n = 1;
        while (n < cap) {
            n <<= 1;
        }
        std::vector<uint64_t> oldKeys = std::move(keys_);
        std::vector<V> oldVals = std::move(values_);
        keys_.assign(n, kEmpty);
        values_.resize(n);
        mask_ = n - 1;
        size_ = 0;
        for (std::size_t i = 0; i < oldKeys.size(); ++i) {
            if (oldKeys[i] != kEmpty) {
                insertNoResize(oldKeys[i], oldVals[i]);
            }
        }
    }
};

struct GameNode {
    bool isNumber = false;
    bool reduced = false;
    Dyadic val;
    std::vector<int> L;
    std::vector<int> R;
};

struct VecKey {
    std::vector<int> L;
    std::vector<int> R;

    bool operator==(const VecKey& other) const { return L == other.L && R == other.R; }
};

struct VecKeyHash {
    std::size_t operator()(const VecKey& key) const noexcept {
        uint64_t h = 1469598103934665603ULL;
        auto mix = [&](uint64_t x) {
            h ^= x;
            h *= 1099511628211ULL;
        };
        for (int x : key.L) {
            mix(static_cast<uint64_t>(x) + 0x9e3779b97f4a7c15ULL);
        }
        mix(0xA5A5A5A5A5A5A5A5ULL);
        for (int x : key.R) {
            mix(static_cast<uint64_t>(x) + 0xD1B54A32D192ED03ULL);
        }
        return static_cast<std::size_t>(h);
    }
};

struct StopPair {
    Dyadic L;
    Dyadic R;
};

struct NumKey {
    long long num;
    int exp;

    bool operator==(const NumKey& other) const { return num == other.num && exp == other.exp; }
};

struct NumKeyHash {
    std::size_t operator()(const NumKey& key) const noexcept {
        uint64_t h = 1469598103934665603ULL;
        h ^= static_cast<uint64_t>(key.num) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= static_cast<uint64_t>(key.exp) + 0xD1B54A32D192ED03ULL + (h << 6) + (h >> 2);
        return static_cast<std::size_t>(h);
    }
};

class GameEngine {
public:
    GameEngine() { zeroId_ = makeNumber(Dyadic(0, 0)); }

    int zero() const { return zeroId_; }

    int canonical(int id) { return rep(id); }

    int number(const Dyadic& value) { return makeNumber(value); }

    std::size_t nodeCount() const { return nodes_.size(); }

    const GameNode& node(int id) const { return nodes_[id]; }

    int makeBinary(int l, int r) {
        if (l != -1) {
            l = rep(l);
        }
        if (r != -1) {
            r = rep(r);
        }
        uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(l + 1)) << 32) |
                       static_cast<uint32_t>(r + 1);
        if (int* it = binMap_.find(key)) {
            return rep(*it);
        }

        int id = static_cast<int>(nodes_.size());
        nodes_.push_back(GameNode());
        parent_.push_back(id);
        if (l != -1) {
            nodes_[id].L.push_back(l);
        }
        if (r != -1) {
            nodes_[id].R.push_back(r);
        }
        int rid = reduceInPlace(id);
        binMap_.insert(key, rid);
        return rid;
    }

    int add(int a, int b) { return addImpl(a, b, true); }

    int addNoMemo(int a, int b) { return addImpl(a, b, false); }

    StopPair stops(int id) {
        id = rep(id);
        ensureStopSize(id);
        if (stopSeen_[static_cast<std::size_t>(id)]) {
            return stopMemo_[static_cast<std::size_t>(id)];
        }
        StopPair res;
        if (nodes_[id].isNumber) {
            res = StopPair{nodes_[id].val, nodes_[id].val};
        } else {
            bool hasL = false;
            bool hasR = false;
            Dyadic Lval(0, 0);
            Dyadic Rval(0, 0);
            for (int l : nodes_[id].L) {
                StopPair sp = stops(l);
                if (!hasL || cmp_dyadic(sp.R, Lval) > 0) {
                    Lval = sp.R;
                    hasL = true;
                }
            }
            for (int r : nodes_[id].R) {
                StopPair sp = stops(r);
                if (!hasR || cmp_dyadic(sp.L, Rval) < 0) {
                    Rval = sp.L;
                    hasR = true;
                }
            }
            res = StopPair{Lval, Rval};
        }
        stopMemo_[static_cast<std::size_t>(id)] = res;
        stopSeen_[static_cast<std::size_t>(id)] = 1;
        return res;
    }

    bool leftFirstWin(int g) {
        g = rep(g);
        ensureOutcomeSize(g);
        if (leftMemo_[g] != -1) {
            return leftMemo_[g] == 1;
        }

        if (nodes_[g].isNumber) {
            leftMemo_[g] = (nodes_[g].val.num > 0) ? 1 : 0;
            rightMemo_[g] = (nodes_[g].val.num < 0) ? 1 : 0;
            return leftMemo_[g] == 1;
        }

        for (int l : nodes_[g].L) {
            if (!rightFirstWin(l)) {
                leftMemo_[g] = 1;
                return true;
            }
        }
        leftMemo_[g] = 0;
        return false;
    }

    bool rightFirstWin(int g) {
        g = rep(g);
        ensureOutcomeSize(g);
        if (rightMemo_[g] != -1) {
            return rightMemo_[g] == 1;
        }

        if (nodes_[g].isNumber) {
            leftMemo_[g] = (nodes_[g].val.num > 0) ? 1 : 0;
            rightMemo_[g] = (nodes_[g].val.num < 0) ? 1 : 0;
            return rightMemo_[g] == 1;
        }

        for (int r : nodes_[g].R) {
            if (!leftFirstWin(r)) {
                rightMemo_[g] = 1;
                return true;
            }
        }
        rightMemo_[g] = 0;
        return false;
    }

private:
    std::vector<GameNode> nodes_;
    std::vector<int> parent_;

    std::unordered_map<NumKey, int, NumKeyHash> numMap_;
    FlatHash64<int> binMap_;
    std::unordered_map<VecKey, int, VecKeyHash> canonMap_;
    FlatHash64<int> addMap_;
    FlatHash64<uint8_t> leqMemo_;
    std::vector<int8_t> leftMemo_;
    std::vector<int8_t> rightMemo_;
    std::vector<StopPair> stopMemo_;
    std::vector<char> stopSeen_;

    int zeroId_ = -1;

    int rep(int x) {
        while (parent_[x] != x) {
            parent_[x] = parent_[parent_[x]];
            x = parent_[x];
        }
        return x;
    }

    void aliasTo(int from, int to) {
        from = rep(from);
        to = rep(to);
        if (from != to) {
            parent_[from] = to;
        }
    }

    void ensureStopSize(int id) {
        if (stopMemo_.size() <= static_cast<std::size_t>(id)) {
            stopMemo_.resize(static_cast<std::size_t>(id) + 1);
            stopSeen_.resize(static_cast<std::size_t>(id) + 1, 0);
        }
    }

    int makeNumber(const Dyadic& in) {
        Dyadic d = in;
        d.normalize();
        NumKey key{d.num, d.exp};
        auto it = numMap_.find(key);
        if (it != numMap_.end()) {
            return it->second;
        }

        int id = static_cast<int>(nodes_.size());
        nodes_.push_back(GameNode());
        parent_.push_back(id);
        nodes_[id].isNumber = true;
        nodes_[id].reduced = true;
        nodes_[id].val = d;
        numMap_[key] = id;

        if (d.num == 0) {
            return id;
        }
        if (d.exp == 0) {
            if (d.num > 0) {
                nodes_[id].L.push_back(makeNumber(Dyadic(d.num - 1, 0)));
            } else {
                nodes_[id].R.push_back(makeNumber(Dyadic(d.num + 1, 0)));
            }
        } else {
            nodes_[id].L.push_back(makeNumber(Dyadic(d.num - 1, d.exp)));
            nodes_[id].R.push_back(makeNumber(Dyadic(d.num + 1, d.exp)));
        }
        return id;
    }

    void ensureOutcomeSize(int g) {
        if (static_cast<int>(leftMemo_.size()) <= g) {
            leftMemo_.resize(g + 1, -1);
            rightMemo_.resize(g + 1, -1);
        }
    }

    bool leq(int a, int b) {
        a = rep(a);
        b = rep(b);
        if (a == b) {
            return true;
        }
        if (nodes_[a].isNumber && nodes_[b].isNumber) {
            return cmp_dyadic(nodes_[a].val, nodes_[b].val) <= 0;
        }

        uint64_t key = (static_cast<uint64_t>(a) << 32) | static_cast<uint32_t>(b);
        if (uint8_t* it = leqMemo_.find(key)) {
            return *it == 1;
        }

        for (int al : nodes_[a].L) {
            if (leq(b, al)) {
                leqMemo_.insert(key, 2);
                return false;
            }
        }
        for (int br : nodes_[b].R) {
            if (leq(br, a)) {
                leqMemo_.insert(key, 2);
                return false;
            }
        }
        leqMemo_.insert(key, 1);
        return true;
    }

    void removeDominatedLeft(std::vector<int>& L) {
        std::vector<char> kill(L.size(), 0);
        for (std::size_t i = 0; i < L.size(); ++i) {
            if (kill[i]) {
                continue;
            }
            for (std::size_t j = 0; j < L.size(); ++j) {
                if (i == j || kill[j]) {
                    continue;
                }
                if (leq(L[i], L[j]) && !leq(L[j], L[i])) {
                    kill[i] = 1;
                    break;
                }
            }
        }
        std::vector<int> out;
        out.reserve(L.size());
        for (std::size_t i = 0; i < L.size(); ++i) {
            if (!kill[i]) {
                out.push_back(L[i]);
            }
        }
        L.swap(out);
    }

    void removeDominatedRight(std::vector<int>& R) {
        std::vector<char> kill(R.size(), 0);
        for (std::size_t i = 0; i < R.size(); ++i) {
            if (kill[i]) {
                continue;
            }
            for (std::size_t j = 0; j < R.size(); ++j) {
                if (i == j || kill[j]) {
                    continue;
                }
                if (leq(R[j], R[i]) && !leq(R[i], R[j])) {
                    kill[i] = 1;
                    break;
                }
            }
        }
        std::vector<int> out;
        out.reserve(R.size());
        for (std::size_t i = 0; i < R.size(); ++i) {
            if (!kill[i]) {
                out.push_back(R[i]);
            }
        }
        R.swap(out);
    }

    int reduceInPlace(int id) {
        id = rep(id);
        if (nodes_[id].isNumber || nodes_[id].reduced) {
            return id;
        }

        auto& L = nodes_[id].L;
        auto& R = nodes_[id].R;

        for (int& x : L) {
            x = rep(x);
        }
        for (int& x : R) {
            x = rep(x);
        }
        std::sort(L.begin(), L.end());
        L.erase(std::unique(L.begin(), L.end()), L.end());
        std::sort(R.begin(), R.end());
        R.erase(std::unique(R.begin(), R.end()), R.end());

        removeDominatedLeft(L);
        removeDominatedRight(R);

        bool changed = true;
        while (changed) {
            changed = false;

            for (std::size_t i = 0; i < L.size() && !changed; ++i) {
                int l = rep(L[i]);
                for (int lr : nodes_[l].R) {
                    if (leq(lr, id)) {
                        L.erase(L.begin() + static_cast<long long>(i));
                        int lr_rep = rep(lr);
                        for (int x : nodes_[lr_rep].L) {
                            L.push_back(rep(x));
                        }
                        changed = true;
                        break;
                    }
                }
            }
            if (changed) {
                std::sort(L.begin(), L.end());
                L.erase(std::unique(L.begin(), L.end()), L.end());
                removeDominatedLeft(L);
                continue;
            }

            for (std::size_t i = 0; i < R.size() && !changed; ++i) {
                int r = rep(R[i]);
                for (int rl : nodes_[r].L) {
                    if (leq(id, rl)) {
                        R.erase(R.begin() + static_cast<long long>(i));
                        int rl_rep = rep(rl);
                        for (int x : nodes_[rl_rep].R) {
                            R.push_back(rep(x));
                        }
                        changed = true;
                        break;
                    }
                }
            }
            if (changed) {
                std::sort(R.begin(), R.end());
                R.erase(std::unique(R.begin(), R.end()), R.end());
                removeDominatedRight(R);
            }
        }

        if (L.empty() && R.empty()) {
            aliasTo(id, zeroId_);
            return zeroId_;
        }

        auto allNumbers = [&](const std::vector<int>& v) -> bool {
            for (int x : v) {
                if (!nodes_[rep(x)].isNumber) {
                    return false;
                }
            }
            return true;
        };

        bool lNum = allNumbers(L);
        bool rNum = allNumbers(R);

        if (!L.empty() && !R.empty() && lNum && rNum) {
            Dyadic maxL = nodes_[rep(L[0])].val;
            for (int x : L) {
                Dyadic v = nodes_[rep(x)].val;
                if (cmp_dyadic(v, maxL) > 0) {
                    maxL = v;
                }
            }
            Dyadic minR = nodes_[rep(R[0])].val;
            for (int x : R) {
                Dyadic v = nodes_[rep(x)].val;
                if (cmp_dyadic(v, minR) < 0) {
                    minR = v;
                }
            }
            if (cmp_dyadic(maxL, minR) < 0) {
                Dyadic mid = simplest_between(maxL, minR);
                int nid = makeNumber(mid);
                aliasTo(id, nid);
                return nid;
            }
        } else if (lNum && R.empty()) {
            Dyadic maxL = nodes_[rep(L[0])].val;
            for (int x : L) {
                Dyadic v = nodes_[rep(x)].val;
                if (cmp_dyadic(v, maxL) > 0) {
                    maxL = v;
                }
            }
            long long fl = floor_dyadic(maxL);
            int nid = makeNumber(Dyadic(fl + 1, 0));
            aliasTo(id, nid);
            return nid;
        } else if (rNum && L.empty()) {
            Dyadic minR = nodes_[rep(R[0])].val;
            for (int x : R) {
                Dyadic v = nodes_[rep(x)].val;
                if (cmp_dyadic(v, minR) < 0) {
                    minR = v;
                }
            }
            long long ce = ceil_dyadic(minR);
            int nid = makeNumber(Dyadic(ce - 1, 0));
            aliasTo(id, nid);
            return nid;
        }

        VecKey key{L, R};
        auto it = canonMap_.find(key);
        if (it != canonMap_.end()) {
            aliasTo(id, it->second);
            return rep(it->second);
        }
        canonMap_.emplace(std::move(key), id);
        nodes_[id].reduced = true;
        return id;
    }

    int makeGame(std::vector<int> L, std::vector<int> R) {
        for (int& x : L) {
            x = rep(x);
        }
        for (int& x : R) {
            x = rep(x);
        }
        std::sort(L.begin(), L.end());
        L.erase(std::unique(L.begin(), L.end()), L.end());
        std::sort(R.begin(), R.end());
        R.erase(std::unique(R.begin(), R.end()), R.end());

        if (L.size() <= 1 && R.size() <= 1) {
            int l = L.empty() ? -1 : L[0];
            int r = R.empty() ? -1 : R[0];
            return makeBinary(l, r);
        }

        int id = static_cast<int>(nodes_.size());
        nodes_.push_back(GameNode());
        parent_.push_back(id);
        nodes_[id].L = std::move(L);
        nodes_[id].R = std::move(R);
        return reduceInPlace(id);
    }

    int addImpl(int a, int b, bool cacheThis) {
        a = rep(a);
        b = rep(b);
        if (a == zeroId_) {
            return b;
        }
        if (b == zeroId_) {
            return a;
        }
        if (nodes_[a].isNumber && nodes_[b].isNumber) {
            return makeNumber(add_dyadic(nodes_[a].val, nodes_[b].val));
        }

        int x = a;
        int y = b;
        if (x > y) {
            std::swap(x, y);
        }
        uint64_t key = (static_cast<uint64_t>(x) << 32) | static_cast<uint32_t>(y);
        if (int* it = addMap_.find(key)) {
            return rep(*it);
        }

        const std::vector<int> aL = nodes_[a].L;
        const std::vector<int> aR = nodes_[a].R;
        const std::vector<int> bL = nodes_[b].L;
        const std::vector<int> bR = nodes_[b].R;
        std::vector<int> L;
        std::vector<int> R;
        L.reserve(aL.size() + bL.size());
        R.reserve(aR.size() + bR.size());

        if (!nodes_[a].isNumber) {
            for (int al : aL) {
                L.push_back(addImpl(al, b, true));
            }
            for (int ar : aR) {
                R.push_back(addImpl(ar, b, true));
            }
        }
        if (!nodes_[b].isNumber) {
            for (int bl : bL) {
                L.push_back(addImpl(a, bl, true));
            }
            for (int br : bR) {
                R.push_back(addImpl(a, br, true));
            }
        }

        int res = makeGame(std::move(L), std::move(R));
        if (cacheThis) {
            addMap_.insert(key, res);
        }
        return res;
    }
};

struct StaircaseFreq {
    std::unordered_map<int, int> numFreq;  // value -> count
    std::vector<std::unordered_map<int, int>> switchFreq;  // weight -> (R -> count)
    int maxWeight = 0;
};

static StaircaseFreq collect_staircase_freq(int max_w, GameEngine& eng) {
    StaircaseFreq freq;
    freq.switchFreq.resize(static_cast<std::size_t>(max_w + 1));

    for (int a = 1; a <= max_w - 2; ++a) {
        for (int b = 1; b <= max_w - a - 1; ++b) {
            int K = max_w - a - b;
            if (K < 1) {
                continue;
            }

            int A = a;
            int B = b;
            std::vector<int> dp(static_cast<std::size_t>((K + 1) * A * B), eng.zero());
            auto idx = [&](int k, int i, int j) { return (k * A + i) * B + j; };

            for (int k = 1; k <= K; ++k) {
                for (int i = A - 1; i >= 0; --i) {
                    for (int j = B - 1; j >= 0; --j) {
                        int left = -1;
                        int right = -1;

                        if (k == 1 && j == B - 1) {
                            left = -1;
                        } else if (j < B - 1) {
                            left = dp[static_cast<std::size_t>(idx(k, i, j + 1))];
                        } else {
                            left = dp[static_cast<std::size_t>(idx(k - 1, i, 0))];
                        }

                        if (k == 1 && i == A - 1) {
                            right = -1;
                        } else if (i < A - 1) {
                            right = dp[static_cast<std::size_t>(idx(k, i + 1, j))];
                        } else {
                            right = dp[static_cast<std::size_t>(idx(k - 1, 0, j))];
                        }

                        dp[static_cast<std::size_t>(idx(k, i, j))] = eng.makeBinary(left, right);
                    }
                }
            }

            for (int k = 1; k <= K; ++k) {
                int id = eng.canonical(dp[static_cast<std::size_t>(idx(k, 0, 0))]);
                StopPair sp = eng.stops(id);
                if (sp.L.exp != 0 || sp.R.exp != 0) {
                    std::cerr << "Non-integer stops encountered.\n";
                    std::exit(1);
                }
                int L = static_cast<int>(sp.L.num);
                int R = static_cast<int>(sp.R.num);
                const auto& node = eng.node(id);
                if (node.isNumber) {
                    int& slot = freq.numFreq[L];
                    ++slot;
                    if (slot >= kMod) {
                        slot -= kMod;
                    }
                } else {
                    if (node.L.size() != 1 || node.R.size() != 1) {
                        std::cerr << "Non-switch game encountered.\n";
                        std::exit(1);
                    }
                    int l = eng.canonical(node.L[0]);
                    int r = eng.canonical(node.R[0]);
                    if (!eng.node(l).isNumber || !eng.node(r).isNumber) {
                        std::cerr << "Switch options are non-numbers.\n";
                        std::exit(1);
                    }
                    int weight = L - R;
                    if (weight < 0) {
                        std::cerr << "Negative switch weight encountered.\n";
                        std::exit(1);
                    }
                    if (weight >= static_cast<int>(freq.switchFreq.size())) {
                        freq.switchFreq.resize(static_cast<std::size_t>(weight + 1));
                    }
                    int& slot = freq.switchFreq[static_cast<std::size_t>(weight)][R];
                    ++slot;
                    if (slot >= kMod) {
                        slot -= kMod;
                    }
                    freq.maxWeight = std::max(freq.maxWeight, weight);
                }
            }
        }
    }

    return freq;
}

static std::vector<std::unordered_map<int, int>> compute_dist(
    const std::unordered_map<int, int>& base_freq,
    int max_len) {
    std::vector<std::unordered_map<int, int>> dist(static_cast<std::size_t>(max_len + 1));
    dist[0][0] = 1;
    if (base_freq.empty()) {
        return dist;
    }

    for (int len = 1; len <= max_len; ++len) {
        auto& cur = dist[len];
        auto& prev = dist[len - 1];
        for (const auto& kv : prev) {
            int sum = kv.first;
            int cnt = kv.second;
            for (const auto& bv : base_freq) {
                int val = bv.first;
                int freq = bv.second;
                int& slot = cur[sum + val];
                long long add = static_cast<long long>(cnt) * freq % kMod;
                int nv = slot + static_cast<int>(add);
                if (nv >= kMod) {
                    nv -= kMod;
                }
                slot = nv;
            }
        }
    }
    return dist;
}

static uint64_t pack_key(int sumR, int sumOdd) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(sumR)) << 32) |
           static_cast<uint32_t>(sumOdd);
}

static int unpack_sumR(uint64_t key) { return static_cast<int32_t>(key >> 32); }
static int unpack_sumOdd(uint64_t key) { return static_cast<int>(key & 0xFFFFFFFFu); }

static int compute_S(int m, int max_w, GameEngine& eng) {
    StaircaseFreq freq = collect_staircase_freq(max_w, eng);

    std::vector<std::vector<int>> binom(m + 1, std::vector<int>(m + 1, 0));
    for (int n = 0; n <= m; ++n) {
        binom[n][0] = binom[n][n] = 1;
        for (int k = 1; k < n; ++k) {
            int nv = binom[n - 1][k - 1] + binom[n - 1][k];
            if (nv >= kMod) {
                nv -= kMod;
            }
            binom[n][k] = nv;
        }
    }

    auto num_dist = compute_dist(freq.numFreq, m);

    std::vector<std::vector<std::unordered_map<int, int>>> dist_by_weight;
    dist_by_weight.resize(static_cast<std::size_t>(freq.maxWeight + 1));
    for (int w = 0; w <= freq.maxWeight; ++w) {
        if (freq.switchFreq[static_cast<std::size_t>(w)].empty()) {
            continue;
        }
        dist_by_weight[static_cast<std::size_t>(w)] =
            compute_dist(freq.switchFreq[static_cast<std::size_t>(w)], m);
    }

    using StateMap = std::unordered_map<uint64_t, int>;
    std::vector<std::array<StateMap, 2>> dp(static_cast<std::size_t>(m + 1));
    dp[0][0][pack_key(0, 0)] = 1;

    for (int w = freq.maxWeight; w >= 0; --w) {
        const auto& dist_w = dist_by_weight[static_cast<std::size_t>(w)];
        if (dist_w.empty()) {
            continue;
        }
        std::vector<std::array<StateMap, 2>> next(static_cast<std::size_t>(m + 1));

        for (int s = 0; s <= m; ++s) {
            for (int parity = 0; parity <= 1; ++parity) {
                const auto& cur = dp[s][parity];
                if (cur.empty()) {
                    continue;
                }
                for (const auto& kv : cur) {
                    int sumR = unpack_sumR(kv.first);
                    int sumOdd = unpack_sumOdd(kv.first);
                    int cnt = kv.second;
                    for (int c = 0; c <= m - s; ++c) {
                        const auto& dist_c = dist_w[static_cast<std::size_t>(c)];
                        if (dist_c.empty()) {
                            continue;
                        }
                        int leftCount = (c + (parity == 0 ? 1 : 0)) / 2;
                        int newParity = parity ^ (c & 1);
                        int comb = binom[s + c][c];
                        for (const auto& rv : dist_c) {
                            int newS = s + c;
                            int newSumR = sumR + rv.first;
                            int newSumOdd = sumOdd + leftCount * w;
                            long long ways = static_cast<long long>(cnt) * rv.second % kMod;
                            ways = ways * comb % kMod;
                            auto& slot = next[newS][newParity][pack_key(newSumR, newSumOdd)];
                            int nv = slot + static_cast<int>(ways);
                            if (nv >= kMod) {
                                nv -= kMod;
                            }
                            slot = nv;
                        }
                    }
                }
            }
        }
        dp.swap(next);
    }

    long long ans = 0;
    for (int s = 0; s <= m; ++s) {
        int n = m - s;
        int interleave = binom[m][s];
        for (int parity = 0; parity <= 1; ++parity) {
            const auto& cur = dp[s][parity];
            if (cur.empty()) {
                continue;
            }
            for (const auto& kv : cur) {
                int sumR = unpack_sumR(kv.first);
                int sumOdd = unpack_sumOdd(kv.first);
                int cnt = kv.second;
                const auto& num_map = num_dist[static_cast<std::size_t>(n)];
                for (const auto& nv : num_map) {
                    long long total = static_cast<long long>(sumR) + sumOdd + nv.first;
                    bool win;
                    if (total > 0) {
                        win = true;
                    } else if (total < 0) {
                        win = false;
                    } else {
                        win = (s % 2 == 1);
                    }
                    if (!win) {
                        continue;
                    }
                    long long ways = static_cast<long long>(cnt) * nv.second % kMod;
                    ways = ways * interleave % kMod;
                    ans += ways;
                    if (ans >= kMod) {
                        ans %= kMod;
                    }
                }
            }
        }
    }

    return static_cast<int>(ans % kMod);
}

}  // namespace

int main() {
    GameEngine eng;

    int s24 = compute_S(2, 4, eng);
    if (s24 != 7) {
        std::cerr << "Validation failed: S(2,4) != 7, got " << s24 << "\n";
        return 1;
    }

    int s39 = compute_S(3, 9, eng);
    if (s39 != 315319) {
        std::cerr << "Validation failed: S(3,9) != 315319, got " << s39 << "\n";
        return 1;
    }

    std::cerr << "Validation checkpoints passed.\n";

    int s864 = compute_S(8, 64, eng);
    std::cout << s864 << "\n";
    return 0;
}