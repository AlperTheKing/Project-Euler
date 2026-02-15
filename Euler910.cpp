#include <cstdint>
#include <memory>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

static constexpr uint32_t MOD = 1000000000u;

static uint32_t mul_mod(uint32_t a, uint32_t b) {
    return static_cast<uint32_t>((static_cast<uint64_t>(a) * b) % MOD);
}

static uint32_t pow_mod(uint32_t base, uint32_t exp) {
    uint64_t result = 1;
    uint64_t cur = base % MOD;
    while (exp > 0) {
        if (exp & 1u) {
            result = (result * cur) % MOD;
        }
        cur = (cur * cur) % MOD;
        exp >>= 1u;
    }
    return static_cast<uint32_t>(result);
}

struct Func {
    virtual uint32_t eval(uint32_t x) = 0;
    virtual ~Func() = default;
};

struct FlatMap {
    vector<uint32_t> keys;
    vector<uint32_t> values;
    size_t used = 0;

    static constexpr uint32_t kEmpty = 0xFFFFFFFFu;

    static uint32_t hash(uint32_t x) {
        x ^= x >> 16;
        x *= 0x7feb352d;
        x ^= x >> 15;
        x *= 0x846ca68b;
        x ^= x >> 16;
        return x;
    }

    void rehash(size_t cap) {
        vector<uint32_t> old_keys = std::move(keys);
        vector<uint32_t> old_values = std::move(values);
        keys.assign(cap, kEmpty);
        values.assign(cap, 0);
        used = 0;
        for (size_t i = 0; i < old_keys.size(); ++i) {
            if (old_keys[i] != kEmpty) {
                insert(old_keys[i], old_values[i]);
            }
        }
    }

    void reserve(size_t cap) {
        size_t pow2 = 1;
        while (pow2 < cap) {
            pow2 <<= 1;
        }
        if (pow2 > keys.size()) {
            rehash(pow2);
        }
    }

    bool get(uint32_t key, uint32_t& out) const {
        if (keys.empty()) {
            return false;
        }
        const size_t mask = keys.size() - 1;
        size_t idx = hash(key) & mask;
        while (true) {
            const uint32_t cur = keys[idx];
            if (cur == kEmpty) {
                return false;
            }
            if (cur == key) {
                out = values[idx];
                return true;
            }
            idx = (idx + 1) & mask;
        }
    }

    void insert(uint32_t key, uint32_t value) {
        if (keys.empty() || (used + 1) * 10 >= keys.size() * 7) {
            const size_t next = keys.empty() ? 64 : keys.size() * 2;
            rehash(next);
        }
        const size_t mask = keys.size() - 1;
        size_t idx = hash(key) & mask;
        while (true) {
            const uint32_t cur = keys[idx];
            if (cur == kEmpty) {
                keys[idx] = key;
                values[idx] = value;
                ++used;
                return;
            }
            if (cur == key) {
                values[idx] = value;
                return;
            }
            idx = (idx + 1) & mask;
        }
    }
};

struct IterFunc {
    Func* func;
    uint32_t steps;
    int maxbit;
    vector<FlatMap> jump;

    IterFunc(Func* f, uint32_t steps) : func(f), steps(steps) {
        maxbit = 0;
        uint32_t tmp = steps;
        while (tmp > 0) {
            ++maxbit;
            tmp >>= 1u;
        }
        jump.resize(maxbit);
    }

    uint32_t jump_k(int k, uint32_t x) {
        auto& mp = jump[k];
        uint32_t cached;
        if (mp.get(x, cached)) {
            return cached;
        }
        uint32_t y;
        if (k == 0) {
            y = func->eval(x);
        } else {
            uint32_t mid = jump_k(k - 1, x);
            y = jump_k(k - 1, mid);
        }
        mp.insert(x, y);
        return y;
    }

    uint32_t apply(uint32_t x) {
        if (steps == 0) {
            return x;
        }
        uint32_t mask = steps;
        int bit = 0;
        while (mask > 0) {
            if (mask & 1u) {
                x = jump_k(bit, x);
            }
            mask >>= 1u;
            ++bit;
        }
        return x;
    }
};

struct FuncC : Func {
    uint32_t c;
    FlatMap cache;
    explicit FuncC(uint32_t c) : c(c) {}

    uint32_t eval(uint32_t x) override {
        uint32_t cached;
        if (cache.get(x, cached)) {
            return cached;
        }
        uint32_t res = mul_mod(x + 1u, pow_mod(x, c));
        cache.insert(x, res);
        return res;
    }
};

struct FuncLevel0 : Func {
    FuncC* base;
    IterFunc* iter_base;
    FlatMap cache;

    FuncLevel0(FuncC* base, IterFunc* iter_base) : base(base), iter_base(iter_base) {}

    uint32_t eval(uint32_t x) override {
        uint32_t cached;
        if (cache.get(x, cached)) {
            return cached;
        }
        uint32_t fx = base->eval(x);
        uint32_t start = mul_mod(x, fx);
        uint32_t w = iter_base->apply(start);
        uint32_t res = base->eval(w);
        cache.insert(x, res);
        return res;
    }
};

struct FuncLevelN : Func {
    Func* prev;
    IterFunc* iter_prev;
    FlatMap cache;

    FuncLevelN(Func* prev, IterFunc* iter_prev) : prev(prev), iter_prev(iter_prev) {}

    uint32_t eval(uint32_t x) override {
        uint32_t cached;
        if (cache.get(x, cached)) {
            return cached;
        }
        uint32_t fx = prev->eval(x);
        uint32_t start = mul_mod(x, fx);
        uint32_t res = iter_prev->apply(start);
        cache.insert(x, res);
        return res;
    }
};

class Evaluator {
public:
    Evaluator(uint32_t a, uint32_t b, uint32_t c) : a_(a), b_(b), c_(c), base_(c_), iter_base_(&base_, b_) {
        funcs_.reserve(a_ + 1);
        iters_.reserve(a_ + 1);
        base_.cache.reserve(1 << 12);

        funcs_.push_back(make_unique<FuncLevel0>(&base_, &iter_base_));
        iters_.push_back(make_unique<IterFunc>(funcs_[0].get(), b_));

        for (uint32_t i = 1; i <= a_; ++i) {
            funcs_.push_back(make_unique<FuncLevelN>(funcs_[i - 1].get(), iters_[i - 1].get()));
            iters_.push_back(make_unique<IterFunc>(funcs_[i].get(), b_));
        }
    }

    uint32_t compute(uint32_t d) {
        return funcs_[a_]->eval(d);
    }

private:
    uint32_t a_;
    uint32_t b_;
    uint32_t c_;
    FuncC base_;
    IterFunc iter_base_;
    vector<unique_ptr<Func>> funcs_;
    vector<unique_ptr<IterFunc>> iters_;
};

static bool run_validation() {
    struct Case {
        uint32_t a, b, c, d;
        uint32_t expected;
    };

    const Case cases[] = {
        {0, 0, 0, 2, 7},
        {0, 1, 0, 2, 8},
        {0, 0, 1, 1, 6},
        {1, 0, 0, 2, 14},
        {1, 1, 0, 2, 274},
    };

    for (const auto& tc : cases) {
        Evaluator eval(tc.a, tc.b, tc.c);
        uint32_t got = eval.compute(tc.d);
        if (got != tc.expected) {
            cerr << "Validation failed for a=" << tc.a << " b=" << tc.b << " c=" << tc.c
                 << " d=" << tc.d << " expected=" << tc.expected << " got=" << got << "\n";
            return false;
        }
    }

    return true;
}

int main() {
    const uint32_t a = 12;
    const uint32_t b = 345678;
    const uint32_t c = 9012345;
    const uint32_t d = 678;
    const uint32_t e = 90;

    if (!run_validation()) {
        return 1;
    }

    Evaluator eval(a, b, c);
    uint32_t k = eval.compute(d);
    uint32_t result = (k + e) % MOD;

    cout << setw(9) << setfill('0') << result << "\n";
    return 0;
}
