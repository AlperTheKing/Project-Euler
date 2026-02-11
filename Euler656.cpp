#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;

constexpr u64 kMod = 1'000'000'000'000'000ULL;

struct SqrtCfTail {
    int D;
    int a0;
    int m;
    int d;
    int a;

    explicit SqrtCfTail(int value) : D(value), a0(0), m(0), d(1), a(0) {
        while ((a0 + 1) * (a0 + 1) <= D) ++a0;
        a = a0;
    }

    int next() {
        m = d * a - m;
        d = (D - m * m) / d;
        a = (a0 + m) / d;
        return a;
    }
};

u64 H_mod(int beta, int g) {
    SqrtCfTail cf(beta);
    std::vector<int> a_terms;
    a_terms.reserve(256);

    auto ensure_a = [&](int idx) {
        while (static_cast<int>(a_terms.size()) < idx) {
            a_terms.push_back(cf.next());
        }
    };

    std::vector<u64> q;
    q.reserve(512);
    q.push_back(0);  // q_{-1}
    q.push_back(1);  // q_0
    int q_max_idx = 0;

    auto ensure_q = [&](int idx) {
        while (q_max_idx < idx) {
            const int next_idx = q_max_idx + 1;
            ensure_a(next_idx);
            const int a_next = a_terms[static_cast<std::size_t>(next_idx - 1)];
            const u64 q_next =
                (static_cast<u64>(a_next) * q[static_cast<std::size_t>(q_max_idx + 1)] +
                 q[static_cast<std::size_t>(q_max_idx)]) %
                kMod;
            q.push_back(q_next);
            ++q_max_idx;
        }
    };

    u64 sum = 0;
    int produced = 0;
    int block = 0;

    while (produced < g) {
        const int i = 2 * block + 1;
        ensure_q(i - 1);
        ensure_a(i);

        const int ai = a_terms[static_cast<std::size_t>(i - 1)];
        const u64 base = q[static_cast<std::size_t>(i - 1)];  // q_{i-2}
        const u64 step = q[static_cast<std::size_t>(i)];      // q_{i-1}

        const int take = std::min(ai, g - produced);
        u64 term = base;
        for (int t = 0; t < take; ++t) {
            term += step;
            if (term >= kMod) term -= kMod;
            sum += term;
            if (sum >= kMod) sum -= kMod;
        }

        produced += take;
        ++block;
    }

    return sum;
}

}  // namespace

int main() {
    assert(H_mod(31, 20) == 150'243'655ULL);
    assert(H_mod(2, 8) == 780ULL);
    assert(H_mod(3, 8) == 14'840ULL);

    u64 total = 0;
    for (int beta = 2; beta <= 1000; ++beta) {
        int r = 0;
        while ((r + 1) * (r + 1) <= beta) ++r;
        if (r * r == beta) continue;
        total += H_mod(beta, 100);
        total %= kMod;
    }

    std::cout << std::setw(15) << std::setfill('0') << total << "\n";
    return 0;
}
