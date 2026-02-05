#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_set>
#include <vector>

using namespace std;

namespace {

struct Entry {
    long double value;
    uint16_t a;
    uint16_t c;
};

struct ArrayHash {
    size_t operator()(const array<int, 8>& v) const noexcept {
        size_t h = 0;
        for (int x : v) {
            h = h * 181u + static_cast<unsigned>(x);
        }
        return h;
    }
};

array<int, 8> canonical(const array<int, 8>& seq) {
    array<pair<int, int>, 4> pairs = {
        make_pair(seq[0], seq[1]),
        make_pair(seq[2], seq[3]),
        make_pair(seq[4], seq[5]),
        make_pair(seq[6], seq[7]),
    };

    array<int, 8> best{};
    bool first = true;

    auto consider = [&](const array<pair<int, int>, 4>& ps, int rot) {
        array<int, 8> s;
        for (int i = 0; i < 4; ++i) {
            const auto& p = ps[(rot + i) & 3];
            s[2 * i] = p.first;
            s[2 * i + 1] = p.second;
        }
        if (first || s < best) {
            best = s;
            first = false;
        }
    };

    for (int r = 0; r < 4; ++r) consider(pairs, r);

    array<pair<int, int>, 4> ref = {
        make_pair(pairs[0].second, pairs[0].first),
        make_pair(pairs[3].second, pairs[3].first),
        make_pair(pairs[2].second, pairs[2].first),
        make_pair(pairs[1].second, pairs[1].first),
    };
    for (int r = 0; r < 4; ++r) consider(ref, r);

    return best;
}

uint64_t solve() {
    const long double pi = acosl(-1.0L);
    array<long double, 181> sinv{};
    for (int d = 0; d <= 180; ++d) {
        sinv[d] = sinl(pi * static_cast<long double>(d) / 180.0L);
    }

    vector<vector<long double>> g(181);
    for (int n = 2; n <= 178; ++n) {
        g[n].resize(n);
        for (int k = 1; k <= n - 1; ++k) {
            g[n][k] = sinv[n - k] / sinv[k];
        }
    }

    vector<vector<Entry>> ratio_list(181);
    for (int U = 2; U <= 178; ++U) {
        int V = 180 - U;
        auto& vec = ratio_list[U];
        vec.reserve(static_cast<size_t>(U - 1) * static_cast<size_t>(V - 1));
        for (int c = 1; c <= U - 1; ++c) {
            long double gU = g[U][c];
            for (int a = 1; a <= V - 1; ++a) {
                long double ratio = gU / g[V][a];
                vec.push_back({ratio, static_cast<uint16_t>(a), static_cast<uint16_t>(c)});
            }
        }
        sort(vec.begin(), vec.end(), [](const Entry& lhs, const Entry& rhs) {
            return lhs.value < rhs.value;
        });
    }

    unordered_set<array<int, 8>, ArrayHash> unique;
    unique.reserve(200000);

    const long double search_eps = 1e-10L;
    const long double verify_eps = 1e-12L;

    for (int p = 1; p <= 178; ++p) {
        for (int q = 1; q <= 178 - p; ++q) {
            int U = p + q;
            int V = 180 - U;
            auto& vec = ratio_list[U];
            long double sinp = sinv[p];
            long double sinq = sinv[q];
            for (int r = 1; r <= V - 1; ++r) {
                int s = V - r;
                long double K = (sinp * sinv[r]) / (sinq * sinv[s]);
                long double low = K - search_eps;
                long double high = K + search_eps;

                auto it1 = lower_bound(vec.begin(), vec.end(), low,
                                       [](const Entry& e, long double v) {
                                           return e.value < v;
                                       });
                auto it2 = upper_bound(vec.begin(), vec.end(), high,
                                       [](long double v, const Entry& e) {
                                           return v < e.value;
                                       });
                for (auto it = it1; it != it2; ++it) {
                    int a = it->a;
                    int c = it->c;
                    long double left = sinq * sinv[s] * sinv[U - c] * sinv[a];
                    long double right = sinp * sinv[r] * sinv[c] * sinv[V - a];
                    long double diff = fabsl(left - right);
                    long double scale = max(1.0L, fabsl(left));
                    if (diff > verify_eps * scale) continue;

                    array<int, 8> seq = {a, p, q, r, s, c, U - c, V - a};
                    unique.insert(canonical(seq));
                }
            }
        }
    }

    return unique.size();
}

}  // namespace

int main() {
    cout << solve() << "\n";
    return 0;
}
