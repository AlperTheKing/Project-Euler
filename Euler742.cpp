#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using namespace std;

using ll = long long;

struct Vector {
    int u;
    int v;
};

bool slope_less(const Vector& lhs, const Vector& rhs) {
    ll left = static_cast<ll>(lhs.v) * rhs.u;
    ll right = static_cast<ll>(rhs.v) * lhs.u;
    if (left != right) return left < right;
    if (lhs.u != rhs.u) return lhs.u < rhs.u;
    return lhs.v < rhs.v;
}

vector<Vector> generate_candidates(int limit) {
    vector<Vector> candidates;
    candidates.reserve(limit * limit);
    for (int u = 1; u <= limit; ++u) {
        for (int v = 1; v <= limit; ++v) {
            if (std::gcd(u, v) == 1) {
                candidates.push_back({u, v});
            }
        }
    }
    sort(candidates.begin(), candidates.end(), slope_less);
    return candidates;
}

ll solve_with_axes(int N, int limit) {
    if (N < 4 || (N - 4) % 4 != 0) return -1;
    int K = (N - 4) / 4;
    if (K == 0) return 1;

    auto candidates = generate_candidates(limit);
    if (static_cast<int>(candidates.size()) < K) return -1;

    int max_w = limit * K;
    const ll kInf = (1LL << 62);

    vector<vector<ll>> dp(K + 1, vector<ll>(max_w + 1, kInf));
    vector<int> max_sum(K + 1, -1);
    dp[0][0] = 1;
    max_sum[0] = 0;

    for (const auto& vec : candidates) {
        int u = vec.u;
        int v = vec.v;

        ll local_cost = 2LL * u * v + 2LL * u + 2LL * v;
        ll b4 = 4LL * v;

        for (int k = K - 1; k >= 0; --k) {
            int limit_sum = max_sum[k];
            if (limit_sum < 0) continue;
            int max_here = min(limit_sum, max_w - u);
            auto& cur = dp[k];
            auto& nxt = dp[k + 1];
            for (int w = 0; w <= max_here; ++w) {
                ll cur_area = cur[w];
                if (cur_area == kInf) continue;
                int new_w = w + u;
                ll new_area = cur_area + local_cost + b4 * w;
                if (new_area < nxt[new_w]) {
                    nxt[new_w] = new_area;
                    if (new_w > max_sum[k + 1]) max_sum[k + 1] = new_w;
                }
            }
        }
    }

    ll best = kInf;
    for (ll area : dp[K]) best = min(best, area);
    return best;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const int limit = 80;

    if (solve_with_axes(4, limit) != 1) {
        cerr << "Validation failed: A(4) != 1\n";
        return 1;
    }
    if (solve_with_axes(8, limit) != 7) {
        cerr << "Validation failed: A(8) != 7\n";
        return 1;
    }
    if (solve_with_axes(40, limit) != 1039) {
        cerr << "Validation failed: A(40) != 1039\n";
        return 1;
    }
    if (solve_with_axes(100, limit) != 17473) {
        cerr << "Validation failed: A(100) != 17473\n";
        return 1;
    }

    cout << solve_with_axes(1000, limit) << "\n";
    return 0;
}
