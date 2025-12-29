#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstdint>
#include <deque>
#include <iostream>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

namespace {

using u128 = unsigned __int128;

u128 euclid_sum(u128 a, u128 b) {
    if (a < b) swap(a, b);
    u128 sum = 0;
    while (b != 0) {
        sum += a / b;
        u128 r = a % b;
        a = b;
        b = r;
    }
    return sum;
}

string to_string_u128(u128 v) {
    if (v == 0) return "0";
    string s;
    while (v > 0) {
        int digit = static_cast<int>(v % 10);
        s.push_back(static_cast<char>('0' + digit));
        v /= 10;
    }
    reverse(s.begin(), s.end());
    return s;
}

struct Triple {
    int a;
    int b;
    int c;
};

struct TripleHash {
    size_t operator()(const Triple& t) const noexcept {
        size_t h1 = std::hash<int>{}(t.a);
        size_t h2 = std::hash<int>{}(t.b);
        size_t h3 = std::hash<int>{}(t.c);
        size_t h = h1;
        h ^= h2 + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= h3 + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

bool operator==(const Triple& x, const Triple& y) {
    return x.a == y.a && x.b == y.b && x.c == y.c;
}

int min_steps_bfs(int a, int b, int c, int bound, int depth_limit) {
    if (a == 0 || b == 0 || c == 0) return 0;
    deque<pair<Triple, int>> q;
    unordered_set<Triple, TripleHash> seen;

    q.push_back({{a, b, c}, 0});
    seen.insert({a, b, c});

    while (!q.empty()) {
        auto [cur, dist] = q.front();
        q.pop_front();
        if (dist >= depth_limit) continue;

        const int aa = cur.a;
        const int bb = cur.b;
        const int cc = cur.c;

        const Triple nexts[3] = {
            {2 * (bb + cc) - aa, bb, cc},
            {aa, 2 * (cc + aa) - bb, cc},
            {aa, bb, 2 * (aa + bb) - cc}
        };

        for (const auto& nxt : nexts) {
            if (nxt.a == 0 || nxt.b == 0 || nxt.c == 0) return dist + 1;
            if (max({abs(nxt.a), abs(nxt.b), abs(nxt.c)}) > bound) continue;
            if (seen.insert(nxt).second) {
                q.push_back({nxt, dist + 1});
            }
        }
    }
    return 0;
}

u128 compute_F_bruteforce(int a, int b) {
    const int N = 4 * a * b;
    unordered_set<int> candidates;
    for (int d = 1; d * d <= N; ++d) {
        if (N % d != 0) continue;
        const int p = d;
        const int q = N / d;
        if ((p + q) & 1) continue;
        const int u = (p + q) / 2;
        const int c_pos = a + b - u;
        if (c_pos > 0) candidates.insert(c_pos);
        const int c_neg = a + b + u;
        if (c_neg > 0) candidates.insert(c_neg);
    }

    u128 total = 0;
    for (int c : candidates) {
        int f = min_steps_bfs(a, b, c, 500, 12);
        total += static_cast<u128>(f);
    }
    return total;
}

u128 compute_F_power(int k,
                     const vector<u128>& pow2,
                     const vector<u128>& pow3,
                     const vector<u128>& pow5) {
    const u128 a = pow2[k] * pow3[k];
    const u128 b = pow2[k] * pow5[k];
    const u128 two_a = a * 2;
    const u128 two_b = b * 2;
    const u128 two_ab = (a + b) * 2;

    const u128 N = pow2[2 * k + 2] * pow3[k] * pow5[k];

    u128 total = 0;

    for (int e2 = 0; e2 <= 2 * k + 2; ++e2) {
        const u128 v2 = pow2[e2];
        for (int e3 = 0; e3 <= k; ++e3) {
            const u128 v23 = v2 * pow3[e3];
            for (int e5 = 0; e5 <= k; ++e5) {
                const u128 q = v23 * pow5[e5];
                const u128 p = N / q;
                if (q > p) continue;
                if ((p + q) & 1) continue;

                // For each divisor pair p*q=4ab, q is the smaller factor.
                // The step count to zero a or b equals the sum of Euclidean
                // quotients for (2a,q) or (2b,q); take the minimum. For the
                // positive-u branch (c <= a+b), we need one fewer step.
                const u128 s1 = euclid_sum(max(two_a, q), min(two_a, q));
                const u128 s2 = euclid_sum(max(two_b, q), min(two_b, q));
                const u128 f = min(s1, s2);

                total += f; // u negative branch is always valid.
                if (p + q < two_ab) {
                    total += f - 1; // u positive branch, if c>0.
                }
            }
        }
    }
    return total;
}

void run_validations(const vector<u128>& pow2,
                     const vector<u128>& pow3,
                     const vector<u128>& pow5) {
    const u128 f1 = compute_F_power(1, pow2, pow3, pow5);
    const u128 f2 = compute_F_power(2, pow2, pow3, pow5);
    if (f1 != 17 || f2 != 179) {
        cerr << "Validation failed for given samples: F(6,10)="
             << to_string_u128(f1) << ", F(36,100)="
             << to_string_u128(f2) << "\n";
        exit(1);
    }

    const vector<pair<int, int>> small_cases = {
        {2, 3},
        {3, 5},
        {4, 6}
    };

    for (const auto& [a, b] : small_cases) {
        const u128 brute = compute_F_bruteforce(a, b);
        // Compute formula result via brute divisor enumeration (small N).
        const int N = 4 * a * b;
        u128 formula = 0;
        for (int d = 1; d * d <= N; ++d) {
            if (N % d != 0) continue;
            const int q = d;
            const int p = N / d;
            if ((p + q) & 1) continue;
            const u128 s1 = euclid_sum(max<u128>(2ULL * a, q), min<u128>(2ULL * a, q));
            const u128 s2 = euclid_sum(max<u128>(2ULL * b, q), min<u128>(2ULL * b, q));
            const u128 f = min(s1, s2);
            formula += f;
            if (p + q < 2 * (a + b)) {
                formula += f - 1;
            }
        }
        if (brute != formula) {
            cerr << "Validation failed for small case (" << a << "," << b
                 << "): brute=" << to_string_u128(brute)
                 << ", formula=" << to_string_u128(formula) << "\n";
            exit(1);
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int maxK = 18;
    unsigned threads = thread::hardware_concurrency();
    bool validate = true;

    // Optional CLI: ./a.out [maxK] [threads] [validate(0/1)]
    if (argc >= 2) maxK = stoi(argv[1]);
    if (argc >= 3) threads = static_cast<unsigned>(stoul(argv[2]));
    if (argc >= 4) validate = (stoi(argv[3]) != 0);

    if (maxK <= 0) {
        cout << "0\n";
        return 0;
    }

    const int max_pow2 = 2 * maxK + 2;
    vector<u128> pow2(max_pow2 + 1, 1);
    for (int i = 1; i <= max_pow2; ++i) pow2[i] = pow2[i - 1] * 2;

    vector<u128> pow3(maxK + 1, 1);
    vector<u128> pow5(maxK + 1, 1);
    for (int i = 1; i <= maxK; ++i) {
        pow3[i] = pow3[i - 1] * 3;
        pow5[i] = pow5[i - 1] * 5;
    }

    if (validate) {
        run_validations(pow2, pow3, pow5);
    }

    if (threads == 0) threads = 1;
    threads = min<unsigned>(threads, static_cast<unsigned>(maxK));

    vector<u128> results(maxK + 1, 0);
    atomic<int> nextK(1);

    auto worker = [&]() {
        while (true) {
            int k = nextK.fetch_add(1);
            if (k > maxK) break;
            results[k] = compute_F_power(k, pow2, pow3, pow5);
        }
    };

    vector<thread> pool;
    pool.reserve(threads);
    for (unsigned i = 0; i < threads; ++i) {
        pool.emplace_back(worker);
    }
    for (auto& th : pool) th.join();

    u128 total = 0;
    for (int k = 1; k <= maxK; ++k) {
        total += results[k];
    }

    cout << to_string_u128(total) << "\n";
    return 0;
}
