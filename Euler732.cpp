#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;

constexpr i64 kMod = 1'000'000'007LL;

struct Job {
    int deadline;
    int p;
    int q;
};

int ceil_div_by_sqrt2(const int hsum) {
    const i64 target = static_cast<i64>(hsum) * static_cast<i64>(hsum);
    int lo = 0;
    int hi = hsum;
    while (lo < hi) {
        const int mid = lo + (hi - lo) / 2;
        const i64 lhs = 2LL * static_cast<i64>(mid) * static_cast<i64>(mid);
        if (lhs >= target) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    return lo;
}

int Q(const int n) {
    std::vector<int> h(static_cast<std::size_t>(n));
    std::vector<int> l(static_cast<std::size_t>(n));
    std::vector<int> q(static_cast<std::size_t>(n));

    i64 pow5 = 1;
    for (int idx = 0; idx < 3 * n; ++idx) {
        if (idx > 0) {
            pow5 = (pow5 * 5LL) % kMod;
        }
        const int r = static_cast<int>(pow5 % 101LL) + 50;
        const int troll = idx / 3;
        if (idx % 3 == 0) {
            h[static_cast<std::size_t>(troll)] = r;
        } else if (idx % 3 == 1) {
            l[static_cast<std::size_t>(troll)] = r;
        } else {
            q[static_cast<std::size_t>(troll)] = r;
        }
    }

    int hsum = 0;
    for (const int x : h) {
        hsum += x;
    }

    const int ceil_h_over_sqrt2 = ceil_div_by_sqrt2(hsum);

    std::vector<Job> jobs;
    jobs.reserve(static_cast<std::size_t>(n));

    int max_deadline = 0;
    for (int i = 0; i < n; ++i) {
        const int start_deadline = hsum + l[static_cast<std::size_t>(i)] - ceil_h_over_sqrt2;
        const int completion_deadline = start_deadline + h[static_cast<std::size_t>(i)];
        jobs.push_back(Job{completion_deadline, h[static_cast<std::size_t>(i)], q[static_cast<std::size_t>(i)]});
        if (completion_deadline > max_deadline) {
            max_deadline = completion_deadline;
        }
    }

    std::sort(jobs.begin(), jobs.end(), [](const Job& a, const Job& b) {
        if (a.deadline != b.deadline) {
            return a.deadline < b.deadline;
        }
        return a.p < b.p;
    });

    constexpr int kNeg = -1'000'000'000;
    std::vector<int> dp(static_cast<std::size_t>(max_deadline + 1), kNeg);
    dp[0] = 0;

    for (const Job& job : jobs) {
        for (int t = job.deadline; t >= job.p; --t) {
            const int prev = dp[static_cast<std::size_t>(t - job.p)];
            if (prev == kNeg) {
                continue;
            }
            const int cand = prev + job.q;
            if (cand > dp[static_cast<std::size_t>(t)]) {
                dp[static_cast<std::size_t>(t)] = cand;
            }
        }
    }

    int best = 0;
    for (const int v : dp) {
        if (v > best) {
            best = v;
        }
    }
    return best;
}

}  // namespace

int main() {
    assert(Q(5) == 401);
    assert(Q(15) == 941);

    std::cout << Q(1000) << '\n';
    return 0;
}
