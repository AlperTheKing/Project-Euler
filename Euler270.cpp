#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = uint32_t;
using u64 = uint64_t;

constexpr u32 kMod = 100000000u;
constexpr int kDefaultN = 30;

struct Point {
    int x = 0;
    int y = 0;
};

std::vector<Point> build_boundary(int n) {
    std::vector<Point> pts;
    pts.reserve(4 * n);
    for (int x = 0; x <= n; ++x) pts.push_back({x, 0});
    for (int y = 1; y <= n; ++y) pts.push_back({n, y});
    for (int x = n - 1; x >= 0; --x) pts.push_back({x, n});
    for (int y = n - 1; y >= 1; --y) pts.push_back({0, y});
    return pts;
}

bool same_side(const Point& a, const Point& b, int n) {
    if (a.y == 0 && b.y == 0) return true;
    if (a.y == n && b.y == n) return true;
    if (a.x == 0 && b.x == 0) return true;
    if (a.x == n && b.x == n) return true;
    return false;
}

u32 count_triangulations(int n, int threads) {
    if (n <= 0) return 0;
    std::vector<Point> pts = build_boundary(n);
    const int m = static_cast<int>(pts.size());

    // Maximal non-crossing cuts correspond to triangulations using only allowed diagonals.
    std::vector<std::vector<uint8_t>> allowed(m, std::vector<uint8_t>(m, 0));
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            bool ok = false;
            if (j == i + 1 || (i == 0 && j == m - 1)) {
                ok = true;
            } else if (!same_side(pts[i], pts[j], n)) {
                ok = true;
            }
            allowed[i][j] = allowed[j][i] = static_cast<uint8_t>(ok);
        }
    }

    std::vector<std::vector<u32>> dp(m, std::vector<u32>(m, 0));
    for (int i = 0; i + 1 < m; ++i) dp[i][i + 1] = 1;

    auto compute_range = [&](int length, int start, int end) {
        for (int idx = start; idx < end; ++idx) {
            int i = idx;
            int j = i + length;
            if (!allowed[i][j]) {
                dp[i][j] = 0;
                continue;
            }
            u64 sum = 0;
            for (int k = i + 1; k < j; ++k) {
                u32 left = dp[i][k];
                if (left == 0) continue;
                u32 right = dp[k][j];
                if (right == 0) continue;
                sum += static_cast<u64>(left) * right;
                if (sum >= static_cast<u64>(kMod) * kMod) sum %= kMod;
            }
            dp[i][j] = static_cast<u32>(sum % kMod);
        }
    };

    for (int length = 2; length <= m - 1; ++length) {
        int tasks = m - length;
        int use_threads = threads;
        if (use_threads < 2 || tasks < 32) {
            compute_range(length, 0, tasks);
            continue;
        }
        use_threads = std::min(use_threads, tasks);
        int chunk = (tasks + use_threads - 1) / use_threads;
        std::vector<std::thread> pool;
        pool.reserve(use_threads);
        for (int t = 0; t < use_threads; ++t) {
            int start = t * chunk;
            int end = std::min(start + chunk, tasks);
            if (start >= end) break;
            pool.emplace_back([&, length, start, end]() {
                compute_range(length, start, end);
            });
        }
        for (auto& th : pool) th.join();
    }

    return dp[0][m - 1];
}

void run_validation() {
    struct Check {
        int n;
        u32 expected;
    };
    const Check checks[] = {
        {1, 2},
        {2, 30},
        {3, 604},
        {4, 12168},
    };
    for (const auto& chk : checks) {
        u32 got = count_triangulations(chk.n, 1);
        if (got != chk.expected) {
            std::cerr << "Validation failed for N=" << chk.n << ": got " << got
                      << ", expected " << chk.expected << '\n';
            std::exit(1);
        }
    }
}

void usage(const char* argv0) {
    std::cerr << "Usage: " << argv0 << " [-n N] [-t THREADS] [--no-validate]\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int n = kDefaultN;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads <= 0) threads = 1;
    bool validate = true;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-n" && i + 1 < argc) {
            n = std::stoi(argv[++i]);
        } else if (arg == "-t" && i + 1 < argc) {
            threads = std::max(1, std::stoi(argv[++i]));
        } else if (arg == "--no-validate") {
            validate = false;
        } else if (arg == "--validate") {
            validate = true;
        } else {
            usage(argv[0]);
            return 1;
        }
    }

    if (validate) run_validation();

    u32 answer = count_triangulations(n, threads);
    std::cout << answer << '\n';
    return 0;
}
