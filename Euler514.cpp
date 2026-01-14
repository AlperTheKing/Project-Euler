#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

using namespace std;

struct Direction {
    int dx;
    int dy;
};

struct Point {
    int x;
    int y;
};

static vector<Direction> build_directions(int N) {
    vector<Direction> dirs;
    dirs.reserve(4 * N * N + 4);
    for (int dx = -N; dx <= N; ++dx) {
        for (int dy = -N; dy <= N; ++dy) {
            if (dx == 0 && dy == 0) continue;
            int g = std::gcd(std::abs(dx), std::abs(dy));
            if (g != 1) continue;
            dirs.push_back({dx, dy});
        }
    }
    return dirs;
}

static long double process_direction(int N,
                                     const Direction& dir,
                                     const vector<long double>& pow_q,
                                     long double p,
                                     int total_points) {
    // For a line with direction (dx, dy), each oriented edge contributes via
    // P(edge) = p^2 * q^R * q^O * (1 - q^L), with R/L counts of points on the
    // right/left half-planes and O collinear points outside the segment.
    const int dx = dir.dx;
    const int dy = dir.dy;
    const int n1 = dy;
    const int n2 = -dx;

    auto eval = [&](int x, int y) -> int {
        return n1 * x + n2 * y;
    };

    const int s00 = 0;
    const int sN0 = eval(N, 0);
    const int s0N = eval(0, N);
    const int sNN = eval(N, N);
    const int min_s = min(min(s00, sN0), min(s0N, sNN));
    const int max_s = max(max(s00, sN0), max(s0N, sNN));
    const int range = max_s - min_s + 1;

    vector<int> counts(range, 0);
    vector<Point> starts;
    starts.reserve((N + 1) * (std::abs(dx) + std::abs(dy) + 1));

    for (int x = 0; x <= N; ++x) {
        for (int y = 0; y <= N; ++y) {
            int s = eval(x, y);
            counts[s - min_s]++;
            int px = x - dx;
            int py = y - dy;
            if (px < 0 || px > N || py < 0 || py > N) {
                starts.push_back({x, y});
            }
        }
    }

    vector<int> prefix(range + 1, 0);
    for (int i = 0; i < range; ++i) prefix[i + 1] = prefix[i] + counts[i];

    const long double p2 = p * p;
    long double sum_dir = 0.0L;
    vector<Point> line;

    for (const auto& start : starts) {
        line.clear();
        int x = start.x;
        int y = start.y;
        while (x >= 0 && x <= N && y >= 0 && y <= N) {
            line.push_back({x, y});
            x += dx;
            y += dy;
        }
        int k = static_cast<int>(line.size());
        if (k < 2) continue;

        int s0 = eval(start.x, start.y);
        int idx = s0 - min_s;
        int L = prefix[idx];
        if (L == 0) continue;
        int k_on = counts[idx];
        int R = total_points - prefix[idx] - k_on;

        long double factor = p2 * pow_q[R] * (1.0L - pow_q[L]);
        if (factor == 0.0L) continue;

        long double sum_line = 0.0L;
        for (int i = 0; i < k - 1; ++i) {
            long double wi = pow_q[i];
            int xi = line[i].x;
            int yi = line[i].y;
            for (int j = i + 1; j < k; ++j) {
                long double wj = pow_q[k - 1 - j];
                long double cross = static_cast<long double>(xi) * line[j].y
                                    - static_cast<long double>(line[j].x) * yi;
                sum_line += cross * wi * wj;
            }
        }
        sum_dir += factor * sum_line;
    }
    return sum_dir;
}

static long double compute_expected_area(int N, int threads) {
    const int total_points = (N + 1) * (N + 1);
    const long double p = 1.0L / static_cast<long double>(N + 1);
    const long double q = 1.0L - p;

    vector<long double> pow_q(total_points + 1, 1.0L);
    for (int i = 1; i <= total_points; ++i) pow_q[i] = pow_q[i - 1] * q;

    vector<Direction> dirs = build_directions(N);
    if (threads < 1) threads = 1;
    if (threads > static_cast<int>(dirs.size())) {
        threads = static_cast<int>(dirs.size());
        if (threads < 1) threads = 1;
    }

    vector<long double> partial(threads, 0.0L);
    vector<thread> pool;
    pool.reserve(threads);

    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            long double local = 0.0L;
            for (size_t idx = t; idx < dirs.size(); idx += threads) {
                local += process_direction(N, dirs[idx], pow_q, p, total_points);
            }
            partial[t] = local;
        });
    }

    for (auto& th : pool) th.join();

    long double sum = 0.0L;
    for (long double v : partial) sum += v;
    return 0.5L * sum;
}

static long double round5(long double x) {
    return floor(x * 100000.0L + 0.5L) / 100000.0L;
}

static bool run_validation(int threads) {
    struct Check {
        int N;
        long double expected;
    };

    const Check checks[] = {
        {1, 0.18750L},
        {2, 0.94335L},
        {10, 55.03013L},
    };

    bool ok = true;
    for (const auto& c : checks) {
        long double got = compute_expected_area(c.N, threads);
        long double rounded = round5(got);
        if (fabsl(rounded - c.expected) > 1e-9L) {
            cerr << "Validation failed for N=" << c.N
                 << ": got " << fixed << setprecision(5) << rounded
                 << ", expected " << fixed << setprecision(5) << c.expected
                 << " (raw " << setprecision(10) << got << ")\n";
            ok = false;
        }
    }

    if (ok) {
        cerr << "Validation checkpoints passed." << '\n';
    }
    return ok;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N = 100;
    unsigned hw = thread::hardware_concurrency();
    int threads = hw ? static_cast<int>(hw) : 1;
    threads = max(1, min(threads, 8));
    bool validate = true;

    // Optional CLI: ./a.out [N] [threads] [validate(0/1)]
    if (argc >= 2) N = stoi(argv[1]);
    if (argc >= 3) threads = max(1, stoi(argv[2]));
    if (argc >= 4) validate = (stoi(argv[3]) != 0);

    if (validate && !run_validation(min(threads, 4))) {
        return 1;
    }

    long double answer = compute_expected_area(N, threads);
    cout << fixed << setprecision(5) << answer << "\n";
    return 0;
}
