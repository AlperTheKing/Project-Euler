#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using i64 = long long;
using i128 = __int128_t;
using u64 = std::uint64_t;

struct Options {
    int n = 12;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    int parsed = 0;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1;
}

struct Point {
    i64 x;
    i64 y;
    i64 z;
};

struct Key {
    i64 u;
    i64 v;
};

bool operator<(const Key& a, const Key& b) {
    if (a.u != b.u) {
        return a.u < b.u;
    }
    return a.v < b.v;
}

inline i64 abs64(i64 x) { return x < 0 ? -x : x; }

u64 compute_T(int n) {
    i64 L = 1;
    for (int d = 1; d <= n; ++d) {
        L = std::lcm(L, static_cast<i64>(d));
    }

    std::vector<i64> coords;
    coords.reserve(static_cast<std::size_t>(2 * n * n + 1));
    for (int d = 1; d <= n; ++d) {
        const i64 q = L / d;
        for (int num = -d + 1; num <= d - 1; ++num) {
            coords.push_back(static_cast<i64>(num) * q);
        }
    }
    std::sort(coords.begin(), coords.end());
    coords.erase(std::unique(coords.begin(), coords.end()), coords.end());

    const i64 L2 = L * L;
    std::vector<Point> pts;
    pts.reserve(coords.size() * coords.size() / 2U);

    for (i64 x : coords) {
        const i64 x2 = x * x;
        for (i64 y : coords) {
            const i64 y2 = y * y;
            if (x2 + y2 >= L2) {
                continue;
            }
            pts.push_back(Point{x * L, y * L, x2 + y2 + L2});
        }
    }

    u64 count = 0;
    std::vector<Key> keys;

    for (std::size_t i = 0; i + 2 < pts.size(); ++i) {
        const Point& p = pts[i];
        keys.clear();
        keys.reserve(pts.size() - i - 1);

        for (std::size_t j = i + 1; j < pts.size(); ++j) {
            const Point& q = pts[j];
            i128 u;
            i128 v;

            if (p.x == 0) {
                u = q.x;
                v = static_cast<i128>(p.z) * q.y - static_cast<i128>(p.y) * q.z;
            } else {
                u = static_cast<i128>(p.y) * q.x - static_cast<i128>(p.x) * q.y;
                v = static_cast<i128>(p.z) * q.x - static_cast<i128>(p.x) * q.z;
            }

            i64 uu = static_cast<i64>(u);
            i64 vv = static_cast<i64>(v);
            i64 g = std::gcd(abs64(uu), abs64(vv));
            if (g == 0) {
                g = 1;
            }

            uu /= g;
            vv /= g;

            if (uu < 0 || (uu == 0 && vv < 0)) {
                uu = -uu;
                vv = -vv;
            }

            keys.push_back(Key{uu, vv});
        }

        std::sort(keys.begin(), keys.end());
        for (std::size_t s = 0; s < keys.size();) {
            std::size_t t = s + 1;
            while (t < keys.size() && keys[t].u == keys[s].u && keys[t].v == keys[s].v) {
                ++t;
            }
            const u64 cnt = static_cast<u64>(t - s);
            count += cnt * (cnt - 1) / 2;
            s = t;
        }
    }

    return 6 * count;
}

bool run_checkpoints() {
    if (compute_T(2) != 24ULL) {
        std::cerr << "Checkpoint failed: T(2)=24" << '\n';
        return false;
    }
    if (compute_T(3) != 1296ULL) {
        std::cerr << "Checkpoint failed: T(3)=1296" << '\n';
        return false;
    }
    if (compute_T(4) != 5052ULL) {
        std::cerr << "Checkpoint failed: T(4)=5052" << '\n';
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }
    std::cout << compute_T(options.n) << '\n';
    return 0;
}
