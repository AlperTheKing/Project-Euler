#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

struct Options {
    int n = 12;
    int threads = 0;  // 0 => auto
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
        if (parse_int_after_prefix(arg, "--n=", options.n) ||
            parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1 && options.threads >= 0;
}

struct Fraction {
    i64 num{0};
    i64 den{1};

    Fraction() = default;
    Fraction(i64 n, i64 d = 1) : num(n), den(d) {
        normalize();
    }

    void normalize() {
        if (den == 0) {
            throw std::runtime_error("Zero denominator");
        }
        if (den < 0) {
            num = -num;
            den = -den;
        }
        if (num == 0) {
            den = 1;
            return;
        }
        const i64 g = std::gcd(std::llabs(num), den);
        num /= g;
        den /= g;
    }

    bool is_zero() const {
        return num == 0;
    }

    Fraction operator+(const Fraction& other) const {
        return Fraction(num * other.den + other.num * den, den * other.den);
    }
    Fraction operator-(const Fraction& other) const {
        return Fraction(num * other.den - other.num * den, den * other.den);
    }
    Fraction operator*(const Fraction& other) const {
        return Fraction(num * other.num, den * other.den);
    }
    Fraction operator/(const Fraction& other) const {
        if (other.num == 0) {
            throw std::runtime_error("Division by zero");
        }
        return Fraction(num * other.den, den * other.num);
    }

    bool operator==(const Fraction& other) const {
        return num == other.num && den == other.den;
    }
};

struct FractionHash {
    std::size_t operator()(const Fraction& f) const noexcept {
        const std::size_t h1 = std::hash<i64>{}(f.num);
        const std::size_t h2 = std::hash<i64>{}(f.den);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

int compare_fraction(const Fraction& a, const Fraction& b) {
    const __int128 lhs = static_cast<__int128>(a.num) * b.den;
    const __int128 rhs = static_cast<__int128>(b.num) * a.den;
    if (lhs < rhs) {
        return -1;
    }
    if (lhs > rhs) {
        return 1;
    }
    return 0;
}

struct Point {
    Fraction x;
    Fraction y;
};

struct GeoKey {
    bool diameter{false};
    bool vertical{false};
    Fraction value1{};
    Fraction value2{};
};

struct GeoKeyHash {
    std::size_t operator()(const GeoKey& k) const noexcept {
        std::size_t h = std::hash<bool>{}(k.diameter);
        h ^= (std::hash<bool>{}(k.vertical) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
        FractionHash fh;
        if (k.diameter) {
            if (!k.vertical) {
                h ^= (fh(k.value1) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
            }
        } else {
            h ^= (fh(k.value1) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
            h ^= (fh(k.value2) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
        }
        return h;
    }
};

struct GeoKeyEq {
    bool operator()(const GeoKey& a, const GeoKey& b) const noexcept {
        if (a.diameter != b.diameter || a.vertical != b.vertical) {
            return false;
        }
        if (a.diameter) {
            if (a.vertical) {
                return true;
            }
            return a.value1 == b.value1;
        }
        return a.value1 == b.value1 && a.value2 == b.value2;
    }
};

using PointSet = std::unordered_set<int>;
using GeoMap = std::unordered_map<GeoKey, PointSet, GeoKeyHash, GeoKeyEq>;

Fraction norm2(const Point& p) {
    return p.x * p.x + p.y * p.y;
}

Fraction cross(const Point& a, const Point& b) {
    return a.x * b.y - a.y * b.x;
}

GeoKey geodesic_key(const Point& p, const Point& q) {
    const Fraction det = cross(p, q);
    if (det.is_zero()) {
        bool vertical = false;
        Fraction slope(0, 1);
        if (!p.x.is_zero()) {
            slope = p.y / p.x;
        } else if (!q.x.is_zero()) {
            slope = q.y / q.x;
        } else {
            vertical = true;
        }
        return GeoKey{true, vertical, slope, Fraction(0, 1)};
    }

    const Fraction a = (norm2(p) + Fraction(1, 1)) / Fraction(2, 1);
    const Fraction b = (norm2(q) + Fraction(1, 1)) / Fraction(2, 1);
    const Fraction cx = (a * q.y - p.y * b) / det;
    const Fraction cy = (p.x * b - a * q.x) / det;
    return GeoKey{false, false, cx, cy};
}

std::vector<Fraction> generate_rationals(const int n) {
    std::vector<Fraction> values;
    values.reserve(static_cast<std::size_t>(n * n));
    std::unordered_set<Fraction, FractionHash> seen;
    seen.reserve(static_cast<std::size_t>(n * n * 2));

    for (int d = 1; d <= n; ++d) {
        for (int x = -d + 1; x <= d - 1; ++x) {
            if (std::gcd(std::abs(x), d) != 1) {
                continue;
            }
            Fraction f(x, d);
            if (seen.insert(f).second) {
                values.push_back(f);
            }
        }
    }
    Fraction zero(0, 1);
    if (seen.insert(zero).second) {
        values.push_back(zero);
    }
    return values;
}

std::vector<Point> build_points(const int n) {
    const std::vector<Fraction> rats = generate_rationals(n);
    std::vector<Point> points;
    points.reserve(rats.size() * rats.size() / 2U);

    const Fraction one(1, 1);
    for (const Fraction& x : rats) {
        for (const Fraction& y : rats) {
            Point p{x, y};
            if (compare_fraction(norm2(p), one) < 0) {
                points.push_back(p);
            }
        }
    }
    return points;
}

void process_pairs(const std::size_t start_i, const std::size_t end_i, const std::vector<Point>& pts, GeoMap& local) {
    const std::size_t n = pts.size();
    for (std::size_t i = start_i; i < end_i; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const GeoKey key = geodesic_key(pts[i], pts[j]);
            auto& bucket = local[key];
            bucket.insert(static_cast<int>(i));
            bucket.insert(static_cast<int>(j));
        }
    }
}

u64 compute_T(const int n, int threads) {
    if (n <= 0) {
        return 0ULL;
    }
    if (threads <= 0) {
        threads = static_cast<int>(std::max(1U, std::thread::hardware_concurrency()));
    }

    const std::vector<Point> points = build_points(n);
    const std::size_t count = points.size();
    if (count < 3) {
        return 0ULL;
    }

    threads = std::min<int>(threads, static_cast<int>(count));
    const std::size_t block = (count + static_cast<std::size_t>(threads) - 1U) / static_cast<std::size_t>(threads);

    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(threads));
    std::vector<GeoMap> locals(static_cast<std::size_t>(threads));

    for (int t = 0; t < threads; ++t) {
        const std::size_t start = static_cast<std::size_t>(t) * block;
        if (start >= count) {
            break;
        }
        const std::size_t end = std::min(count, start + block);
        workers.emplace_back(process_pairs, start, end, std::cref(points), std::ref(locals[static_cast<std::size_t>(t)]));
    }
    for (std::thread& th : workers) {
        th.join();
    }

    GeoMap merged;
    for (GeoMap& mp : locals) {
        for (auto& row : mp) {
            auto& target = merged[row.first];
            target.insert(row.second.begin(), row.second.end());
        }
    }

    u64 total = 0ULL;
    for (const auto& row : merged) {
        const std::size_t m = row.second.size();
        if (m >= 3U) {
            total += static_cast<u64>(m) * static_cast<u64>(m - 1U) * static_cast<u64>(m - 2U);
        }
    }
    return total;
}

bool run_checkpoints(const int threads) {
    if (compute_T(2, 1) != 24ULL) {
        std::cerr << "Checkpoint failed: T(2)=24" << '\n';
        return false;
    }
    if (compute_T(3, 1) != 1296ULL) {
        std::cerr << "Checkpoint failed: T(3)=1296" << '\n';
        return false;
    }
    if (threads > 1) {
        if (compute_T(4, 1) != compute_T(4, threads)) {
            std::cerr << "Checkpoint failed: serial/parallel determinism" << '\n';
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    const int threads = (options.threads > 0)
                            ? options.threads
                            : static_cast<int>(std::max(1U, std::thread::hardware_concurrency()));
    if (options.run_checkpoints && !run_checkpoints(threads)) {
        return 2;
    }
    std::cout << compute_T(options.n, threads) << '\n';
    return 0;
}
