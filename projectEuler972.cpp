#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using i64 = long long;

struct Fraction {
    i64 num{0};
    i64 den{1};

    Fraction() = default;
    Fraction(i64 n, i64 d = 1) : num(n), den(d) { normalize(); }

    void normalize() {
        if (den == 0) throw std::runtime_error("Zero denominator");
        if (den < 0) {
            num = -num;
            den = -den;
        }
        if (num == 0) {
            den = 1;
            return;
        }
        i64 g = std::gcd(std::abs(num), den);
        num /= g;
        den /= g;
    }

    bool isZero() const { return num == 0; }

    Fraction operator+(const Fraction& other) const {
        i64 n = num * other.den + other.num * den;
        i64 d = den * other.den;
        return Fraction(n, d);
    }
    Fraction operator-(const Fraction& other) const {
        i64 n = num * other.den - other.num * den;
        i64 d = den * other.den;
        return Fraction(n, d);
    }
    Fraction operator*(const Fraction& other) const {
        i64 n = num * other.num;
        i64 d = den * other.den;
        return Fraction(n, d);
    }
    Fraction operator/(const Fraction& other) const {
        if (other.num == 0) throw std::runtime_error("Division by zero");
        i64 n = num * other.den;
        i64 d = den * other.num;
        return Fraction(n, d);
    }

    bool operator==(const Fraction& other) const {
        return num == other.num && den == other.den;
    }
};

struct FractionHash {
    std::size_t operator()(const Fraction& f) const noexcept {
        std::size_t h1 = std::hash<i64>{}(f.num);
        std::size_t h2 = std::hash<i64>{}(f.den);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

int compare(const Fraction& a, const Fraction& b) {
    __int128 lhs = static_cast<__int128>(a.num) * b.den;
    __int128 rhs = static_cast<__int128>(b.num) * a.den;
    if (lhs < rhs) return -1;
    if (lhs > rhs) return 1;
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
        if (a.diameter != b.diameter) return false;
        if (a.vertical != b.vertical) return false;
        if (a.diameter) {
            if (a.vertical) return true;
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
    Fraction det = cross(p, q);
    if (det.isZero()) {
        bool vertical = false;
        Fraction slope(0, 1);
        if (!p.x.isZero()) {
            slope = p.y / p.x;
        } else if (!q.x.isZero()) {
            slope = q.y / q.x;
        } else {
            vertical = true;
        }
        return GeoKey{true, vertical, slope, Fraction(0, 1)};
    }

    Fraction a = (norm2(p) + Fraction(1, 1)) / Fraction(2, 1);
    Fraction b = (norm2(q) + Fraction(1, 1)) / Fraction(2, 1);

    Fraction cx = (a * q.y - p.y * b) / det;
    Fraction cy = (p.x * b - a * q.x) / det;

    return GeoKey{false, false, cx, cy};
}

std::vector<Fraction> generate_rationals(int N) {
    std::vector<Fraction> values;
    values.reserve(N * N);
    std::unordered_set<Fraction, FractionHash> seen;
    seen.reserve(N * N * 2);

    for (int d = 1; d <= N; ++d) {
        for (int n = -d + 1; n <= d - 1; ++n) {
            if (std::gcd(std::abs(n), d) != 1) continue;
            Fraction f(n, d);
            if (seen.insert(f).second) values.push_back(f);
        }
    }

    Fraction zero(0, 1);
    if (seen.insert(zero).second) values.push_back(zero);

    return values;
}

std::vector<Point> build_points(int N) {
    auto rats = generate_rationals(N);
    std::vector<Point> pts;
    pts.reserve(rats.size() * rats.size() / 2);

    Fraction one(1, 1);
    for (const auto& x : rats) {
        for (const auto& y : rats) {
            Point p{x, y};
            if (compare(norm2(p), one) < 0) {
                pts.push_back(p);
            }
        }
    }
    return pts;
}

void process_pairs(size_t start_i, size_t end_i, const std::vector<Point>& pts, GeoMap& local) {
    const size_t n = pts.size();
    for (size_t i = start_i; i < end_i; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            GeoKey key = geodesic_key(pts[i], pts[j]);
            auto& bucket = local[key];
            bucket.insert(static_cast<int>(i));
            bucket.insert(static_cast<int>(j));
        }
    }
}

unsigned long long compute_T(int N, unsigned int thread_count = 0) {
    if (N <= 0) return 0;
    if (thread_count == 0) thread_count = std::thread::hardware_concurrency();
    if (thread_count == 0) thread_count = 4;

    auto pts = build_points(N);
    const size_t n = pts.size();
    if (n < 3) return 0;

    thread_count = std::min<unsigned int>(thread_count, static_cast<unsigned int>(n));

    std::vector<std::thread> workers;
    std::vector<GeoMap> locals(thread_count);

    size_t block = (n + thread_count - 1) / thread_count;
    for (unsigned int t = 0; t < thread_count; ++t) {
        size_t start_i = t * block;
        if (start_i >= n) break;
        size_t end_i = std::min(n, start_i + block);
        workers.emplace_back(process_pairs, start_i, end_i, std::cref(pts), std::ref(locals[t]));
    }

    for (auto& th : workers) th.join();

    GeoMap merged;
    for (auto& mp : locals) {
        for (auto& [key, s] : mp) {
            auto& target = merged[key];
            target.insert(s.begin(), s.end());
        }
    }

    unsigned long long total = 0;
    for (auto& [_, s] : merged) {
        const size_t m = s.size();
        if (m >= 3) {
            total += static_cast<unsigned long long>(m) * (m - 1) * (m - 2);
        }
    }
    return total;
}

bool validate() {
    bool ok = true;
    unsigned long long t2 = compute_T(2);
    std::cout << "Validation T(2) = " << t2;
    if (t2 == 24) std::cout << " [OK]\n";
    else {
        std::cout << " [FAIL]\n";
        ok = false;
    }

    unsigned long long t3 = compute_T(3);
    std::cout << "Validation T(3) = " << t3;
    if (t3 == 1296) std::cout << " [OK]\n";
    else {
        std::cout << " [FAIL]\n";
        ok = false;
    }

    return ok;
}

int main(int argc, char** argv) {
    int N = 12;
    if (argc >= 2) {
        N = std::max(1, std::atoi(argv[1]));
    }

    unsigned int threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 4;

    std::cout << "Running validation checkpoints...\n";
    if (!validate()) {
        std::cerr << "Validation failed; aborting.\n";
        return 1;
    }

    std::cout << "Computing T(" << N << ") with " << threads << " threads...\n";
    unsigned long long answer = compute_T(N, threads);
    std::cout << "T(" << N << ") = " << answer << "\n";
    return 0;
}