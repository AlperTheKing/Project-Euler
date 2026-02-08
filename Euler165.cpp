#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace {

using i64 = long long;

struct Options {
    int segments = 5000;
    bool run_checkpoints = true;
};

struct Point {
    i64 x = 0;
    i64 y = 0;
};

struct Segment {
    Point a;
    Point b;
};

i64 cross(const Point& a, const Point& b, const Point& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

int sgn(const i64 x) {
    return (x > 0) - (x < 0);
}

struct RationalPoint {
    i64 x_num = 0;
    i64 x_den = 1;
    i64 y_num = 0;
    i64 y_den = 1;

    bool operator<(const RationalPoint& other) const {
        return std::tie(x_num, x_den, y_num, y_den) < std::tie(other.x_num, other.x_den, other.y_num, other.y_den);
    }
};

RationalPoint intersection_point(const Segment& s1, const Segment& s2) {
    const i64 a1 = s1.b.y - s1.a.y;
    const i64 b1 = s1.a.x - s1.b.x;
    const i64 c1 = a1 * s1.a.x + b1 * s1.a.y;

    const i64 a2 = s2.b.y - s2.a.y;
    const i64 b2 = s2.a.x - s2.b.x;
    const i64 c2 = a2 * s2.a.x + b2 * s2.a.y;

    i64 den = a1 * b2 - a2 * b1;
    i64 x_num = c1 * b2 - c2 * b1;
    i64 y_num = a1 * c2 - a2 * c1;

    if (den < 0) {
        den = -den;
        x_num = -x_num;
        y_num = -y_num;
    }

    const i64 gx = std::gcd(std::llabs(x_num), den);
    const i64 gy = std::gcd(std::llabs(y_num), den);

    RationalPoint out;
    out.x_num = x_num / gx;
    out.x_den = den / gx;
    out.y_num = y_num / gy;
    out.y_den = den / gy;
    return out;
}

std::vector<Segment> generate_segments(const int n) {
    std::vector<i64> s(static_cast<std::size_t>(4 * n + 1), 0);
    s[0] = 290797;
    for (int i = 0; i < 4 * n; ++i) {
        s[static_cast<std::size_t>(i + 1)] = (s[static_cast<std::size_t>(i)] * s[static_cast<std::size_t>(i)]) % 50515093;
    }

    std::vector<Segment> segs(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        segs[static_cast<std::size_t>(i)].a.x = s[static_cast<std::size_t>(4 * i + 1)] % 500;
        segs[static_cast<std::size_t>(i)].a.y = s[static_cast<std::size_t>(4 * i + 2)] % 500;
        segs[static_cast<std::size_t>(i)].b.x = s[static_cast<std::size_t>(4 * i + 3)] % 500;
        segs[static_cast<std::size_t>(i)].b.y = s[static_cast<std::size_t>(4 * i + 4)] % 500;
    }
    return segs;
}

i64 count_true_intersections(const std::vector<Segment>& segs) {
    std::set<RationalPoint> points;
    const int segments = static_cast<int>(segs.size());

    for (int i = 0; i < segments; ++i) {
        for (int j = i + 1; j < segments; ++j) {
            const i64 c1 = cross(segs[static_cast<std::size_t>(i)].a, segs[static_cast<std::size_t>(i)].b, segs[static_cast<std::size_t>(j)].a);
            const i64 c2 = cross(segs[static_cast<std::size_t>(i)].a, segs[static_cast<std::size_t>(i)].b, segs[static_cast<std::size_t>(j)].b);
            const i64 c3 = cross(segs[static_cast<std::size_t>(j)].a, segs[static_cast<std::size_t>(j)].b, segs[static_cast<std::size_t>(i)].a);
            const i64 c4 = cross(segs[static_cast<std::size_t>(j)].a, segs[static_cast<std::size_t>(j)].b, segs[static_cast<std::size_t>(i)].b);

            if (sgn(c1) * sgn(c2) >= 0 || sgn(c3) * sgn(c4) >= 0) {
                continue;
            }

            points.insert(intersection_point(segs[static_cast<std::size_t>(i)], segs[static_cast<std::size_t>(j)]));
        }
    }

    return static_cast<i64>(points.size());
}

i64 solve(const int segments) {
    const std::vector<Segment> segs = generate_segments(segments);
    return count_true_intersections(segs);
}

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
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
        if (parse_int_after_prefix(arg, "--segments=", options.segments)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.segments >= 2;
}

bool run_checkpoints() {
    const std::vector<Segment> sample = generate_segments(1);
    if (sample[0].a.x != 27 || sample[0].a.y != 144 || sample[0].b.x != 12 || sample[0].b.y != 232) {
        std::cerr << "Checkpoint failed for BBS seed" << '\n';
        return false;
    }

    std::vector<Segment> stated(3);
    stated[0] = Segment{Point{27, 44}, Point{12, 32}};
    stated[1] = Segment{Point{46, 53}, Point{17, 62}};
    stated[2] = Segment{Point{46, 70}, Point{22, 40}};
    if (count_true_intersections(stated) != 1LL) {
        std::cerr << "Checkpoint failed for stated sample segments" << '\n';
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

    std::cout << solve(options.segments) << '\n';
    return 0;
}
