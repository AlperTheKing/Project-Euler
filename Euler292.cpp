#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int perimeter_limit = 120;
    bool run_checkpoints = true;
    unsigned requested_threads = 0U;
};

struct Direction {
    int length = 0;
    int dx = 0;
    int dy = 0;
    double angle = 0.0;
};

int gcd2(int x, int y) {
    x = std::abs(x);
    y = std::abs(y);
    if (y == 0) return x;
    if (x == 0) return y;
    while (y != 0) {
        const int r = x % y;
        x = y;
        y = r;
    }
    return x;
}

bool parse_int_after_prefix(const std::string& arg,
                            const std::string& prefix,
                            int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    long long parsed = 0;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10LL + static_cast<long long>(c - '0');
        if (parsed > static_cast<long long>(std::numeric_limits<int>::max())) {
            return false;
        }
    }

    value = static_cast<int>(parsed);
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const std::string& prefix,
                                 unsigned& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
        if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
            return false;
        }
    }

    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--n=", options.perimeter_limit)) {
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", options.requested_threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

std::vector<Direction> init_vectors(const int n) {
    std::vector<Direction> vectors;
    vectors.reserve(static_cast<std::size_t>(4 * n * n + 16));

    for (int y = -n; y <= n; ++y) {
        for (int x = -n; x <= n; ++x) {
            const int t = x * x + y * y;
            const int t2 = static_cast<int>(std::sqrt(static_cast<double>(t)) + 0.5);
            if (gcd2(x, y) != 1 || t2 * t2 != t) {
                continue;
            }

            double ang = std::atan2(static_cast<double>(y), static_cast<double>(x));
            if (ang < 0.0) {
                ang += 6.2831853071795862;
            }
            vectors.push_back(Direction{t2, x, y, ang});
        }
    }

    std::sort(vectors.begin(), vectors.end(),
              [](const Direction& a, const Direction& b) {
                  if (a.angle != b.angle) {
                      return a.angle < b.angle;
                  }
                  if (a.dx != b.dx) {
                      return a.dx < b.dx;
                  }
                  return a.dy < b.dy;
              });

    return vectors;
}

void add_transitions(const std::unordered_map<int, u64>& src,
                     const int dx,
                     const int dy,
                     const int nrl,
                     const int n,
                     std::unordered_map<int, u64>& dst) {
    const int base = n + 1;
    const int nrl2 = nrl * nrl;

    for (const auto& kv : src) {
        const int packed = kv.first;
        const u64 count = kv.second;

        const int cy = packed % base;
        const int cx = packed / base - n;
        const int nx = cx + dx;
        const int ny = cy + dy;
        const int np = nx * nx + ny * ny;

        if (ny >= 0 && np <= nrl2) {
            const int key = (nx + n) * base + ny;
            dst[key] += count;
        }
    }
}

u64 count_polygons(const int n) {
    if (n < 3) {
        return 0ULL;
    }

    const std::vector<Direction> vectors = init_vectors((n - 1) / 2);

    std::vector<std::unordered_map<int, u64>> arr(static_cast<std::size_t>(n + 1));
    for (auto& h : arr) {
        h.reserve(8);
    }

    const int origin_key = n * (n + 1);
    arr[static_cast<std::size_t>(n)][origin_key] = 1ULL;

    for (const Direction& v : vectors) {
        for (int i = v.length; i <= n; ++i) {
            const auto& src = arr[static_cast<std::size_t>(i)];
            if (src.empty()) {
                continue;
            }
            for (int g = 1; g <= i / v.length; ++g) {
                const int nrl = i - v.length * g;
                add_transitions(src, v.dx * g, v.dy * g, nrl, n, arr[static_cast<std::size_t>(nrl)]);
            }
        }
    }

    long long diagonal = 0;
    for (const Direction& v : vectors) {
        diagonal += n / (2 * v.length);
    }
    diagonal /= 2;

    long long answer = -diagonal - 1;
    for (int i = 0; i <= n; ++i) {
        const auto it = arr[static_cast<std::size_t>(i)].find(origin_key);
        if (it != arr[static_cast<std::size_t>(i)].end()) {
            answer += static_cast<long long>(it->second);
        }
    }

    return static_cast<u64>(answer);
}

void run_checkpoints() {
    struct Checkpoint {
        int n;
        u64 expected;
    };

    const std::vector<Checkpoint> checkpoints = {
        {4, 1ULL},
        {30, 3655ULL},
        {60, 891045ULL},
    };

    for (const Checkpoint& cp : checkpoints) {
        const u64 got = count_polygons(cp.n);
        if (got != cp.expected) {
            throw std::runtime_error("Checkpoint failed for P(" +
                                     std::to_string(cp.n) + "): got " +
                                     std::to_string(got) + ", expected " +
                                     std::to_string(cp.expected));
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        Options options;
        if (!parse_arguments(argc, argv, options)) {
            return 1;
        }

        if (options.perimeter_limit < 0) {
            throw std::invalid_argument("--n must be non-negative");
        }

        if (options.run_checkpoints) {
            run_checkpoints();
        }

        const u64 answer = count_polygons(options.perimeter_limit);
        std::cout << answer << '\n';
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
