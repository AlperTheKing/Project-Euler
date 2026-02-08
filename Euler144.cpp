#include <cmath>
#include <iostream>
#include <string>

namespace {

struct Options {
    bool run_checkpoints = true;
};

bool nearly_equal(const double a, const double b, const double eps = 1e-9) {
    return std::abs(a - b) <= eps;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

struct Point {
    double x = 0.0;
    double y = 0.0;
};

Point reflect_direction(const Point& incoming, const Point& point) {
    const double nx = 8.0 * point.x;
    const double ny = 2.0 * point.y;

    const double dot = incoming.x * nx + incoming.y * ny;
    const double norm2 = nx * nx + ny * ny;

    return {
        incoming.x - 2.0 * dot * nx / norm2,
        incoming.y - 2.0 * dot * ny / norm2,
    };
}

Point next_intersection(const Point& current, const Point& direction) {
    const double num = -(8.0 * current.x * direction.x + 2.0 * current.y * direction.y);
    const double den = 4.0 * direction.x * direction.x + direction.y * direction.y;
    const double t = num / den;

    return {
        current.x + t * direction.x,
        current.y + t * direction.y,
    };
}

int solve() {
    Point previous{0.0, 10.1};
    Point current{1.4, -9.6};

    int reflections = 0;
    while (true) {
        const Point incoming{current.x - previous.x, current.y - previous.y};
        const Point outgoing = reflect_direction(incoming, current);
        const Point next = next_intersection(current, outgoing);

        ++reflections;
        previous = current;
        current = next;

        if (std::abs(current.x) <= 0.01 && current.y > 0.0) {
            break;
        }
    }

    return reflections;
}

bool run_checkpoints() {
    const Point previous{0.0, 10.1};
    const Point current{1.4, -9.6};
    const Point incoming{current.x - previous.x, current.y - previous.y};
    const Point outgoing = reflect_direction(incoming, current);
    const Point next = next_intersection(current, outgoing);

    if (!nearly_equal(next.x, -3.990597619361616, 1e-9) ||
        !nearly_equal(next.y, -6.024991498863843, 1e-9)) {
        std::cerr << "Checkpoint failed for first reflected intersection" << '\n';
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

    std::cout << solve() << '\n';
    return 0;
}
