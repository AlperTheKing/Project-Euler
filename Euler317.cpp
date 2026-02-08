#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct Options {
    long double height = 100.0L;
    long double speed = 20.0L;
    long double gravity = 9.81L;
    bool run_checkpoints = true;
};

bool parse_ld_after_prefix(const std::string& arg, const std::string& prefix, long double& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    char* end_ptr = nullptr;
    value = std::strtold(tail.c_str(), &end_ptr);
    return end_ptr != nullptr && *end_ptr == '\0';
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_ld_after_prefix(arg, "--height=", options.height) ||
            parse_ld_after_prefix(arg, "--speed=", options.speed) ||
            parse_ld_after_prefix(arg, "--gravity=", options.gravity)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.height >= 0.0L && options.speed > 0.0L && options.gravity > 0.0L;
}

long double solve(const long double h, const long double v, const long double g) {
    // Envelope in cylindrical coordinates:
    // y(r) = h + v^2/(2g) - g*r^2/(2v^2), for 0 <= r <= r_max where y=0.
    const long double a = h + (v * v) / (2.0L * g);
    const long double b = g / (2.0L * v * v);
    const long double r_max2 = a / b;

    // Volume of revolution: V = pi * integral_0^{r_max} r^2 * 2? No:
    // V = pi * integral_0^{y_max} r(y)^2 dy = pi * integral_0^{r_max} 2r y(r) dr.
    // Directly with y(r)=a-br^2:
    // V = pi * integral_0^{r_max} (a r^2 - b r^4)'? easiest via shells:
    // V = 2*pi * integral_0^{r_max} r*(a-br^2) dr = pi * a^2 / (2b).
    const long double pi = std::acos(-1.0L);
    return pi * a * a / (2.0L * b);
}

bool run_checkpoints() {
    // For h=0, the envelope is centered at the launch point and the volume becomes:
    // pi * (v^6) / (4 g^3).
    {
        const long double v = 10.0L;
        const long double g = 5.0L;
        const long double expected =
            std::acos(-1.0L) * v * v * v * v * v * v / (4.0L * g * g * g);
        const long double got = solve(0.0L, v, g);
        if (std::fabsl(got - expected) > 1e-12L) {
            std::cerr << "Checkpoint failed for h=0 closed form" << '\n';
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
    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }
    const long double ans = solve(options.height, options.speed, options.gravity);
    std::cout << std::fixed << std::setprecision(4) << static_cast<double>(ans) << '\n';
    return 0;
}
