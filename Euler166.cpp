#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = long long;

struct Options {
    int max_digit = 9;
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
        if (parse_int_after_prefix(arg, "--max-digit=", options.max_digit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.max_digit >= 0 && options.max_digit <= 9;
}

i64 solve(const int dmax) {
    i64 answer = 0;

    for (int a = 0; a <= dmax; ++a) {
        for (int b = 0; b <= dmax; ++b) {
            for (int c = 0; c <= dmax; ++c) {
                for (int g = 0; g <= dmax; ++g) {
                    for (int j = 0; j <= dmax; ++j) {
                        const int m = a + b + c - g - j;
                        if (m < 0 || m > dmax) {
                            continue;
                        }

                        for (int e = 0; e <= dmax; ++e) {
                            for (int i = 0; i <= dmax; ++i) {
                                const int S = a + e + i + m;

                                const int d = S - a - b - c;
                                if (d < 0 || d > dmax) {
                                    continue;
                                }

                                for (int f = 0; f <= dmax; ++f) {
                                    const int h = S - e - f - g;
                                    if (h < 0 || h > dmax) {
                                        continue;
                                    }

                                    const int n = S - b - f - j;
                                    if (n < 0 || n > dmax) {
                                        continue;
                                    }

                                    for (int k = 0; k <= dmax; ++k) {
                                        const int o = S - c - g - k;
                                        if (o < 0 || o > dmax) {
                                            continue;
                                        }

                                        const int l = S - i - j - k;
                                        if (l < 0 || l > dmax) {
                                            continue;
                                        }

                                        const int p = S - a - f - k;
                                        if (p < 0 || p > dmax) {
                                            continue;
                                        }

                                        if (m + n + o + p != S) {
                                            continue;
                                        }
                                        if (d + h + l + p != S) {
                                            continue;
                                        }

                                        ++answer;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return answer;
}

i64 brute_small(const int dmax) {
    i64 count = 0;
    for (int a = 0; a <= dmax; ++a)
    for (int b = 0; b <= dmax; ++b)
    for (int c = 0; c <= dmax; ++c)
    for (int d = 0; d <= dmax; ++d)
    for (int e = 0; e <= dmax; ++e)
    for (int f = 0; f <= dmax; ++f)
    for (int g = 0; g <= dmax; ++g)
    for (int h = 0; h <= dmax; ++h)
    for (int i = 0; i <= dmax; ++i)
    for (int j = 0; j <= dmax; ++j)
    for (int k = 0; k <= dmax; ++k)
    for (int l = 0; l <= dmax; ++l)
    for (int m = 0; m <= dmax; ++m)
    for (int n = 0; n <= dmax; ++n)
    for (int o = 0; o <= dmax; ++o)
    for (int p = 0; p <= dmax; ++p) {
        const int S = a + b + c + d;
        if (e + f + g + h != S) continue;
        if (i + j + k + l != S) continue;
        if (m + n + o + p != S) continue;
        if (a + e + i + m != S) continue;
        if (b + f + j + n != S) continue;
        if (c + g + k + o != S) continue;
        if (d + h + l + p != S) continue;
        if (a + f + k + p != S) continue;
        if (d + g + j + m != S) continue;
        ++count;
    }
    return count;
}

bool run_checkpoints() {
    if (solve(1) != 34LL) {
        std::cerr << "Checkpoint failed for max digit 1" << '\n';
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

    std::cout << solve(options.max_digit) << '\n';
    return 0;
}
