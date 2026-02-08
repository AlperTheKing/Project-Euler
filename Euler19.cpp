#include <iostream>
#include <string>

namespace {

struct Options {
    int start_year = 1901;
    int end_year = 2000;
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
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }

    value = parsed;
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--start-year=", options.start_year)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--end-year=", options.end_year)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.start_year >= 1900 && options.end_year >= options.start_year;
}

bool is_leap_year(const int year) {
    if (year % 400 == 0) {
        return true;
    }
    if (year % 100 == 0) {
        return false;
    }
    return (year % 4 == 0);
}

int days_in_month(const int year, const int month) {
    static const int month_days[] = {
        0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
    };
    if (month == 2 && is_leap_year(year)) {
        return 29;
    }
    return month_days[month];
}

int weekday_monday_zero(const int year, const int month, const int day) {
    // 1900-01-01 is Monday -> 0.
    int days = 0;
    for (int y = 1900; y < year; ++y) {
        days += is_leap_year(y) ? 366 : 365;
    }
    for (int m = 1; m < month; ++m) {
        days += days_in_month(year, m);
    }
    days += day - 1;

    return days % 7;
}

int solve(const int start_year, const int end_year) {
    int count = 0;
    for (int year = start_year; year <= end_year; ++year) {
        for (int month = 1; month <= 12; ++month) {
            const int weekday = weekday_monday_zero(year, month, 1);
            if (weekday == 6) {  // Sunday
                ++count;
            }
        }
    }
    return count;
}

bool run_checkpoints() {
    if (weekday_monday_zero(1900, 1, 1) != 0) {
        std::cerr << "Checkpoint failed for 1900-01-01 weekday" << '\n';
        return false;
    }
    if (weekday_monday_zero(1901, 1, 1) != 1) {
        std::cerr << "Checkpoint failed for 1901-01-01 weekday" << '\n';
        return false;
    }
    if (solve(1901, 1901) != 2) {
        std::cerr << "Checkpoint failed for year 1901 count" << '\n';
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

    std::cout << solve(options.start_year, options.end_year) << '\n';
    return 0;
}
