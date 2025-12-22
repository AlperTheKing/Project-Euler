#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace std;

static constexpr int MAX_STEPS = 200;
static constexpr long double STOP_DEPTH = 80.0L;

static bool is_increasing(long double d1, long double *min_diff_out = nullptr) {
    long double prev = 0.0L;
    long double curr = d1;
    long double min_diff = numeric_limits<long double>::infinity();
    for (int i = 0; i < MAX_STEPS; ++i) {
        long double diff = curr - prev;
        if (diff <= 0.0L) {
            if (min_diff_out) {
                *min_diff_out = diff;
            }
            return false;
        }
        if (diff < min_diff) {
            min_diff = diff;
        }
        if (curr > STOP_DEPTH) {
            break;
        }
        long double next = expl(curr - prev);
        prev = curr;
        curr = next;
    }
    if (min_diff_out) {
        *min_diff_out = min_diff;
    }
    return true;
}

static long double expected_cost(long double d1) {
    long double prev = 0.0L;
    long double curr = d1;
    long double sum = 0.0L;
    for (int i = 0; i < MAX_STEPS; ++i) {
        sum += curr * expl(-prev);
        if (curr > STOP_DEPTH) {
            break;
        }
        long double next = expl(curr - prev);
        prev = curr;
        curr = next;
    }
    return sum;
}

static long double find_high_threshold() {
    const long double scan_start = 0.6L;
    const long double scan_end = 1.5L;
    const long double scan_step = 0.01L;

    bool prev_inc = is_increasing(scan_start);
    long double low = scan_start;
    long double high = scan_start;
    for (long double d = scan_start + scan_step; d <= scan_end; d += scan_step) {
        bool inc = is_increasing(d);
        if (!prev_inc && inc) {
            low = d - scan_step;
            high = d;
            break;
        }
        prev_inc = inc;
    }
    if (high == scan_start) {
        return -1.0L;
    }

    for (int i = 0; i < 90; ++i) {
        long double mid = (low + high) * 0.5L;
        if (is_increasing(mid)) {
            high = mid;
        } else {
            low = mid;
        }
    }
    return high;
}

static bool run_validation(long double d1, long double answer) {
    bool ok = true;
    if (!(d1 > 0.7L && d1 < 0.8L)) {
        cerr << "Validation failed: threshold out of expected bracket.\n";
        ok = false;
    }
    if (!is_increasing(d1)) {
        cerr << "Validation failed: threshold not increasing.\n";
        ok = false;
    }
    long double eps = 1e-6L;
    if (is_increasing(d1 - eps)) {
        cerr << "Validation failed: threshold not minimal.\n";
        ok = false;
    }
    long double answer_up = expected_cost(d1 + eps);
    if (!(answer_up >= answer - 1e-12L)) {
        cerr << "Validation failed: expected cost not monotone near threshold.\n";
        ok = false;
    }

    long double best_low = numeric_limits<long double>::infinity();
    for (long double x = 0.02L; x <= 0.3L; x += 0.001L) {
        if (!is_increasing(x)) {
            continue;
        }
        long double v = expected_cost(x);
        if (v < best_low) {
            best_low = v;
        }
    }
    if (!(best_low > answer + 1e-4L)) {
        cerr << "Validation failed: low-branch minimum not higher than optimum.\n";
        ok = false;
    }
    if (!(answer > 2.35L && answer < 2.38L)) {
        cerr << "Validation failed: answer outside sanity range.\n";
        ok = false;
    }
    if (ok) {
        cerr << "Validation checkpoints passed.\n";
    }
    return ok;
}

int main() {
    long double d1 = find_high_threshold();
    if (d1 <= 0.0L) {
        cerr << "Failed to locate threshold.\n";
        return 1;
    }
    long double answer = expected_cost(d1);

    if (!run_validation(d1, answer)) {
        return 1;
    }

    cout.setf(ios::fixed);
    cout << setprecision(9) << answer << "\n";
    return 0;
}
