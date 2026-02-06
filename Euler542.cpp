#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

namespace {

struct Candidate {
    u64 n = 0;          // Jump period: n = u^p.
    u64 d = 0;          // Jump amount: d = u^(p+1) - (u-1)^(p+1).
    u64 relevance = 0;  // Proven upper bound where candidate can still be optimal.
};

u64 pow_with_cap(const u64 base, const int exp, const u64 cap) {
    u128 value = 1;
    for (int i = 0; i < exp; ++i) {
        value *= static_cast<u128>(base);
        if (value > static_cast<u128>(cap)) {
            return cap + 1;
        }
    }
    return static_cast<u64>(value);
}

u64 int_root(const u64 n, const int p) {
    u64 lo = 1;
    u64 hi = 2;
    while (pow_with_cap(hi, p, n) <= n) {
        if (hi > n / 2) {
            hi = n;
            break;
        }
        hi *= 2;
    }
    while (lo + 1 < hi) {
        const u64 mid = lo + (hi - lo) / 2;
        if (pow_with_cap(mid, p, n) <= n) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

u64 delta_value(const u64 u, const int p) {
    // Computes u^(p+1) - (u-1)^(p+1) with a stable recurrence:
    // D_m = u * D_{m-1} + (u-1)^(m-1), D_1 = 1.
    const u64 v = u - 1;
    u128 D = 1;
    u128 vpow = 1;
    for (int m = 2; m <= p + 1; ++m) {
        vpow *= static_cast<u128>(v);
        D = static_cast<u128>(u) * D + vpow;
    }
    return static_cast<u64>(D);
}

int max_power_for_limit(const u64 n) {
    int p = 0;
    u64 x = 1;
    while (x <= n / 2) {
        x *= 2;
        ++p;
    }
    return p;
}

std::vector<Candidate> generate_record_candidates(const u64 limit) {
    if (limit < 4) {
        return {};
    }

    const int p_max = max_power_for_limit(limit);
    std::vector<u64> u_max(static_cast<std::size_t>(p_max + 1), 0);
    for (int p = 2; p <= p_max; ++p) {
        u_max[static_cast<std::size_t>(p)] = int_root(limit, p);
    }

    std::vector<Candidate> records;
    u64 current_best_d = 0;

    while (true) {
        bool found = false;
        u64 best_n = 0;
        u64 best_d = 0;

        for (int p = 2; p <= p_max; ++p) {
            const u64 umax = u_max[static_cast<std::size_t>(p)];
            if (umax < 2) {
                continue;
            }
            if (delta_value(umax, p) <= current_best_d) {
                continue;
            }

            u64 lo = 2;
            u64 hi = umax;
            while (lo < hi) {
                const u64 mid = lo + (hi - lo) / 2;
                if (delta_value(mid, p) > current_best_d) {
                    hi = mid;
                } else {
                    lo = mid + 1;
                }
            }

            const u64 u = lo;
            const u64 n = pow_with_cap(u, p, limit);
            const u64 d = delta_value(u, p);
            if (!found || n < best_n || (n == best_n && d > best_d)) {
                found = true;
                best_n = n;
                best_d = d;
            }
        }

        if (!found) {
            break;
        }

        records.push_back({best_n, best_d, limit});
        current_best_d = best_d;
    }

    return records;
}

void compute_relevance_bounds(std::vector<Candidate>& candidates, const u64 limit) {
    const std::size_t n = candidates.size();
    for (std::size_t i = 0; i < n; ++i) {
        u64 best = limit;
        const u64 ni = candidates[i].n;
        const u64 di = candidates[i].d;

        for (std::size_t j = 0; j < n; ++j) {
            const u64 nj = candidates[j].n;
            const u64 dj = candidates[j].d;

            const u128 lhs = static_cast<u128>(dj) * static_cast<u128>(ni);
            const u128 rhs = static_cast<u128>(di) * static_cast<u128>(nj);
            if (lhs <= rhs) {
                continue;
            }

            // For k > num/den, candidate j is guaranteed > candidate i:
            // num = dj*nj*ni, den = dj*ni - di*nj.
            const cpp_int num = cpp_int(dj) * cpp_int(nj) * cpp_int(ni);
            const cpp_int den = cpp_int(dj) * cpp_int(ni) - cpp_int(di) * cpp_int(nj);
            const cpp_int r = num / den;  // floor(num / den)

            if (r < cpp_int(best)) {
                best = static_cast<u64>(r);
            }
        }

        candidates[i].relevance = best;
    }
}

std::vector<u64> collect_events(const std::vector<Candidate>& candidates, const u64 limit) {
    std::vector<u64> events;
    events.reserve(16384);

    for (const Candidate& c : candidates) {
        const u64 lim = std::min(limit, c.relevance);
        if (lim < c.n) {
            continue;
        }

        u64 t = c.n;
        while (t <= lim) {
            events.push_back(t);
            if (t > lim - c.n) {
                break;
            }
            t += c.n;
        }
    }

    std::sort(events.begin(), events.end());
    events.erase(std::unique(events.begin(), events.end()), events.end());
    return events;
}

i128 to_i128(const u64 x) {
    return static_cast<i128>(x);
}

std::string i128_to_string(i128 value) {
    if (value == 0) {
        return "0";
    }
    bool neg = value < 0;
    if (neg) {
        value = -value;
    }
    std::string out;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    if (neg) {
        out.push_back('-');
    }
    std::reverse(out.begin(), out.end());
    return out;
}

int alternating_sign_sum(const u64 l, const u64 r) {
    if (l > r) {
        return 0;
    }
    const u64 len = r - l + 1;
    if ((len & 1ULL) == 0ULL) {
        return 0;
    }
    return (l & 1ULL) ? -1 : 1;
}

class Solver542 {
public:
    explicit Solver542(const u64 limit) : limit_(limit) {
        candidates_ = generate_record_candidates(limit_);
        compute_relevance_bounds(candidates_, limit_);
        events_ = collect_events(candidates_, limit_);
    }

    u64 S(const u64 k) const {
        if (k < 4) {
            return 0;
        }
        u128 best = 0;
        for (const Candidate& c : candidates_) {
            if (c.n > k) {
                break;
            }
            if (k > c.relevance) {
                continue;
            }
            const u128 value = static_cast<u128>(k / c.n) * static_cast<u128>(c.d);
            if (value > best) {
                best = value;
            }
        }
        return static_cast<u64>(best);
    }

    i128 T(const u64 n) const {
        if (n < 4) {
            return 0;
        }
        if (n > limit_) {
            throw std::runtime_error("Requested n exceeds solver limit.");
        }

        i128 answer = 0;
        u64 current = 4;
        u64 s_value = S(4);
        std::size_t event_idx = 0;

        while (current <= n) {
            const u64 next_event =
                (event_idx < events_.size() && events_[event_idx] <= n) ? events_[event_idx] : (n + 1);

            if (next_event > current) {
                const u64 r = std::min(n, next_event - 1);
                const int sign = alternating_sign_sum(current, r);
                if (sign != 0) {
                    answer += static_cast<i128>(sign) * to_i128(s_value);
                }
                current = r + 1;
                if (current > n) {
                    break;
                }
            }

            // At event points, at least one candidate jumps; recompute exact S(current).
            s_value = S(current);
            ++event_idx;
        }

        return answer;
    }

private:
    u64 limit_;
    std::vector<Candidate> candidates_;
    std::vector<u64> events_;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0) != 0) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    u64 parsed = 0;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(ch - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }
    value = parsed;
    return true;
}

void run_checkpoints() {
    Solver542 check_solver(1000);
    if (check_solver.S(4) != 7) {
        throw std::runtime_error("Checkpoint failed: S(4) != 7");
    }
    if (check_solver.S(10) != 19) {
        throw std::runtime_error("Checkpoint failed: S(10) != 19");
    }
    if (check_solver.S(12) != 21) {
        throw std::runtime_error("Checkpoint failed: S(12) != 21");
    }
    if (check_solver.S(1000) != 3439) {
        throw std::runtime_error("Checkpoint failed: S(1000) != 3439");
    }
    if (check_solver.T(1000) != 2268) {
        throw std::runtime_error("Checkpoint failed: T(1000) != 2268");
    }
}

}  // namespace

int main(int argc, char** argv) {
    u64 n = 100000000000000000ULL;
    bool run_checks = true;

    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            run_checks = false;
            continue;
        }
        u64 parsed = 0;
        if (parse_u64_after_prefix(arg, "--n=", parsed)) {
            n = parsed;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return 1;
    }

    try {
        if (run_checks) {
            run_checkpoints();
        }
        Solver542 solver(n);
        const i128 answer = solver.T(n);
        std::cout << i128_to_string(answer) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
