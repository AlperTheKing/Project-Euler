#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
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
    int dx = 0;
    int dy = 0;
    int length = 0;
    int max_multiple = 0;
};

struct DecodedState {
    int x = 0;
    int y = 0;
    int perimeter = 0;
    int edges_bucket = 0;
    u64 count = 0ULL;
};

using StateMap = std::unordered_map<u64, u64>;

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

unsigned pick_thread_count(const unsigned requested) {
    if (requested > 0U) {
        return requested;
    }

    // Hash-map dominated transitions are memory-bandwidth bound here; the
    // threaded path is available via --threads=N for machines where it helps.
    return 1U;
}

int integer_square_root_if_square(const int value) {
    if (value < 0) {
        return -1;
    }

    int root = static_cast<int>(std::sqrt(static_cast<long double>(value)));
    while ((static_cast<long long>(root) + 1LL) * (static_cast<long long>(root) + 1LL) <=
           static_cast<long long>(value)) {
        ++root;
    }
    while (static_cast<long long>(root) * static_cast<long long>(root) >
           static_cast<long long>(value)) {
        --root;
    }

    if (static_cast<long long>(root) * static_cast<long long>(root) ==
        static_cast<long long>(value)) {
        return root;
    }
    return -1;
}

std::vector<Direction> build_directions(const int n) {
    std::vector<Direction> directions;
    directions.reserve(static_cast<std::size_t>(4 * n));

    struct PolarDirection {
        Direction direction;
        long double angle = 0.0L;
    };

    std::vector<PolarDirection> raw;
    raw.reserve(static_cast<std::size_t>(4 * n));

    for (int dx = -n; dx <= n; ++dx) {
        for (int dy = -n; dy <= n; ++dy) {
            if (dx == 0 && dy == 0) {
                continue;
            }
            if (std::gcd(std::abs(dx), std::abs(dy)) != 1) {
                continue;
            }

            const int length_sq = dx * dx + dy * dy;
            const int length = integer_square_root_if_square(length_sq);
            if (length <= 0 || length > n) {
                continue;
            }

            PolarDirection pd;
            pd.direction.dx = dx;
            pd.direction.dy = dy;
            pd.direction.length = length;
            pd.direction.max_multiple = n / length;
            pd.angle = std::atan2(static_cast<long double>(dy),
                                  static_cast<long double>(dx));
            raw.push_back(pd);
        }
    }

    std::sort(raw.begin(),
              raw.end(),
              [](const PolarDirection& lhs, const PolarDirection& rhs) {
                  if (lhs.angle != rhs.angle) {
                      return lhs.angle < rhs.angle;
                  }
                  if (lhs.direction.dx != rhs.direction.dx) {
                      return lhs.direction.dx < rhs.direction.dx;
                  }
                  return lhs.direction.dy < rhs.direction.dy;
              });

    directions.reserve(raw.size());
    for (const PolarDirection& pd : raw) {
        directions.push_back(pd.direction);
    }

    return directions;
}

u64 encode_state(const int x,
                 const int y,
                 const int perimeter,
                 const int edges_bucket,
                 const int n) {
    const u64 width = static_cast<u64>(2 * n + 1);
    const u64 perimeter_base = static_cast<u64>(n + 1);

    const u64 x_off = static_cast<u64>(x + n);
    const u64 y_off = static_cast<u64>(y + n);

    return (((static_cast<u64>(edges_bucket) * perimeter_base +
              static_cast<u64>(perimeter)) *
                 width +
             x_off) *
                width +
            y_off);
}

DecodedState decode_state(const u64 key, const int n) {
    const u64 width = static_cast<u64>(2 * n + 1);
    const u64 perimeter_base = static_cast<u64>(n + 1);

    DecodedState decoded;

    u64 tmp = key;
    const u64 y_off = tmp % width;
    tmp /= width;
    const u64 x_off = tmp % width;
    tmp /= width;

    decoded.perimeter = static_cast<int>(tmp % perimeter_base);
    tmp /= perimeter_base;

    decoded.edges_bucket = static_cast<int>(tmp);
    decoded.x = static_cast<int>(x_off) - n;
    decoded.y = static_cast<int>(y_off) - n;

    return decoded;
}

void add_single_direction_transitions(const std::vector<DecodedState>& states,
                                      const std::size_t begin,
                                      const std::size_t end,
                                      const Direction& direction,
                                      const int n,
                                      StateMap& out) {
    for (std::size_t i = begin; i < end; ++i) {
        const DecodedState& state = states[i];

        const int next_edges_bucket = (state.edges_bucket >= 3)
                                          ? 3
                                          : (state.edges_bucket + 1);

        int x = state.x + direction.dx;
        int y = state.y + direction.dy;
        int perimeter = state.perimeter + direction.length;

        for (int mult = 1;
             mult <= direction.max_multiple && perimeter <= n;
             ++mult) {
            if (std::abs(x) <= n && std::abs(y) <= n) {
                const u64 key =
                    encode_state(x, y, perimeter, next_edges_bucket, n);
                out[key] += state.count;
            }

            x += direction.dx;
            y += direction.dy;
            perimeter += direction.length;
        }
    }
}

StateMap apply_direction(const StateMap& current,
                         const Direction& direction,
                         const int n,
                         const unsigned thread_count) {
    std::vector<DecodedState> decoded_states;
    decoded_states.reserve(current.size());
    for (const auto& [key, count] : current) {
        DecodedState decoded = decode_state(key, n);
        decoded.count = count;
        decoded_states.push_back(decoded);
    }

    StateMap next = current;

    const u64 estimated_work = static_cast<u64>(decoded_states.size()) *
                               static_cast<u64>(direction.max_multiple);
    const bool use_parallel =
        thread_count > 1U && decoded_states.size() >= 50000U &&
        estimated_work >= 600000U;

    if (!use_parallel) {
        add_single_direction_transitions(
            decoded_states, 0U, decoded_states.size(), direction, n, next);
        return next;
    }

    const unsigned workers =
        std::min<unsigned>(thread_count,
                           static_cast<unsigned>(decoded_states.size()));
    const std::size_t chunk_size =
        (decoded_states.size() + static_cast<std::size_t>(workers) - 1U) /
        static_cast<std::size_t>(workers);

    std::vector<StateMap> locals(static_cast<std::size_t>(workers));
    std::vector<std::thread> threads;
    threads.reserve(static_cast<std::size_t>(workers));

    for (unsigned worker = 0U; worker < workers; ++worker) {
        const std::size_t begin = static_cast<std::size_t>(worker) * chunk_size;
        const std::size_t end = std::min(decoded_states.size(), begin + chunk_size);
        if (begin >= end) {
            continue;
        }

        threads.emplace_back([
            &decoded_states,
            &direction,
            n,
            begin,
            end,
            &locals,
            worker
        ]() {
            StateMap& local = locals[static_cast<std::size_t>(worker)];
            local.reserve((end - begin) * static_cast<std::size_t>(2));
            add_single_direction_transitions(
                decoded_states, begin, end, direction, n, local);
        });
    }

    for (std::thread& t : threads) {
        t.join();
    }

    std::size_t additional_capacity = 0U;
    for (const StateMap& local : locals) {
        additional_capacity += local.size();
    }
    next.reserve(next.size() + additional_capacity);

    for (const StateMap& local : locals) {
        for (const auto& [key, count] : local) {
            next[key] += count;
        }
    }

    return next;
}

u64 count_polygons(const int n, const unsigned thread_count) {
    if (n < 3) {
        return 0ULL;
    }

    const std::vector<Direction> directions = build_directions(n);

    StateMap dp;
    dp.reserve(1024U);
    dp[encode_state(0, 0, 0, 0, n)] = 1ULL;

    for (const Direction& direction : directions) {
        dp = apply_direction(dp, direction, n, thread_count);
    }

    u64 answer = 0ULL;
    for (int p = 0; p <= n; ++p) {
        const u64 key = encode_state(0, 0, p, 3, n);
        const auto it = dp.find(key);
        if (it != dp.end()) {
            answer += it->second;
        }
    }

    return answer;
}

void run_checkpoints(const unsigned thread_count) {
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
        const u64 got = count_polygons(cp.n, thread_count);
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

        const unsigned thread_count = pick_thread_count(options.requested_threads);

        if (options.run_checkpoints) {
            run_checkpoints(thread_count);
        }

        const u64 answer = count_polygons(options.perimeter_limit, thread_count);
        std::cout << answer << '\n';
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
