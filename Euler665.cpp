#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>

namespace {

using u64 = std::uint64_t;

class SuccessorSet {
public:
    int next_free(int value) {
        if (empty_ || value < lo_ || value > hi_) {
            return value;
        }
        int& slot = parent_[static_cast<std::size_t>(value - lo_)];
        if (slot == kNone) {
            return value;
        }
        slot = next_free(slot);
        return slot;
    }

    void insert(int value) {
        ensure_contains(value);
        int& slot = parent_[static_cast<std::size_t>(value - lo_)];
        if (slot != kNone) {
            return;
        }
        slot = next_free(value + 1);
    }

    bool contains(int value) const {
        if (empty_ || value < lo_ || value > hi_) {
            return false;
        }
        return parent_[static_cast<std::size_t>(value - lo_)] != kNone;
    }

private:
    void ensure_contains(int value) {
        if (empty_) {
            lo_ = hi_ = value;
            parent_.assign(1, kNone);
            empty_ = false;
            return;
        }
        if (value >= lo_ && value <= hi_) {
            return;
        }

        int new_lo = lo_;
        int new_hi = hi_;
        while (value < new_lo || value > new_hi) {
            const int span = new_hi - new_lo + 1;
            if (value < new_lo) {
                new_lo -= span;
            }
            if (value > new_hi) {
                new_hi += span;
            }
        }

        std::vector<int> expanded(static_cast<std::size_t>(new_hi - new_lo + 1), kNone);
        const std::size_t shift = static_cast<std::size_t>(lo_ - new_lo);
        std::copy(parent_.begin(), parent_.end(),
                  expanded.begin() + static_cast<std::ptrdiff_t>(shift));
        parent_.swap(expanded);
        lo_ = new_lo;
        hi_ = new_hi;
    }

    static constexpr int kNone = std::numeric_limits<int>::min();
    bool empty_ = true;
    int lo_ = 0;
    int hi_ = -1;
    std::vector<int> parent_;
};

struct Event {
    int activate_x;
    int t_value;
    int s_value;
};

struct EventMinHeapOrder {
    bool operator()(const Event& lhs, const Event& rhs) const {
        return lhs.activate_x > rhs.activate_x;
    }
};

class Solver {
public:
    u64 f(int M) {
        SuccessorSet used;
        SuccessorSet dif_set;
        SuccessorSet low_set;
        SuccessorSet high_even_set;
        SuccessorSet high_odd_set;
        std::priority_queue<Event, std::vector<Event>, EventMinHeapOrder> pending;

        auto insert_s = [&](int s) {
            if ((s & 1) == 0) {
                high_even_set.insert(s / 2);
            } else {
                high_odd_set.insert((s - 1) / 2);
            }
        };

        dif_set.insert(0);
        low_set.insert(0);
        insert_s(0);

        auto activate_events = [&](int x) {
            while (!pending.empty() && pending.top().activate_x <= x) {
                const Event e = pending.top();
                pending.pop();
                low_set.insert(e.t_value);
                insert_s(e.s_value);
            }
        };

        u64 total = 0;
        int x = 1;

        while (x <= M) {
            activate_events(x);
            x = used.next_free(x);
            if (x > M) {
                break;
            }

            int y = x + 1;
            while (true) {
                const int y1 = used.next_free(y);
                const int y2 = x + dif_set.next_free(y1 - x);
                const int y3 = 2 * x + low_set.next_free(y2 - 2 * x);

                const int shift = ((x & 1) == 0) ? (x / 2) : (x / 2 + 1);
                const int y4 = ((x & 1) == 0)
                                   ? (shift + high_even_set.next_free(y3 - shift))
                                   : (shift + high_odd_set.next_free(y3 - shift));
                if (y4 == y) {
                    break;
                }
                y = y4;
            }

            used.insert(x);
            used.insert(y);

            if (static_cast<u64>(x) + static_cast<u64>(y) <= static_cast<u64>(M)) {
                total += static_cast<u64>(x + y);
            }

            dif_set.insert(y - x);
            low_set.insert(y - 2 * x);
            insert_s(2 * y - x);
            pending.push(Event{y + 1, x - 2 * y, 2 * x - y});

            ++x;
        }

        return total;
    }
};

u64 brute_f(int M) {
    const int N = M;
    std::vector<unsigned char> win(static_cast<std::size_t>(N + 1) * static_cast<std::size_t>(N + 1),
                                   0);
    auto at = [&](int a, int b) -> unsigned char& {
        return win[static_cast<std::size_t>(N + 1) * static_cast<std::size_t>(a) +
                   static_cast<std::size_t>(b)];
    };

    u64 total = 0;

    for (int sum = 0; sum <= M; ++sum) {
        for (int a = 0; a <= sum / 2; ++a) {
            const int b = sum - a;
            bool is_win = false;

            for (int n = 1; n <= a && !is_win; ++n) {
                int x = a - n;
                int y = b;
                if (x > y) std::swap(x, y);
                if (at(x, y) == 0) is_win = true;
            }
            for (int n = 1; n <= b && !is_win; ++n) {
                int x = a;
                int y = b - n;
                if (x > y) std::swap(x, y);
                if (at(x, y) == 0) is_win = true;
            }
            for (int n = 1; n <= a && !is_win; ++n) {
                int x = a - n;
                int y = b - n;
                if (x > y) std::swap(x, y);
                if (at(x, y) == 0) is_win = true;
            }
            for (int n = 1; n <= std::min(a, b / 2) && !is_win; ++n) {
                int x = a - n;
                int y = b - 2 * n;
                if (x > y) std::swap(x, y);
                if (at(x, y) == 0) is_win = true;
            }
            for (int n = 1; n <= std::min(a / 2, b) && !is_win; ++n) {
                int x = a - 2 * n;
                int y = b - n;
                if (x > y) std::swap(x, y);
                if (at(x, y) == 0) is_win = true;
            }

            at(a, b) = static_cast<unsigned char>(is_win ? 1 : 0);
            if (!is_win && sum > 0) {
                total += static_cast<u64>(sum);
            }
        }
    }

    return total;
}

}  // namespace

int main() {
    Solver solver;

    assert(brute_f(60) == solver.f(60));
    assert(solver.f(10) == 21ULL);
    assert(solver.f(100) == 1164ULL);
    assert(solver.f(1000) == 117002ULL);

    std::cout << solver.f(10'000'000) << "\n";
    return 0;
}
