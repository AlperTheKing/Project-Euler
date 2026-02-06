#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

namespace {

constexpr int LFG_MOD = 1000000;
constexpr int POS_MOD = 10000;
constexpr int LEN_MOD = 399;

struct Cuboid {
    int x0;
    int x1;
    int y0;
    int y1;
    int z0;
    int z1;
};

struct RectYZ {
    std::uint16_t y0;
    std::uint16_t y1;
    std::uint16_t z0;
    std::uint16_t z1;
};

struct Event {
    int y;
    int z0;
    int z1;
    int delta;

    bool operator<(const Event& other) const {
        return y < other.y;
    }
};

class AreaUnion {
public:
    explicit AreaUnion(int max_z)
        : max_z_(max_z), cover_(static_cast<std::size_t>(4 * max_z + 8), 0),
          length_(static_cast<std::size_t>(4 * max_z + 8), 0) {}

    std::uint64_t compute(const std::vector<RectYZ>& rects) {
        if (rects.empty()) return 0;

        events_.clear();
        events_.reserve(rects.size() * 2);
        for (const RectYZ& r : rects) {
            events_.push_back({static_cast<int>(r.y0), static_cast<int>(r.z0), static_cast<int>(r.z1), +1});
            events_.push_back({static_cast<int>(r.y1), static_cast<int>(r.z0), static_cast<int>(r.z1), -1});
        }
        std::sort(events_.begin(), events_.end());

        std::fill(cover_.begin(), cover_.end(), 0);
        std::fill(length_.begin(), length_.end(), 0);

        std::uint64_t area = 0;
        int prev_y = events_.front().y;
        std::size_t i = 0;

        while (i < events_.size()) {
            const int y = events_[i].y;
            area += static_cast<std::uint64_t>(length_[1]) * static_cast<std::uint64_t>(y - prev_y);

            while (i < events_.size() && events_[i].y == y) {
                update(1, 0, max_z_, events_[i].z0, events_[i].z1, events_[i].delta);
                ++i;
            }
            prev_y = y;
        }

        return area;
    }

private:
    void update(int node, int l, int r, int ql, int qr, int delta) {
        if (qr <= l || r <= ql) return;

        if (ql <= l && r <= qr) {
            cover_[static_cast<std::size_t>(node)] += delta;
        } else {
            const int mid = l + (r - l) / 2;
            update(node * 2, l, mid, ql, qr, delta);
            update(node * 2 + 1, mid, r, ql, qr, delta);
        }

        if (cover_[static_cast<std::size_t>(node)] > 0) {
            length_[static_cast<std::size_t>(node)] = r - l;
        } else if (r - l == 1) {
            length_[static_cast<std::size_t>(node)] = 0;
        } else {
            length_[static_cast<std::size_t>(node)] =
                length_[static_cast<std::size_t>(node * 2)] + length_[static_cast<std::size_t>(node * 2 + 1)];
        }
    }

    int max_z_;
    std::vector<int> cover_;
    std::vector<int> length_;
    std::vector<Event> events_;
};

std::vector<int> build_lfg_sequence(int needed) {
    std::vector<int> s(static_cast<std::size_t>(needed + 1), 0);

    for (int k = 1; k <= std::min(55, needed); ++k) {
        long long value = 100003LL - 200003LL * k + 300007LL * k * k * k;
        value %= LFG_MOD;
        if (value < 0) value += LFG_MOD;
        s[static_cast<std::size_t>(k)] = static_cast<int>(value);
    }

    for (int k = 56; k <= needed; ++k) {
        s[static_cast<std::size_t>(k)] =
            (s[static_cast<std::size_t>(k - 24)] + s[static_cast<std::size_t>(k - 55)]) % LFG_MOD;
    }

    return s;
}

std::vector<Cuboid> generate_cuboids(int n) {
    const int needed = 6 * n;
    std::vector<int> s = build_lfg_sequence(needed);

    std::vector<Cuboid> cuboids;
    cuboids.reserve(static_cast<std::size_t>(n));

    for (int i = 1; i <= n; ++i) {
        const int base = 6 * i;
        const int x0 = s[static_cast<std::size_t>(base - 5)] % POS_MOD;
        const int y0 = s[static_cast<std::size_t>(base - 4)] % POS_MOD;
        const int z0 = s[static_cast<std::size_t>(base - 3)] % POS_MOD;

        const int dx = 1 + (s[static_cast<std::size_t>(base - 2)] % LEN_MOD);
        const int dy = 1 + (s[static_cast<std::size_t>(base - 1)] % LEN_MOD);
        const int dz = 1 + (s[static_cast<std::size_t>(base)] % LEN_MOD);

        cuboids.push_back({x0, x0 + dx, y0, y0 + dy, z0, z0 + dz});
    }

    return cuboids;
}

std::uint64_t combined_volume(const std::vector<Cuboid>& cuboids, unsigned thread_count) {
    if (cuboids.empty()) return 0;

    int max_x = 0;
    int max_z = 0;
    for (const Cuboid& c : cuboids) {
        max_x = std::max(max_x, c.x1);
        max_z = std::max(max_z, c.z1);
    }

    std::vector<int> diff(static_cast<std::size_t>(max_x + 1), 0);
    for (const Cuboid& c : cuboids) {
        ++diff[static_cast<std::size_t>(c.x0)];
        --diff[static_cast<std::size_t>(c.x1)];
    }

    std::vector<int> counts(static_cast<std::size_t>(max_x), 0);
    int active = 0;
    for (int x = 0; x < max_x; ++x) {
        active += diff[static_cast<std::size_t>(x)];
        counts[static_cast<std::size_t>(x)] = active;
    }

    std::vector<std::vector<RectYZ>> slabs(static_cast<std::size_t>(max_x));
    for (int x = 0; x < max_x; ++x) {
        slabs[static_cast<std::size_t>(x)].reserve(static_cast<std::size_t>(counts[static_cast<std::size_t>(x)]));
    }

    for (const Cuboid& c : cuboids) {
        const RectYZ r{static_cast<std::uint16_t>(c.y0), static_cast<std::uint16_t>(c.y1),
                       static_cast<std::uint16_t>(c.z0), static_cast<std::uint16_t>(c.z1)};
        for (int x = c.x0; x < c.x1; ++x) {
            slabs[static_cast<std::size_t>(x)].push_back(r);
        }
    }

    if (thread_count == 0) thread_count = 1;
    thread_count = std::min<unsigned>(thread_count, static_cast<unsigned>(max_x));

    std::atomic<int> next_x{0};
    std::vector<std::uint64_t> partial(static_cast<std::size_t>(thread_count), 0);
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(thread_count));

    for (unsigned tid = 0; tid < thread_count; ++tid) {
        workers.emplace_back([&, tid]() {
            AreaUnion area_union(max_z);
            std::uint64_t local = 0;

            while (true) {
                const int x = next_x.fetch_add(1, std::memory_order_relaxed);
                if (x >= max_x) break;
                local += area_union.compute(slabs[static_cast<std::size_t>(x)]);
            }

            partial[static_cast<std::size_t>(tid)] = local;
        });
    }

    for (auto& th : workers) th.join();

    std::uint64_t total = 0;
    for (std::uint64_t v : partial) total += v;
    return total;
}

std::uint64_t brute_force_volume(const std::vector<Cuboid>& cuboids) {
    int max_x = 0;
    int max_y = 0;
    int max_z = 0;
    for (const Cuboid& c : cuboids) {
        max_x = std::max(max_x, c.x1);
        max_y = std::max(max_y, c.y1);
        max_z = std::max(max_z, c.z1);
    }

    const int xy = max_x * max_y;
    std::vector<unsigned char> occ(static_cast<std::size_t>(xy * max_z), 0);

    for (const Cuboid& c : cuboids) {
        for (int x = c.x0; x < c.x1; ++x) {
            for (int y = c.y0; y < c.y1; ++y) {
                const int base = (x * max_y + y) * max_z;
                for (int z = c.z0; z < c.z1; ++z) {
                    occ[static_cast<std::size_t>(base + z)] = 1;
                }
            }
        }
    }

    std::uint64_t total = 0;
    for (unsigned char b : occ) total += b;
    return total;
}

void validation_checks() {
    {
        const std::vector<Cuboid> c = generate_cuboids(2);
        assert(c[0].x0 == 7 && c[0].y0 == 53 && c[0].z0 == 183);
        assert(c[0].x1 - c[0].x0 == 94 && c[0].y1 - c[0].y0 == 369 && c[0].z1 - c[0].z0 == 56);
        assert(c[1].x0 == 2383 && c[1].y0 == 3563 && c[1].z0 == 5079);
        assert(c[1].x1 - c[1].x0 == 42 && c[1].y1 - c[1].y0 == 212 && c[1].z1 - c[1].z0 == 344);
    }

    {
        const std::vector<Cuboid> c100 = generate_cuboids(100);
        assert(combined_volume(c100, 1) == 723581599ULL);
    }

    {
        std::uint32_t seed = 1;
        auto next_rand = [&seed]() {
            seed = seed * 1664525u + 1013904223u;
            return seed;
        };

        for (int tc = 0; tc < 8; ++tc) {
            std::vector<Cuboid> small;
            for (int i = 0; i < 30; ++i) {
                const int x0 = static_cast<int>(next_rand() % 12);
                const int y0 = static_cast<int>(next_rand() % 12);
                const int z0 = static_cast<int>(next_rand() % 12);
                const int dx = 1 + static_cast<int>(next_rand() % 4);
                const int dy = 1 + static_cast<int>(next_rand() % 4);
                const int dz = 1 + static_cast<int>(next_rand() % 4);
                small.push_back({x0, x0 + dx, y0, y0 + dy, z0, z0 + dz});
            }

            const std::uint64_t fast = combined_volume(small, 1);
            const std::uint64_t brute = brute_force_volume(small);
            assert(fast == brute);
        }
    }
}

}  // namespace

int main() {
    validation_checks();

    const std::vector<Cuboid> cuboids = generate_cuboids(50000);

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;

    const std::uint64_t answer = combined_volume(cuboids, threads);
    std::cout << answer << '\n';
    return 0;
}
