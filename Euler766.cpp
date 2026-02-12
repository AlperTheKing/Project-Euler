#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr int DIRS = 4;
constexpr std::array<int, DIRS> DR = {-1, 1, 0, 0};
constexpr std::array<int, DIRS> DC = {0, 0, -1, 1};

std::array<std::array<u64, 8>, 31> build_binom() {
    std::array<std::array<u64, 8>, 31> c{};
    for (int n = 0; n <= 30; ++n) {
        c[n][0] = 1;
        for (int k = 1; k <= 7; ++k) {
            c[n][k] = (k > n) ? 0 : c[n - 1][k - 1] + c[n - 1][k];
        }
    }
    return c;
}

const std::array<std::array<u64, 8>, 31> kBinom = build_binom();

template <std::size_t K>
u64 rank_combination(const std::array<std::uint8_t, K>& a) {
    u64 rank = 0;
    for (std::size_t i = 0; i < K; ++i) {
        rank += kBinom[a[i]][i + 1];
    }
    return rank;
}

struct Shape {
    int h = 0;
    int w = 0;
    int max_r = 0;
    int max_c = 0;
    int rows = 0;
    int cols = 0;
    std::vector<u32> masks;
    std::vector<std::array<int16_t, DIRS>> next;

    int idx_of(int r, int c) const {
        return r * cols + c;
    }
};

Shape build_shape(const int h, const int w, const std::vector<std::pair<int, int>>& cells) {
    Shape s;
    s.h = h;
    s.w = w;
    for (const auto& [r, c] : cells) {
        s.max_r = std::max(s.max_r, r);
        s.max_c = std::max(s.max_c, c);
    }
    s.rows = h - s.max_r;
    s.cols = w - s.max_c;
    const int total = s.rows * s.cols;
    s.masks.resize(total);
    s.next.resize(total);

    for (int ar = 0; ar < s.rows; ++ar) {
        for (int ac = 0; ac < s.cols; ++ac) {
            const int idx = s.idx_of(ar, ac);
            u32 mask = 0;
            for (const auto& [dr, dc] : cells) {
                const int rr = ar + dr;
                const int cc = ac + dc;
                mask |= (1u << (rr * w + cc));
            }
            s.masks[idx] = mask;

            for (int d = 0; d < DIRS; ++d) {
                const int nr = ar + DR[d];
                const int nc = ac + DC[d];
                if (nr >= 0 && nr < s.rows && nc >= 0 && nc < s.cols) {
                    s.next[idx][d] = static_cast<int16_t>(s.idx_of(nr, nc));
                } else {
                    s.next[idx][d] = -1;
                }
            }
        }
    }

    return s;
}

struct TargetState {
    std::array<std::uint8_t, 2> red{};
    std::array<std::uint8_t, 2> green{};
    std::array<std::uint8_t, 2> yellow{};
    std::array<std::uint8_t, 6> pink{};
    std::uint8_t blue = 0;
    std::uint8_t cyan = 0;
};

void canonicalize(TargetState& s) {
    if (s.red[0] > s.red[1]) {
        std::swap(s.red[0], s.red[1]);
    }
    if (s.green[0] > s.green[1]) {
        std::swap(s.green[0], s.green[1]);
    }
    if (s.yellow[0] > s.yellow[1]) {
        std::swap(s.yellow[0], s.yellow[1]);
    }
    std::sort(s.pink.begin(), s.pink.end());
}

u64 encode_target(const TargetState& s) {
    constexpr u64 kR = 190;      // C(20,2)
    constexpr u64 kG = 190;      // C(20,2)
    constexpr u64 kY = 276;      // C(24,2)
    constexpr u64 kP = 593775;   // C(30,6)
    constexpr u64 kB = 20;
    constexpr u64 kC = 25;

    const u64 rr = rank_combination(s.red);
    const u64 gg = rank_combination(s.green);
    const u64 yy = rank_combination(s.yellow);
    const u64 pp = rank_combination(s.pink);

    u64 key = rr;
    key = key * kG + gg;
    key = key * kY + yy;
    key = key * kP + pp;
    key = key * kB + s.blue;
    key = key * kC + s.cyan;
    return key;
}

u64 solve_target() {
    constexpr int H = 5;
    constexpr int W = 6;

    const Shape red_l = build_shape(H, W, {{0, 0}, {0, 1}, {1, 0}});
    const Shape green_l = build_shape(H, W, {{0, 1}, {1, 0}, {1, 1}});
    const Shape v_domino = build_shape(H, W, {{0, 0}, {1, 0}});
    const Shape single = build_shape(H, W, {{0, 0}});
    const Shape blue_sq = build_shape(H, W, {{0, 0}, {0, 1}, {1, 0}, {1, 1}});
    const Shape h_domino = build_shape(H, W, {{0, 0}, {0, 1}});

    TargetState start;
    start.red = {static_cast<std::uint8_t>(red_l.idx_of(0, 1)),
                 static_cast<std::uint8_t>(red_l.idx_of(0, 4))};
    start.green = {static_cast<std::uint8_t>(green_l.idx_of(0, 2)),
                   static_cast<std::uint8_t>(green_l.idx_of(3, 4))};
    start.yellow = {static_cast<std::uint8_t>(v_domino.idx_of(1, 5)),
                    static_cast<std::uint8_t>(v_domino.idx_of(2, 4))};
    start.pink = {static_cast<std::uint8_t>(single.idx_of(2, 0)),
                  static_cast<std::uint8_t>(single.idx_of(2, 1)),
                  static_cast<std::uint8_t>(single.idx_of(3, 0)),
                  static_cast<std::uint8_t>(single.idx_of(3, 1)),
                  static_cast<std::uint8_t>(single.idx_of(4, 0)),
                  static_cast<std::uint8_t>(single.idx_of(4, 1))};
    start.blue = static_cast<std::uint8_t>(blue_sq.idx_of(2, 2));
    start.cyan = static_cast<std::uint8_t>(h_domino.idx_of(4, 2));
    canonicalize(start);

    std::unordered_set<u64> seen;
    seen.reserve(3'000'000);

    std::vector<TargetState> q;
    q.reserve(3'000'000);

    seen.insert(encode_target(start));
    q.push_back(start);

    auto push_state = [&](const TargetState& ns) {
        const u64 key = encode_target(ns);
        if (seen.insert(key).second) {
            q.push_back(ns);
        }
    };

    for (std::size_t head = 0; head < q.size(); ++head) {
        const TargetState cur = q[head];

        const u32 m_r0 = red_l.masks[cur.red[0]];
        const u32 m_r1 = red_l.masks[cur.red[1]];
        const u32 m_g0 = green_l.masks[cur.green[0]];
        const u32 m_g1 = green_l.masks[cur.green[1]];
        const u32 m_y0 = v_domino.masks[cur.yellow[0]];
        const u32 m_y1 = v_domino.masks[cur.yellow[1]];
        const std::array<u32, 6> m_p = {
            single.masks[cur.pink[0]],
            single.masks[cur.pink[1]],
            single.masks[cur.pink[2]],
            single.masks[cur.pink[3]],
            single.masks[cur.pink[4]],
            single.masks[cur.pink[5]],
        };
        const u32 m_b = blue_sq.masks[cur.blue];
        const u32 m_c = h_domino.masks[cur.cyan];

        u32 occ = 0;
        occ |= m_r0;
        occ |= m_r1;
        occ |= m_g0;
        occ |= m_g1;
        occ |= m_y0;
        occ |= m_y1;
        for (u32 m : m_p) {
            occ |= m;
        }
        occ |= m_b;
        occ |= m_c;

        for (int i = 0; i < 2; ++i) {
            const std::uint8_t base = cur.red[i];
            const u32 own = (i == 0) ? m_r0 : m_r1;
            const u32 occ_wo = occ ^ own;
            for (int d = 0; d < DIRS; ++d) {
                int cur_idx = base;
                while (true) {
                    const int next_idx = red_l.next[cur_idx][d];
                    if (next_idx < 0) {
                        break;
                    }
                    const u32 nm = red_l.masks[next_idx];
                    if ((nm & occ_wo) != 0) {
                        break;
                    }
                    TargetState ns = cur;
                    ns.red[i] = static_cast<std::uint8_t>(next_idx);
                    if (ns.red[0] > ns.red[1]) {
                        std::swap(ns.red[0], ns.red[1]);
                    }
                    push_state(ns);
                    cur_idx = next_idx;
                }
            }
        }

        for (int i = 0; i < 2; ++i) {
            const std::uint8_t base = cur.green[i];
            const u32 own = (i == 0) ? m_g0 : m_g1;
            const u32 occ_wo = occ ^ own;
            for (int d = 0; d < DIRS; ++d) {
                int cur_idx = base;
                while (true) {
                    const int next_idx = green_l.next[cur_idx][d];
                    if (next_idx < 0) {
                        break;
                    }
                    const u32 nm = green_l.masks[next_idx];
                    if ((nm & occ_wo) != 0) {
                        break;
                    }
                    TargetState ns = cur;
                    ns.green[i] = static_cast<std::uint8_t>(next_idx);
                    if (ns.green[0] > ns.green[1]) {
                        std::swap(ns.green[0], ns.green[1]);
                    }
                    push_state(ns);
                    cur_idx = next_idx;
                }
            }
        }

        for (int i = 0; i < 2; ++i) {
            const std::uint8_t base = cur.yellow[i];
            const u32 own = (i == 0) ? m_y0 : m_y1;
            const u32 occ_wo = occ ^ own;
            for (int d = 0; d < DIRS; ++d) {
                int cur_idx = base;
                while (true) {
                    const int next_idx = v_domino.next[cur_idx][d];
                    if (next_idx < 0) {
                        break;
                    }
                    const u32 nm = v_domino.masks[next_idx];
                    if ((nm & occ_wo) != 0) {
                        break;
                    }
                    TargetState ns = cur;
                    ns.yellow[i] = static_cast<std::uint8_t>(next_idx);
                    if (ns.yellow[0] > ns.yellow[1]) {
                        std::swap(ns.yellow[0], ns.yellow[1]);
                    }
                    push_state(ns);
                    cur_idx = next_idx;
                }
            }
        }

        for (int i = 0; i < 6; ++i) {
            const std::uint8_t base = cur.pink[i];
            const u32 own = m_p[i];
            const u32 occ_wo = occ ^ own;
            for (int d = 0; d < DIRS; ++d) {
                int cur_idx = base;
                while (true) {
                    const int next_idx = single.next[cur_idx][d];
                    if (next_idx < 0) {
                        break;
                    }
                    const u32 nm = single.masks[next_idx];
                    if ((nm & occ_wo) != 0) {
                        break;
                    }
                    TargetState ns = cur;
                    ns.pink[i] = static_cast<std::uint8_t>(next_idx);
                    std::sort(ns.pink.begin(), ns.pink.end());
                    push_state(ns);
                    cur_idx = next_idx;
                }
            }
        }

        {
            const std::uint8_t base = cur.blue;
            const u32 occ_wo = occ ^ m_b;
            for (int d = 0; d < DIRS; ++d) {
                int cur_idx = base;
                while (true) {
                    const int next_idx = blue_sq.next[cur_idx][d];
                    if (next_idx < 0) {
                        break;
                    }
                    const u32 nm = blue_sq.masks[next_idx];
                    if ((nm & occ_wo) != 0) {
                        break;
                    }
                    TargetState ns = cur;
                    ns.blue = static_cast<std::uint8_t>(next_idx);
                    push_state(ns);
                    cur_idx = next_idx;
                }
            }
        }

        {
            const std::uint8_t base = cur.cyan;
            const u32 occ_wo = occ ^ m_c;
            for (int d = 0; d < DIRS; ++d) {
                int cur_idx = base;
                while (true) {
                    const int next_idx = h_domino.next[cur_idx][d];
                    if (next_idx < 0) {
                        break;
                    }
                    const u32 nm = h_domino.masks[next_idx];
                    if ((nm & occ_wo) != 0) {
                        break;
                    }
                    TargetState ns = cur;
                    ns.cyan = static_cast<std::uint8_t>(next_idx);
                    push_state(ns);
                    cur_idx = next_idx;
                }
            }
        }
    }

    return seen.size();
}

struct SampleState {
    std::uint8_t l = 0;
    std::array<std::uint8_t, 7> s{};
};

void canonicalize(SampleState& st) {
    std::sort(st.s.begin(), st.s.end());
}

u64 encode_sample(const SampleState& st) {
    constexpr u64 kS = 792;  // C(12,7)
    constexpr u64 kL = 6;
    return rank_combination(st.s) * kL + st.l;
}

u64 solve_sample() {
    constexpr int H = 3;
    constexpr int W = 4;

    const Shape l_shape = build_shape(H, W, {{0, 0}, {0, 1}, {1, 0}});
    const Shape single = build_shape(H, W, {{0, 0}});

    SampleState start;
    start.l = static_cast<std::uint8_t>(l_shape.idx_of(0, 0));
    start.s = {
        static_cast<std::uint8_t>(single.idx_of(0, 2)),
        static_cast<std::uint8_t>(single.idx_of(1, 1)),
        static_cast<std::uint8_t>(single.idx_of(1, 2)),
        static_cast<std::uint8_t>(single.idx_of(2, 0)),
        static_cast<std::uint8_t>(single.idx_of(2, 1)),
        static_cast<std::uint8_t>(single.idx_of(2, 2)),
        static_cast<std::uint8_t>(single.idx_of(2, 3)),
    };
    canonicalize(start);

    std::unordered_set<u64> seen;
    seen.reserve(1024);
    std::vector<SampleState> q;
    q.reserve(1024);

    seen.insert(encode_sample(start));
    q.push_back(start);

    auto push_state = [&](const SampleState& ns) {
        const u64 key = encode_sample(ns);
        if (seen.insert(key).second) {
            q.push_back(ns);
        }
    };

    for (std::size_t head = 0; head < q.size(); ++head) {
        const SampleState cur = q[head];

        const u32 m_l = l_shape.masks[cur.l];
        const std::array<u32, 7> m_s = {
            single.masks[cur.s[0]],
            single.masks[cur.s[1]],
            single.masks[cur.s[2]],
            single.masks[cur.s[3]],
            single.masks[cur.s[4]],
            single.masks[cur.s[5]],
            single.masks[cur.s[6]],
        };

        u32 occ = m_l;
        for (u32 m : m_s) {
            occ |= m;
        }

        {
            const u32 occ_wo = occ ^ m_l;
            for (int d = 0; d < DIRS; ++d) {
                int cur_idx = cur.l;
                while (true) {
                    const int next_idx = l_shape.next[cur_idx][d];
                    if (next_idx < 0) {
                        break;
                    }
                    const u32 nm = l_shape.masks[next_idx];
                    if ((nm & occ_wo) != 0) {
                        break;
                    }
                    SampleState ns = cur;
                    ns.l = static_cast<std::uint8_t>(next_idx);
                    push_state(ns);
                    cur_idx = next_idx;
                }
            }
        }

        for (int i = 0; i < 7; ++i) {
            const u32 occ_wo = occ ^ m_s[i];
            for (int d = 0; d < DIRS; ++d) {
                int cur_idx = cur.s[i];
                while (true) {
                    const int next_idx = single.next[cur_idx][d];
                    if (next_idx < 0) {
                        break;
                    }
                    const u32 nm = single.masks[next_idx];
                    if ((nm & occ_wo) != 0) {
                        break;
                    }
                    SampleState ns = cur;
                    ns.s[i] = static_cast<std::uint8_t>(next_idx);
                    canonicalize(ns);
                    push_state(ns);
                    cur_idx = next_idx;
                }
            }
        }
    }

    return seen.size();
}

}  // namespace

int main() {
    assert(solve_sample() == 208ULL);
    std::cout << solve_target() << '\n';
    return 0;
}
