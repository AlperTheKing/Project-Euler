#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;

constexpr u64 kMod = 1'001'001'011ULL;

u64 add_mod(const u64 a, const u64 b) {
    const u64 s = a + b;
    if (s >= kMod) {
        return s - kMod;
    }
    return s;
}

u64 mul_mod(const u64 a, const u64 b) {
    return static_cast<u64>((static_cast<unsigned __int128>(a) * static_cast<unsigned __int128>(b)) %
                            static_cast<unsigned __int128>(kMod));
}

u64 solve(const int n) {
    const int N = 1 << n;

    std::vector<int> succ(static_cast<std::size_t>(N), 0);
    for (int x = 0; x < N; ++x) {
        const int b1 = x & 1;
        const int b2 = (x >> 1) & 1;
        const int b3 = (x >> 2) & 1;
        const int last = b1 & (b2 ^ b3);
        succ[static_cast<std::size_t>(x)] = (x >> 1) | (last << (n - 1));
    }

    std::vector<int> head(static_cast<std::size_t>(N), -1);
    std::vector<int> to(static_cast<std::size_t>(N), 0);
    std::vector<int> next_edge(static_cast<std::size_t>(N), -1);
    int edge_ptr = 0;
    for (int v = 0; v < N; ++v) {
        const int u = succ[static_cast<std::size_t>(v)];
        to[static_cast<std::size_t>(edge_ptr)] = v;
        next_edge[static_cast<std::size_t>(edge_ptr)] = head[static_cast<std::size_t>(u)];
        head[static_cast<std::size_t>(u)] = edge_ptr;
        ++edge_ptr;
    }

    std::vector<std::uint8_t> color(static_cast<std::size_t>(N), 0U);
    std::vector<int> pos(static_cast<std::size_t>(N), -1);
    std::vector<std::uint8_t> in_cycle(static_cast<std::size_t>(N), 0U);
    std::vector<std::vector<int>> cycles;
    cycles.reserve(N / 4);

    std::vector<int> path;
    path.reserve(256);

    for (int start = 0; start < N; ++start) {
        if (color[static_cast<std::size_t>(start)] != 0U) {
            continue;
        }

        path.clear();
        int v = start;
        while (color[static_cast<std::size_t>(v)] == 0U) {
            color[static_cast<std::size_t>(v)] = 1U;
            pos[static_cast<std::size_t>(v)] = static_cast<int>(path.size());
            path.push_back(v);
            v = succ[static_cast<std::size_t>(v)];
        }

        if (color[static_cast<std::size_t>(v)] == 1U) {
            const int cycle_start = pos[static_cast<std::size_t>(v)];
            std::vector<int> cyc;
            cyc.reserve(path.size() - static_cast<std::size_t>(cycle_start));
            for (int i = cycle_start; i < static_cast<int>(path.size()); ++i) {
                const int node = path[static_cast<std::size_t>(i)];
                in_cycle[static_cast<std::size_t>(node)] = 1U;
                cyc.push_back(node);
            }
            cycles.push_back(std::move(cyc));
        }

        for (const int node : path) {
            color[static_cast<std::size_t>(node)] = 2U;
            pos[static_cast<std::size_t>(node)] = -1;
        }
    }

    std::vector<u64> dp0(static_cast<std::size_t>(N), 0ULL);
    std::vector<u64> dp1(static_cast<std::size_t>(N), 0ULL);
    std::vector<std::uint8_t> done(static_cast<std::size_t>(N), 0U);

    auto compute_tree = [&](const int root) {
        if (done[static_cast<std::size_t>(root)] != 0U) {
            return;
        }

        std::vector<std::pair<int, bool>> st;
        st.reserve(128);
        st.push_back({root, false});

        while (!st.empty()) {
            const auto [u, expanded] = st.back();
            st.pop_back();

            if (!expanded) {
                st.push_back({u, true});
                for (int e = head[static_cast<std::size_t>(u)]; e != -1;
                     e = next_edge[static_cast<std::size_t>(e)]) {
                    const int ch = to[static_cast<std::size_t>(e)];
                    if (in_cycle[static_cast<std::size_t>(ch)] != 0U) {
                        continue;
                    }
                    if (done[static_cast<std::size_t>(ch)] == 0U) {
                        st.push_back({ch, false});
                    }
                }
                continue;
            }

            u64 ways0 = 1ULL;
            u64 ways1 = 1ULL;
            for (int e = head[static_cast<std::size_t>(u)]; e != -1;
                 e = next_edge[static_cast<std::size_t>(e)]) {
                const int ch = to[static_cast<std::size_t>(e)];
                if (in_cycle[static_cast<std::size_t>(ch)] != 0U) {
                    continue;
                }
                const u64 sum = add_mod(dp0[static_cast<std::size_t>(ch)], dp1[static_cast<std::size_t>(ch)]);
                ways0 = mul_mod(ways0, sum);
                ways1 = mul_mod(ways1, dp0[static_cast<std::size_t>(ch)]);
            }

            dp0[static_cast<std::size_t>(u)] = ways0;
            dp1[static_cast<std::size_t>(u)] = ways1;
            done[static_cast<std::size_t>(u)] = 1U;
        }
    };

    u64 answer = 1ULL;

    for (const auto& cyc : cycles) {
        const int m = static_cast<int>(cyc.size());
        std::vector<u64> w0(static_cast<std::size_t>(m), 1ULL);
        std::vector<u64> w1(static_cast<std::size_t>(m), 1ULL);

        for (int i = 0; i < m; ++i) {
            const int u = cyc[static_cast<std::size_t>(i)];
            u64 a0 = 1ULL;
            u64 a1 = 1ULL;

            for (int e = head[static_cast<std::size_t>(u)]; e != -1;
                 e = next_edge[static_cast<std::size_t>(e)]) {
                const int ch = to[static_cast<std::size_t>(e)];
                if (in_cycle[static_cast<std::size_t>(ch)] != 0U) {
                    continue;
                }
                compute_tree(ch);
                a0 = mul_mod(a0, add_mod(dp0[static_cast<std::size_t>(ch)], dp1[static_cast<std::size_t>(ch)]));
                a1 = mul_mod(a1, dp0[static_cast<std::size_t>(ch)]);
            }

            w0[static_cast<std::size_t>(i)] = a0;
            w1[static_cast<std::size_t>(i)] = a1;
        }

        u64 component_count = 0ULL;
        if (m == 1) {
            component_count = w0[0];
        } else {
            u64 dp_prev0 = w0[0];
            u64 dp_prev1 = 0ULL;
            for (int i = 1; i < m; ++i) {
                const u64 ndp0 = mul_mod(add_mod(dp_prev0, dp_prev1), w0[static_cast<std::size_t>(i)]);
                const u64 ndp1 = mul_mod(dp_prev0, w1[static_cast<std::size_t>(i)]);
                dp_prev0 = ndp0;
                dp_prev1 = ndp1;
            }
            component_count = add_mod(component_count, add_mod(dp_prev0, dp_prev1));

            dp_prev0 = 0ULL;
            dp_prev1 = w1[0];
            for (int i = 1; i < m; ++i) {
                const u64 ndp0 = mul_mod(add_mod(dp_prev0, dp_prev1), w0[static_cast<std::size_t>(i)]);
                const u64 ndp1 = mul_mod(dp_prev0, w1[static_cast<std::size_t>(i)]);
                dp_prev0 = ndp0;
                dp_prev1 = ndp1;
            }
            component_count = add_mod(component_count, dp_prev0);
        }

        answer = mul_mod(answer, component_count);
    }

    return answer;
}

}  // namespace

int main() {
    assert(solve(3) == 35ULL);
    assert(solve(4) == 2118ULL);

    std::cout << solve(20) << '\n';
    return 0;
}

