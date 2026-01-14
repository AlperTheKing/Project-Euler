#include <algorithm>
#include <array>
#include <cstdint>
#include <future>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

namespace {

constexpr int64_t kMod = 10000000000LL;

using Pairing = vector<pair<int, int>>;
using State = vector<int>;

struct StateHash {
    size_t operator()(const State& s) const noexcept {
        size_t h = s.size();
        for (int v : s) {
            h ^= static_cast<size_t>(v + 1) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};

void generate_pairings_rec(const vector<int>& indices, vector<Pairing>& out) {
    if (indices.empty()) {
        out.push_back({});
        return;
    }
    int first = indices[0];
    for (size_t k = 1; k < indices.size(); k += 2) {
        int second = indices[k];
        vector<int> left(indices.begin() + 1, indices.begin() + static_cast<long>(k));
        vector<int> right(indices.begin() + static_cast<long>(k) + 1, indices.end());
        vector<Pairing> left_pairings;
        vector<Pairing> right_pairings;
        generate_pairings_rec(left, left_pairings);
        generate_pairings_rec(right, right_pairings);
        for (const auto& lp : left_pairings) {
            for (const auto& rp : right_pairings) {
                Pairing combined;
                combined.reserve(1 + lp.size() + rp.size());
                combined.push_back({first, second});
                combined.insert(combined.end(), lp.begin(), lp.end());
                combined.insert(combined.end(), rp.begin(), rp.end());
                out.push_back(std::move(combined));
            }
        }
    }
}

vector<Pairing> generate_pairings(int count) {
    vector<Pairing> out;
    if (count == 0) {
        out.push_back({});
        return out;
    }
    vector<int> indices(count);
    for (int i = 0; i < count; ++i) indices[i] = i;
    generate_pairings_rec(indices, out);
    return out;
}

struct PairingCache {
    vector<Pairing> p2;
    vector<Pairing> p4;
    vector<Pairing> p8;

    PairingCache() {
        p2 = generate_pairings(2);
        p4 = generate_pairings(4);
        p8 = generate_pairings(8);
    }

    const vector<Pairing>& get(int count) const {
        if (count == 2) return p2;
        if (count == 4) return p4;
        return p8;
    }
};

enum Port {
    E_L = 0,
    E_U = 1,
    N_R = 2,
    N_L = 3,
    W_U = 4,
    W_L = 5,
    S_L = 6,
    S_R = 7
};

inline void add_mod(int64_t& a, int64_t b) {
    a += b;
    if (a >= kMod) a -= kMod;
}

bool apply_pair(int a, int b, int old_count, vector<int>& partner, int& closed) {
    bool a_old = (a < old_count);
    bool b_old = (b < old_count);

    if (a_old && b_old) {
        int pa = partner[a];
        int pb = partner[b];
        if (pa < 0 || pb < 0) return false;
        if (pa == b) {
            partner[a] = -1;
            partner[b] = -1;
            ++closed;
            return true;
        }
        partner[a] = -1;
        partner[pa] = -1;
        partner[b] = -1;
        partner[pb] = -1;
        partner[pa] = pb;
        partner[pb] = pa;
        return true;
    }

    if (a_old != b_old) {
        int old = a_old ? a : b;
        int neu = a_old ? b : a;
        int pold = partner[old];
        if (pold < 0) return false;
        partner[old] = -1;
        partner[pold] = -1;
        partner[pold] = neu;
        partner[neu] = pold;
        return true;
    }

    if (partner[a] != -1 || partner[b] != -1) return false;
    partner[a] = b;
    partner[b] = a;
    return true;
}

struct VertexContext {
    int old_count;
    int offset;
    int num_bottom;
    int num_left;
    int total_labels;
    bool is_final;
    array<int, 8> port_label;
    vector<int> active_ports;
    vector<int> outgoing_ports;
    vector<int> new_labels;
    vector<int> label_pos;
    const vector<Pairing>* local_pairings;
};

void process_state(const State& state, int64_t count, const VertexContext& ctx,
                   unordered_map<State, int64_t, StateHash>& out) {
    if (static_cast<int>(state.size()) != ctx.old_count) return;
    vector<int> base_partner(ctx.total_labels, -1);
    for (int i = 0; i < ctx.old_count; ++i) base_partner[i] = state[i];

    for (const auto& pairing : *ctx.local_pairings) {
        vector<int> partner = base_partner;
        int closed = 0;
        bool ok = true;
        for (const auto& pr : pairing) {
            int port_a = ctx.active_ports[pr.first];
            int port_b = ctx.active_ports[pr.second];
            int la = ctx.port_label[port_a];
            int lb = ctx.port_label[port_b];
            if (la < 0 || lb < 0) {
                ok = false;
                break;
            }
            if (!apply_pair(la, lb, ctx.old_count, partner, closed)) {
                ok = false;
                break;
            }
            if (closed > 1) {
                ok = false;
                break;
            }
        }
        if (!ok) continue;

        if (ctx.is_final) {
            if (closed != 1 || !ctx.new_labels.empty()) continue;
            State empty;
            auto it = out.find(empty);
            if (it == out.end()) out.emplace(empty, count);
            else add_mod(it->second, count);
            continue;
        }
        if (closed != 0) continue;

        State next_state(ctx.new_labels.size(), -1);
        bool valid = true;
        for (size_t i = 0; i < ctx.new_labels.size(); ++i) {
            int label = ctx.new_labels[i];
            int pl = partner[label];
            if (pl < 0) {
                valid = false;
                break;
            }
            int pos = ctx.label_pos[pl];
            if (pos < 0) {
                valid = false;
                break;
            }
            next_state[i] = pos;
        }
        if (!valid) continue;

        auto it = out.find(next_state);
        if (it == out.end()) out.emplace(std::move(next_state), count);
        else add_mod(it->second, count);
    }
}

int64_t compute_L(int m, int n, bool use_threads) {
    if (n > m) swap(n, m); // Rotate to minimize boundary width.
    PairingCache cache;

    vector<int> row_wires(n + 1, 0);
    for (int k = 0; k <= n; ++k) {
        row_wires[k] = (k > 0) + (k < n);
    }
    vector<int> prefix(n + 2, 0);
    for (int k = 0; k <= n; ++k) {
        prefix[k + 1] = prefix[k] + row_wires[k];
    }
    int total_row_wires = prefix[n + 1];

    unordered_map<State, int64_t, StateHash> dp;
    dp[State()] = 1;

    unsigned hw = thread::hardware_concurrency();
    if (hw == 0) hw = 2;

    for (int c = 0; c <= m; ++c) {
        for (int r = 0; r <= n; ++r) {
            bool c_ne = (c < m && r < n);
            bool c_nw = (c > 0 && r < n);
            bool c_se = (c < m && r > 0);
            bool c_sw = (c > 0 && r > 0);

            vector<int> active_ports;
            active_ports.reserve(8);
            if (c_se) active_ports.push_back(E_L);
            if (c_ne) active_ports.push_back(E_U);
            if (c_ne) active_ports.push_back(N_R);
            if (c_nw) active_ports.push_back(N_L);
            if (c_nw) active_ports.push_back(W_U);
            if (c_sw) active_ports.push_back(W_L);
            if (c_sw) active_ports.push_back(S_L);
            if (c_se) active_ports.push_back(S_R);

            const auto& local_pairings = cache.get(static_cast<int>(active_ports.size()));

            int sum_right = (c < m) ? prefix[r] : 0;
            int num_bottom = (r > 0) ? ((c < m) + (c > 0)) : 0;
            int num_left = (c > 0) ? ((r > 0) + (r < n)) : 0;
            int sum_left = (c > 0) ? (total_row_wires - prefix[r]) : 0;
            int expected_B = sum_right + num_bottom + sum_left;
            int offset = sum_right;

            vector<int> outgoing_ports;
            outgoing_ports.reserve(4);
            if (c_se) outgoing_ports.push_back(E_L);
            if (c_ne) outgoing_ports.push_back(E_U);
            if (c_ne) outgoing_ports.push_back(N_R);
            if (c_nw) outgoing_ports.push_back(N_L);

            int num_out = static_cast<int>(outgoing_ports.size());
            int total_labels = expected_B + num_out;

            array<int, 8> port_label;
            port_label.fill(-1);

            int idx = offset;
            if (c_se) port_label[S_R] = idx++;
            if (c_sw) port_label[S_L] = idx++;

            idx = offset + num_bottom;
            if (c_sw) port_label[W_L] = idx++;
            if (c_nw) port_label[W_U] = idx++;

            for (int i = 0; i < num_out; ++i) {
                port_label[outgoing_ports[i]] = expected_B + i;
            }

            vector<int> new_labels;
            new_labels.reserve(expected_B - num_bottom - num_left + num_out);
            for (int i = 0; i < offset; ++i) new_labels.push_back(i);
            for (int port : outgoing_ports) new_labels.push_back(port_label[port]);
            int remove_end = offset + num_bottom + num_left;
            for (int i = remove_end; i < expected_B; ++i) new_labels.push_back(i);

            vector<int> label_pos(total_labels, -1);
            for (size_t i = 0; i < new_labels.size(); ++i) label_pos[new_labels[i]] = static_cast<int>(i);

            VertexContext ctx;
            ctx.old_count = expected_B;
            ctx.offset = offset;
            ctx.num_bottom = num_bottom;
            ctx.num_left = num_left;
            ctx.total_labels = total_labels;
            ctx.is_final = (c == m && r == n);
            ctx.port_label = port_label;
            ctx.active_ports = active_ports;
            ctx.outgoing_ports = outgoing_ports;
            ctx.new_labels = new_labels;
            ctx.label_pos = label_pos;
            ctx.local_pairings = &local_pairings;

            unordered_map<State, int64_t, StateHash> next;
            size_t state_count = dp.size();
            if (use_threads && state_count > 2000 && hw > 1) {
                vector<pair<State, int64_t>> items;
                items.reserve(state_count);
                for (const auto& kv : dp) items.push_back(kv);

                size_t threads = min<size_t>(hw, items.size());
                size_t chunk = (items.size() + threads - 1) / threads;
                vector<future<unordered_map<State, int64_t, StateHash>>> futures;
                futures.reserve(threads);
                for (size_t t = 0; t < threads; ++t) {
                    size_t start = t * chunk;
                    if (start >= items.size()) break;
                    size_t end = min(items.size(), start + chunk);
                    futures.push_back(async(launch::async, [start, end, &items, &ctx]() {
                        unordered_map<State, int64_t, StateHash> local;
                        for (size_t i = start; i < end; ++i) {
                            process_state(items[i].first, items[i].second, ctx, local);
                        }
                        return local;
                    }));
                }
                for (auto& fut : futures) {
                    auto local = fut.get();
                    for (auto& kv : local) {
                        auto it = next.find(kv.first);
                        if (it == next.end()) next.emplace(std::move(kv.first), kv.second);
                        else add_mod(it->second, kv.second);
                    }
                }
            } else {
                for (const auto& kv : dp) {
                    process_state(kv.first, kv.second, ctx, next);
                }
            }

            dp.swap(next);
        }
    }

    auto it = dp.find(State());
    if (it == dp.end()) return 0;
    return it->second % kMod;
}

} // namespace

int main() {
    struct Check { int m; int n; int64_t expected; };
    vector<Check> checks = {
        {1, 1, 1},
        {1, 2, 2},
        {2, 2, 37},
        {3, 3, 104290},
    };

    for (const auto& chk : checks) {
        int64_t got = compute_L(chk.m, chk.n, false);
        if (got != chk.expected) {
            cerr << "Validation failed for L(" << chk.m << "," << chk.n << "): got "
                 << got << ", expected " << chk.expected << "\n";
            return 1;
        }
    }

    int64_t answer = compute_L(6, 10, true);
    cout << answer << "\n";
    return 0;
}
