#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

constexpr std::uint32_t kMod = 100'000'007U;
constexpr int kCells = 16;
constexpr std::uint8_t kAsciiL = 76;
constexpr std::uint8_t kAsciiR = 82;
constexpr std::uint8_t kAsciiU = 85;
constexpr std::uint8_t kAsciiD = 68;

constexpr const char* kStartBoard = "WRBBRRBBRRBBRRBB";
constexpr const char* kTargetBoard = "WBRBBRBRRBRBBRBR";
constexpr const char* kExampleAfterLULUR = "RRBBRBBBRWRBRRBB";

using State = std::uint32_t;

struct Move {
    State next_state = 0;
    std::uint8_t ascii = 0;
};

struct Edge {
    int to = -1;
    std::uint8_t ascii = 0;
};

struct AdjList {
    std::array<Edge, 4> edges{};
    int degree = 0;
};

struct SearchSpace {
    std::vector<State> states;
    std::unordered_map<State, int> id_of;
    std::vector<int> dist_from_start;
};

struct SolverData {
    int start_id = -1;
    int target_id = -1;
    int shortest_len = -1;
    std::vector<int> dist_from_start;
    std::vector<int> dist_to_target;
    std::vector<AdjList> adjacency;
};

struct SolveResult {
    std::uint64_t path_count = 0;
    unsigned __int128 checksum_sum = 0;
};

inline int get_cell(State s, int idx) {
    return static_cast<int>((s >> (2 * idx)) & 3U);
}

inline State set_cell(State s, int idx, int value) {
    const State mask = 3U << (2 * idx);
    s &= ~mask;
    s |= (static_cast<State>(value) << (2 * idx));
    return s;
}

State encode_board(const std::string& board) {
    if (board.size() != kCells) return 0;
    State s = 0;
    for (int i = 0; i < kCells; ++i) {
        int value = 0;
        if (board[static_cast<std::size_t>(i)] == 'R') {
            value = 1;
        } else if (board[static_cast<std::size_t>(i)] == 'B') {
            value = 2;
        }
        s = set_cell(s, i, value);
    }
    return s;
}

bool has_expected_tile_counts(State s) {
    int blank = 0;
    int red = 0;
    int blue = 0;
    for (int i = 0; i < kCells; ++i) {
        const int cell = get_cell(s, i);
        if (cell == 0) {
            ++blank;
        } else if (cell == 1) {
            ++red;
        } else if (cell == 2) {
            ++blue;
        } else {
            return false;
        }
    }
    return blank == 1 && red == 7 && blue == 8;
}

int find_blank(State s) {
    for (int i = 0; i < kCells; ++i) {
        if (get_cell(s, i) == 0) return i;
    }
    return -1;
}

Move make_move(State s, int blank_pos, int other_pos, std::uint8_t ascii) {
    Move mv;
    mv.ascii = ascii;
    const int tile = get_cell(s, other_pos);
    s = set_cell(s, blank_pos, tile);
    s = set_cell(s, other_pos, 0);
    mv.next_state = s;
    return mv;
}

int generate_moves(State s, std::array<Move, 4>& out) {
    const int blank = find_blank(s);
    const int r = blank / 4;
    const int c = blank % 4;
    int count = 0;

    // Tile moves left => blank moves right.
    if (c < 3) out[static_cast<std::size_t>(count++)] = make_move(s, blank, blank + 1, kAsciiL);
    // Tile moves right => blank moves left.
    if (c > 0) out[static_cast<std::size_t>(count++)] = make_move(s, blank, blank - 1, kAsciiR);
    // Tile moves up => blank moves down.
    if (r < 3) out[static_cast<std::size_t>(count++)] = make_move(s, blank, blank + 4, kAsciiU);
    // Tile moves down => blank moves up.
    if (r > 0) out[static_cast<std::size_t>(count++)] = make_move(s, blank, blank - 4, kAsciiD);

    return count;
}

std::uint32_t update_checksum(std::uint32_t checksum, std::uint8_t ascii) {
    return static_cast<std::uint32_t>((static_cast<std::uint64_t>(checksum) * 243ULL + ascii) % kMod);
}

std::uint32_t checksum_for_sequence(const std::string& sequence) {
    std::uint32_t checksum = 0;
    for (char ch : sequence) {
        std::uint8_t ascii = 0;
        if (ch == 'L') ascii = kAsciiL;
        if (ch == 'R') ascii = kAsciiR;
        if (ch == 'U') ascii = kAsciiU;
        if (ch == 'D') ascii = kAsciiD;
        checksum = update_checksum(checksum, ascii);
    }
    return checksum;
}

State apply_sequence(State start, const std::string& sequence) {
    State s = start;
    for (char ch : sequence) {
        const int blank = find_blank(s);
        const int r = blank / 4;
        const int c = blank % 4;
        int other = -1;

        if (ch == 'L' && c < 3) other = blank + 1;
        if (ch == 'R' && c > 0) other = blank - 1;
        if (ch == 'U' && r < 3) other = blank + 4;
        if (ch == 'D' && r > 0) other = blank - 4;

        if (other < 0) return 0;
        const int tile = get_cell(s, other);
        s = set_cell(s, blank, tile);
        s = set_cell(s, other, 0);
    }
    return s;
}

SearchSpace build_reachable_space(State start_state) {
    SearchSpace out;
    out.states.reserve(110'000);
    out.id_of.reserve(140'000);
    out.dist_from_start.reserve(110'000);

    out.id_of.emplace(start_state, 0);
    out.states.push_back(start_state);
    out.dist_from_start.push_back(0);

    std::vector<int> queue;
    queue.reserve(110'000);
    queue.push_back(0);

    std::array<Move, 4> moves{};
    for (std::size_t qpos = 0; qpos < queue.size(); ++qpos) {
        const int v = queue[qpos];
        const State s = out.states[static_cast<std::size_t>(v)];
        const int degree = generate_moves(s, moves);

        for (int i = 0; i < degree; ++i) {
            const State ns = moves[static_cast<std::size_t>(i)].next_state;
            auto it_insert = out.id_of.emplace(ns, static_cast<int>(out.states.size()));
            if (it_insert.second) {
                out.states.push_back(ns);
                out.dist_from_start.push_back(-1);
            }

            const int nid = it_insert.first->second;
            if (out.dist_from_start[static_cast<std::size_t>(nid)] == -1) {
                out.dist_from_start[static_cast<std::size_t>(nid)] =
                    out.dist_from_start[static_cast<std::size_t>(v)] + 1;
                queue.push_back(nid);
            }
        }
    }

    return out;
}

std::vector<AdjList> build_adjacency(const std::vector<State>& states,
                                     const std::unordered_map<State, int>& id_of) {
    std::vector<AdjList> adjacency(states.size());
    std::array<Move, 4> moves{};

    for (std::size_t i = 0; i < states.size(); ++i) {
        const int degree = generate_moves(states[i], moves);
        AdjList& row = adjacency[i];
        row.degree = degree;
        for (int j = 0; j < degree; ++j) {
            const Move mv = moves[static_cast<std::size_t>(j)];
            auto it = id_of.find(mv.next_state);
            row.edges[static_cast<std::size_t>(j)] = {it->second, mv.ascii};
        }
    }

    return adjacency;
}

std::vector<int> bfs_dist(const std::vector<AdjList>& adjacency, int source_id) {
    std::vector<int> dist(adjacency.size(), -1);
    dist[static_cast<std::size_t>(source_id)] = 0;

    std::vector<int> queue;
    queue.reserve(adjacency.size());
    queue.push_back(source_id);

    for (std::size_t qpos = 0; qpos < queue.size(); ++qpos) {
        const int v = queue[qpos];
        const AdjList& row = adjacency[static_cast<std::size_t>(v)];
        const int nd = dist[static_cast<std::size_t>(v)] + 1;

        for (int i = 0; i < row.degree; ++i) {
            const int to = row.edges[static_cast<std::size_t>(i)].to;
            if (dist[static_cast<std::size_t>(to)] == -1) {
                dist[static_cast<std::size_t>(to)] = nd;
                queue.push_back(to);
            }
        }
    }

    return dist;
}

SolverData build_solver_data(State start_state, State target_state) {
    SolverData data;
    SearchSpace space = build_reachable_space(start_state);

    data.start_id = 0;
    auto it_target = space.id_of.find(target_state);
    if (it_target == space.id_of.end()) return data;
    data.target_id = it_target->second;
    data.dist_from_start = std::move(space.dist_from_start);
    data.adjacency = build_adjacency(space.states, space.id_of);
    data.dist_to_target = bfs_dist(data.adjacency, data.target_id);
    data.shortest_len = data.dist_from_start[static_cast<std::size_t>(data.target_id)];
    return data;
}

inline bool is_shortest_dag_edge(const SolverData& data, int from, int to) {
    const int d_to = data.dist_from_start[static_cast<std::size_t>(to)];
    return d_to == data.dist_from_start[static_cast<std::size_t>(from)] + 1 &&
           data.dist_to_target[static_cast<std::size_t>(to)] == data.shortest_len - d_to;
}

void dfs_enumerate(const SolverData& data,
                   int node,
                   int depth,
                   std::uint32_t checksum,
                   std::uint64_t& path_count,
                   unsigned __int128& checksum_sum) {
    if (depth == data.shortest_len) {
        if (node == data.target_id) {
            ++path_count;
            checksum_sum += checksum;
        }
        return;
    }

    const AdjList& row = data.adjacency[static_cast<std::size_t>(node)];
    for (int i = 0; i < row.degree; ++i) {
        const Edge& edge = row.edges[static_cast<std::size_t>(i)];
        const int to = edge.to;
        if (!is_shortest_dag_edge(data, node, to)) continue;
        dfs_enumerate(data,
                      to,
                      depth + 1,
                      update_checksum(checksum, edge.ascii),
                      path_count,
                      checksum_sum);
    }
}

struct PrefixTask {
    int node = -1;
    int depth = 0;
    std::uint32_t checksum = 0;
};

std::vector<PrefixTask> build_prefix_tasks(const SolverData& data, std::size_t desired_tasks) {
    std::vector<PrefixTask> tasks(1, PrefixTask{data.start_id, 0, 0});
    if (desired_tasks <= 1 || data.shortest_len <= 0) return tasks;

    while (tasks.size() < desired_tasks) {
        std::vector<PrefixTask> next;
        next.reserve(tasks.size() * 2);

        for (const PrefixTask& t : tasks) {
            if (t.depth >= data.shortest_len) {
                next.push_back(t);
                continue;
            }

            const AdjList& row = data.adjacency[static_cast<std::size_t>(t.node)];
            for (int i = 0; i < row.degree; ++i) {
                const Edge& edge = row.edges[static_cast<std::size_t>(i)];
                const int to = edge.to;
                if (!is_shortest_dag_edge(data, t.node, to)) continue;
                next.push_back(PrefixTask{
                    to,
                    t.depth + 1,
                    update_checksum(t.checksum, edge.ascii),
                });
            }
        }

        if (next.empty() || next.size() == tasks.size()) break;
        tasks.swap(next);
    }

    return tasks;
}

SolveResult solve_shortest_checksum_sum(const SolverData& data, int threads) {
    SolveResult out;
    if (data.start_id < 0 || data.target_id < 0 || data.shortest_len < 0) return out;

    if (threads < 1) threads = 1;
    const std::size_t target_task_count = static_cast<std::size_t>(threads) * 64U;
    const std::vector<PrefixTask> tasks = build_prefix_tasks(data, target_task_count);

    if (threads == 1 || tasks.size() <= 1) {
        for (const PrefixTask& t : tasks) {
            dfs_enumerate(data, t.node, t.depth, t.checksum, out.path_count, out.checksum_sum);
        }
        return out;
    }

    std::atomic<std::size_t> next_task{0};
    std::vector<std::uint64_t> local_counts(static_cast<std::size_t>(threads), 0);
    std::vector<unsigned __int128> local_sums(static_cast<std::size_t>(threads), 0);
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(threads));

    for (int tid = 0; tid < threads; ++tid) {
        pool.emplace_back([&, tid]() {
            std::uint64_t count = 0;
            unsigned __int128 sum = 0;

            while (true) {
                const std::size_t idx = next_task.fetch_add(1, std::memory_order_relaxed);
                if (idx >= tasks.size()) break;
                const PrefixTask& task = tasks[idx];
                dfs_enumerate(data, task.node, task.depth, task.checksum, count, sum);
            }

            local_counts[static_cast<std::size_t>(tid)] = count;
            local_sums[static_cast<std::size_t>(tid)] = sum;
        });
    }

    for (std::thread& th : pool) th.join();
    for (int tid = 0; tid < threads; ++tid) {
        out.path_count += local_counts[static_cast<std::size_t>(tid)];
        out.checksum_sum += local_sums[static_cast<std::size_t>(tid)];
    }

    return out;
}

std::string to_string_u128(unsigned __int128 value) {
    if (value == 0) return "0";
    std::string out;
    while (value > 0) {
        const unsigned digit = static_cast<unsigned>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

bool validate() {
    const State start_state = encode_board(kStartBoard);
    const State target_state = encode_board(kTargetBoard);
    const State example_state = encode_board(kExampleAfterLULUR);

    if (!has_expected_tile_counts(start_state) || !has_expected_tile_counts(target_state)) {
        std::cerr << "Validation failed: tile counts are invalid.\n";
        return false;
    }

    const std::uint32_t checksum_sample = checksum_for_sequence("LULUR");
    if (checksum_sample != 19'761'398U) {
        std::cerr << "Validation failed: checksum(LULUR) expected 19761398, got "
                  << checksum_sample << "\n";
        return false;
    }

    const State moved = apply_sequence(start_state, "LULUR");
    if (moved != example_state) {
        std::cerr << "Validation failed: LULUR did not reach the official example state.\n";
        return false;
    }

    const SolverData data = build_solver_data(start_state, target_state);
    if (data.shortest_len <= 0) {
        std::cerr << "Validation failed: target is not reachable from start.\n";
        return false;
    }

    const SolveResult single = solve_shortest_checksum_sum(data, 1);
    if (single.path_count == 0) {
        std::cerr << "Validation failed: no shortest path found.\n";
        return false;
    }

    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) hw = 2;
    const int threads = static_cast<int>(std::min<unsigned>(hw, 8));
    if (threads > 1) {
        const SolveResult parallel = solve_shortest_checksum_sum(data, threads);
        if (parallel.path_count != single.path_count || parallel.checksum_sum != single.checksum_sum) {
            std::cerr << "Validation failed: single-thread and multi-thread results differ.\n";
            return false;
        }
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (!validate()) return 1;

    int threads = 0;
    if (argc > 1) {
        threads = std::max(1, std::atoi(argv[1]));
    } else {
        unsigned hw = std::thread::hardware_concurrency();
        if (hw == 0) hw = 2;
        threads = static_cast<int>(std::min<unsigned>(hw, 8));
    }

    const State start_state = encode_board(kStartBoard);
    const State target_state = encode_board(kTargetBoard);
    const SolverData data = build_solver_data(start_state, target_state);
    const SolveResult result = solve_shortest_checksum_sum(data, threads);
    std::cout << to_string_u128(result.checksum_sum) << '\n';
    return 0;
}
