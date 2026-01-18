#include <algorithm>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using u64 = unsigned long long;

struct Node {
    int layer = 0;
    int parent_count = 0;
    std::vector<int> adj;
};

int add_node(std::vector<Node>& graph, int layer, int parent_count) {
    graph.push_back(Node{layer, parent_count, {}});
    return static_cast<int>(graph.size() - 1);
}

void add_edge(std::vector<Node>& graph, int u, int v) {
    graph[u].adj.push_back(v);
    graph[v].adj.push_back(u);
}

std::vector<Node> build_graph(int max_layer) {
    std::vector<Node> graph;
    graph.reserve(80000);

    int root = add_node(graph, 0, 0);
    std::vector<int> boundary(7, root);

    for (int depth = 0; depth < max_layer; ++depth) {
        int n = static_cast<int>(boundary.size());
        std::vector<int> generated(n, -1);
        std::vector<char> consumed(n, 0);

        // Corner tiles: transitions between different parents in the boundary.
        for (int i = 0; i < n; ++i) {
            int j = (i + 1) % n;
            if (consumed[i] || consumed[j]) continue;
            if (boundary[i] != boundary[j]) {
                int nid = add_node(graph, depth + 1, 2);
                add_edge(graph, nid, boundary[i]);
                add_edge(graph, nid, boundary[j]);
                generated[i] = nid;
                generated[j] = nid;
                consumed[i] = 1;
                consumed[j] = 1;
            }
        }

        // Side tiles: the remaining boundary edges.
        for (int i = 0; i < n; ++i) {
            if (consumed[i]) continue;
            int nid = add_node(graph, depth + 1, 1);
            add_edge(graph, nid, boundary[i]);
            generated[i] = nid;
            consumed[i] = 1;
        }

        // New layer order: compress runs of identical generated ids.
        std::vector<int> new_layer;
        new_layer.reserve(n);
        for (int i = 0; i < n; ++i) {
            if (i == 0 || generated[i] != generated[i - 1]) {
                new_layer.push_back(generated[i]);
            }
        }
        if (new_layer.size() > 1 && new_layer.front() == new_layer.back()) {
            new_layer.pop_back();
        }

        // Close the cycle of tiles within the layer.
        int m = static_cast<int>(new_layer.size());
        for (int i = 0; i < m; ++i) {
            add_edge(graph, new_layer[i], new_layer[(i + 1) % m]);
        }

        // Build the next boundary: side tiles contribute 4 outer edges, corners 3.
        std::vector<int> next_boundary;
        next_boundary.reserve(m * 4);
        for (int id : new_layer) {
            int out_edges = (graph[id].parent_count == 1) ? 4 : 3;
            for (int k = 0; k < out_edges; ++k) {
                next_boundary.push_back(id);
            }
        }
        boundary.swap(next_boundary);
    }

    return graph;
}

void step_counts(const std::vector<Node>& graph,
                 const std::vector<u64>& cur,
                 std::vector<u64>& next,
                 unsigned int threads) {
    size_t n = graph.size();
    next.assign(n, 0);
    if (n == 0) return;

    const size_t kMinParallel = 20000;
    if (threads <= 1 || n < kMinParallel) {
        for (size_t i = 0; i < n; ++i) {
            u64 sum = 0;
            for (int v : graph[i].adj) sum += cur[v];
            next[i] = sum;
        }
        return;
    }

    threads = std::min<unsigned int>(threads, static_cast<unsigned int>(n));
    size_t block = (n + threads - 1) / threads;
    std::vector<std::thread> workers;
    workers.reserve(threads);
    for (unsigned int t = 0; t < threads; ++t) {
        size_t start = t * block;
        if (start >= n) break;
        size_t end = std::min(n, start + block);
        workers.emplace_back([&graph, &cur, &next, start, end]() {
            for (size_t i = start; i < end; ++i) {
                u64 sum = 0;
                for (int v : graph[i].adj) sum += cur[v];
                next[i] = sum;
            }
        });
    }
    for (auto& th : workers) th.join();
}

bool validate_graph(const std::vector<Node>& graph, int max_layer) {
    for (size_t i = 0; i < graph.size(); ++i) {
        if (graph[i].layer < max_layer && graph[i].adj.size() != 7) {
            std::cerr << "Validation failed: node " << i
                      << " in layer " << graph[i].layer
                      << " has degree " << graph[i].adj.size() << "\n";
            return false;
        }
    }
    return true;
}

bool check_value(int step, u64 got, u64 expected) {
    std::cout << "Validation F(" << step << ") = " << got;
    if (got == expected) {
        std::cout << " [OK]\n";
        return true;
    }
    std::cout << " [FAIL]\n";
    return false;
}

int main(int argc, char** argv) {
    int steps = 20;
    if (argc >= 2) {
        steps = std::max(0, std::atoi(argv[1]));
    }

    const int validate_steps = 4;
    int build_steps = std::max(steps, validate_steps);
    int max_layer = build_steps / 2;

    unsigned int threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 4;

    std::cout << "Building graph up to distance " << max_layer << "...\n";
    auto graph = build_graph(max_layer);
    std::cout << "Nodes in ball: " << graph.size() << "\n";

    if (!validate_graph(graph, max_layer)) {
        return 1;
    }

    std::cout << "Running validation checkpoints...\n";
    std::vector<u64> cur(graph.size(), 0);
    std::vector<u64> next;
    cur[0] = 1;

    bool ok = true;
    u64 answer = (steps == 0) ? 1 : 0;

    for (int step = 1; step <= build_steps; ++step) {
        step_counts(graph, cur, next, threads);
        cur.swap(next);

        if (step == steps) answer = cur[0];
        if (step == 2) ok &= check_value(step, cur[0], 7ULL);
        if (step == 3) ok &= check_value(step, cur[0], 14ULL);
        if (step == 4) ok &= check_value(step, cur[0], 119ULL);
    }

    if (!ok) {
        std::cerr << "Validation failed; aborting.\n";
        return 1;
    }

    std::cout << "F(" << steps << ") = " << answer << "\n";
    return 0;
}
