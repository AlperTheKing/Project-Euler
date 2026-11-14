#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using Edge = std::pair<unsigned, unsigned>;
using Permutation = std::vector<unsigned>;
using Partition = std::vector<u32>;
constexpr u64 FINISHED = u64{1} << 63;
constexpr double EPSILON = 1e-8;

void require(const bool condition, const char* description) {
    if (!condition) throw std::runtime_error(std::string("Check failed: ") + description);
}

unsigned first_bit(const u32 mask) {
    return static_cast<unsigned>(__builtin_ctz(mask));
}

struct Point {
    double x, y, z;

    Point operator+(const Point& other) const { return {x + other.x, y + other.y, z + other.z}; }
    Point operator-(const Point& other) const { return {x - other.x, y - other.y, z - other.z}; }
    Point operator*(const double scale) const { return {x * scale, y * scale, z * scale}; }
};

double dot(const Point& a, const Point& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Point cross(const Point& a, const Point& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

Point unit(const Point& point) {
    return point * (1 / std::sqrt(dot(point, point)));
}

std::array<Point, 3> frame(const Point& a, const Point& b) {
    const Point x = unit(a);
    const Point y = unit(b - x * dot(x, b));
    return {x, y, cross(x, y)};
}

std::vector<Point> cube() {
    std::vector<Point> result;
    for (const double x : {-1.0, 1.0}) {
        for (const double y : {-1.0, 1.0}) {
            for (const double z : {-1.0, 1.0}) result.push_back({x, y, z});
        }
    }
    return result;
}

std::vector<Point> dodecahedron() {
    auto result = cube();
    const double phi = (1 + std::sqrt(5.0)) / 2;
    for (const double a : {-1.0, 1.0}) {
        for (const double b : {-1.0, 1.0}) {
            result.push_back({0, a / phi, b * phi});
            result.push_back({a / phi, b * phi, 0});
            result.push_back({b * phi, 0, a / phi});
        }
    }
    return result;
}

struct Graph {
    unsigned vertices;
    std::vector<Edge> edges;
    std::vector<u32> adjacent;

    Graph(const unsigned size, std::vector<Edge> pairs) : vertices(size), edges(std::move(pairs)) {
        require(vertices <= 20 && edges.size() <= 30, "graph fits the masks");
        adjacent.resize(edges.size());
        for (unsigned i = 0; i < edges.size(); ++i) {
            const auto [a, b] = edges[i];
            require(a < vertices && b < vertices && a != b, "valid edge endpoints");
            for (unsigned j = 0; j < edges.size(); ++j) {
                const auto [c, d] = edges[j];
                if (a == c || a == d || b == c || b == d) adjacent[i] |= u32{1} << j;
            }
        }
    }

    bool connected(const u32 selected) const {
        if (selected == 0) return false;
        u32 reached = u32{1} << first_bit(selected);
        u32 pending = reached;
        while (pending != 0) {
            const unsigned edge = first_bit(pending);
            pending &= pending - 1;
            const u32 next = adjacent[edge] & selected & ~reached;
            reached |= next;
            pending |= next;
        }
        return reached == selected;
    }
};

Graph make_graph(const std::vector<Point>& points) {
    double distance = std::numeric_limits<double>::max();
    for (unsigned i = 0; i < points.size(); ++i) {
        for (unsigned j = i + 1; j < points.size(); ++j) {
            const Point delta = points[i] - points[j];
            distance = std::min(distance, dot(delta, delta));
        }
    }
    std::vector<Edge> edges;
    for (unsigned i = 0; i < points.size(); ++i) {
        for (unsigned j = i + 1; j < points.size(); ++j) {
            const Point delta = points[i] - points[j];
            if (std::abs(dot(delta, delta) - distance) < EPSILON) edges.emplace_back(i, j);
        }
    }
    return Graph(static_cast<unsigned>(points.size()), std::move(edges));
}

std::vector<Permutation> rotations(const std::vector<Point>& points, const Graph& graph) {
    const auto [a, b] = graph.edges.front();
    const auto source = frame(points[a], points[b]);
    std::set<Permutation> unique;
    for (const auto& edge : graph.edges) {
        for (unsigned reverse = 0; reverse < 2; ++reverse) {
            const auto target = reverse == 0 ? frame(points[edge.first], points[edge.second])
                                             : frame(points[edge.second], points[edge.first]);
            Permutation vertices;
            for (const Point& point : points) {
                const Point mapped = target[0] * dot(source[0], point)
                                   + target[1] * dot(source[1], point)
                                   + target[2] * dot(source[2], point);
                unsigned match = graph.vertices;
                for (unsigned i = 0; i < points.size(); ++i) {
                    const Point delta = mapped - points[i];
                    if (dot(delta, delta) < EPSILON * EPSILON) {
                        require(match == graph.vertices, "unique rotated vertex");
                        match = i;
                    }
                }
                require(match < graph.vertices, "rotation preserves the vertex set");
                vertices.push_back(match);
            }
            auto sorted = vertices;
            std::sort(sorted.begin(), sorted.end());
            for (unsigned i = 0; i < sorted.size(); ++i) require(sorted[i] == i, "vertex permutation");
            Permutation permutation;
            for (const auto& [u, v] : graph.edges) {
                const Edge mapped = std::minmax(vertices[u], vertices[v]);
                const auto found = std::find(graph.edges.begin(), graph.edges.end(), mapped);
                require(found != graph.edges.end(), "rotation preserves edges");
                permutation.push_back(static_cast<unsigned>(found - graph.edges.begin()));
            }
            require(unique.insert(permutation).second, "distinct directed edges give distinct rotations");
        }
    }
    return {unique.begin(), unique.end()};
}

std::vector<Edge> sweep_edges(const Graph& graph) {
    std::vector<unsigned> order, rank(graph.vertices, graph.vertices);
    for (unsigned start = 0; start < graph.vertices; ++start) {
        if (rank[start] != graph.vertices) continue;
        rank[start] = static_cast<unsigned>(order.size());
        order.push_back(start);
        for (std::size_t i = rank[start]; i < order.size(); ++i) {
            for (const auto& [u, v] : graph.edges) {
                if (u != order[i] && v != order[i]) continue;
                const unsigned next = u == order[i] ? v : u;
                if (rank[next] == graph.vertices) {
                    rank[next] = static_cast<unsigned>(order.size());
                    order.push_back(next);
                }
            }
        }
    }
    auto result = graph.edges;
    const auto key = [&](const Edge& edge) -> Edge { return std::minmax(rank[edge.first], rank[edge.second]); };
    std::sort(result.begin(), result.end(), [&](const Edge& a, const Edge& b) { return key(a) < key(b); });
    return result;
}

u64 connected_subsets(const Graph& graph) {
    const auto edges = sweep_edges(graph);
    std::vector<unsigned> last(graph.vertices), frontier;
    for (unsigned i = 0; i < edges.size(); ++i) {
        last[edges[i].first] = last[edges[i].second] = i;
    }
    std::unordered_map<u64, u64> current{{0, 1}}, next;
    for (unsigned step = 0; step < edges.size(); ++step) {
        const unsigned previous_size = static_cast<unsigned>(frontier.size());
        const auto position = [&](const unsigned vertex) {
            const auto found = std::find(frontier.begin(), frontier.end(), vertex);
            if (found != frontier.end()) return static_cast<unsigned>(found - frontier.begin());
            frontier.push_back(vertex);
            return static_cast<unsigned>(frontier.size() - 1);
        };
        const unsigned u = position(edges[step].first), v = position(edges[step].second);
        require(frontier.size() <= 15, "frontier fits four-bit component labels");
        std::vector<unsigned> keep;
        for (unsigned i = 0; i < frontier.size(); ++i) {
            if (last[frontier[i]] > step) keep.push_back(i);
        }
        next.clear();
        next.reserve(current.size() * 2);
        for (const auto& [state, count] : current) {
            for (unsigned selected = 0; selected < 2; ++selected) {
                if (selected != 0 && (state & FINISHED) != 0) continue;
                std::array<unsigned, 15> labels{};
                unsigned maximum = 0;
                for (unsigned i = 0; i < previous_size; ++i) {
                    labels[i] = static_cast<unsigned>((state >> (4 * i)) & 15);
                    maximum = std::max(maximum, labels[i]);
                }
                if (selected != 0) {
                    if (labels[u] == 0 && labels[v] == 0) {
                        labels[u] = labels[v] = maximum + 1;
                    } else if (labels[u] == 0) {
                        labels[u] = labels[v];
                    } else if (labels[v] == 0) {
                        labels[v] = labels[u];
                    } else {
                        const unsigned from = labels[v], to = labels[u];
                        for (unsigned& label : labels) {
                            if (label == from) label = to;
                        }
                    }
                }
                u32 used = 0, retained = 0;
                for (const unsigned label : labels) used |= u32{1} << label;
                for (const unsigned i : keep) retained |= u32{1} << labels[i];
                used &= ~u32{1};
                retained &= ~u32{1};
                const u32 closed = used & ~retained;
                u64 encoded = state & FINISHED;
                if (closed != 0) {
                    // A component leaving the frontier can never connect to a later edge.
                    if (retained != 0 || (closed & (closed - 1)) != 0) continue;
                    encoded = FINISHED;
                }
                std::array<unsigned, 16> canonical{};
                unsigned next_label = 1;
                for (unsigned i = 0; i < keep.size(); ++i) {
                    const unsigned label = labels[keep[i]];
                    if (label != 0 && canonical[label] == 0) canonical[label] = next_label++;
                    encoded |= u64{canonical[label]} << (4 * i);
                }
                next[encoded] += count;
            }
        }
        std::vector<unsigned> remaining;
        for (const unsigned i : keep) remaining.push_back(frontier[i]);
        frontier.swap(remaining);
        current.swap(next);
    }
    const auto found = current.find(FINISHED);
    return found == current.end() ? 0 : found->second;
}

u32 transform(u32 mask, const Permutation& permutation) {
    u32 result = 0;
    while (mask != 0) {
        result |= u32{1} << permutation[first_bit(mask)];
        mask &= mask - 1;
    }
    return result;
}

Partition cycles(const Permutation& permutation) {
    Partition result;
    u32 visited = 0;
    for (unsigned i = 0; i < permutation.size(); ++i) {
        if ((visited & (u32{1} << i)) != 0) continue;
        u32 cycle = 0;
        unsigned j = i;
        do {
            cycle |= u32{1} << j;
            j = permutation[j];
        } while (j != i);
        visited |= cycle;
        result.push_back(cycle);
    }
    std::sort(result.begin(), result.end());
    return result;
}

Partition canonical_partition(const Partition& partition, const std::vector<Permutation>& group) {
    Partition result = partition;
    for (const auto& permutation : group) {
        Partition candidate;
        for (const u32 cycle : partition) candidate.push_back(transform(cycle, permutation));
        std::sort(candidate.begin(), candidate.end());
        result = std::min(result, candidate);
    }
    return result;
}

u64 fixed_subsets(const Graph& graph, const Partition& partition) {
    if (partition.size() == graph.edges.size()) return connected_subsets(graph);
    require(partition.size() <= 20, "small nonidentity cycle enumeration");
    std::vector<u32> subsets(std::size_t{1} << partition.size());
    u64 result = 0;
    for (u32 mask = 1; mask < subsets.size(); ++mask) {
        subsets[mask] = subsets[mask & (mask - 1)] | partition[first_bit(mask)];
        result += graph.connected(subsets[mask]);
    }
    return result;
}

u64 skeletons(const Graph& graph, const std::vector<Permutation>& group) {
    std::map<Partition, unsigned> classes;
    for (const auto& permutation : group) ++classes[canonical_partition(cycles(permutation), group)];
    u64 total = 0;
    for (const auto& [partition, multiplicity] : classes) total += multiplicity * fixed_subsets(graph, partition);
    require(total % group.size() == 0, "integral Burnside average");
    return total / group.size();
}

u64 brute_connected(const Graph& graph) {
    u64 result = 0;
    for (u32 mask = 1; mask < (u32{1} << graph.edges.size()); ++mask) result += graph.connected(mask);
    return result;
}

u64 brute_skeletons(const Graph& graph, const std::vector<Permutation>& group) {
    std::set<u32> representatives;
    for (u32 mask = 1; mask < (u32{1} << graph.edges.size()); ++mask) {
        if (!graph.connected(mask)) continue;
        u32 representative = mask;
        for (const auto& permutation : group) representative = std::min(representative, transform(mask, permutation));
        representatives.insert(representative);
    }
    return representatives.size();
}

void check_group(const std::vector<Permutation>& group) {
    const std::set<Permutation> elements(group.begin(), group.end());
    for (const auto& a : group) {
        for (const auto& b : group) {
            Permutation product(a.size());
            for (unsigned i = 0; i < a.size(); ++i) product[i] = a[b[i]];
            require(elements.count(product) == 1, "rotation group closure");
        }
    }
}

void run_tests() {
    std::vector<Edge> complete;
    for (unsigned u = 0; u < 5; ++u) {
        for (unsigned v = u + 1; v < 5; ++v) complete.emplace_back(u, v);
    }
    for (u32 mask = 0; mask < (u32{1} << complete.size()); ++mask) {
        std::vector<Edge> edges;
        for (unsigned i = 0; i < complete.size(); ++i) {
            if ((mask & (u32{1} << i)) != 0) edges.push_back(complete[i]);
        }
        const Graph graph(5, std::move(edges));
        require(connected_subsets(graph) == brute_connected(graph), "all 1024 labelled graphs on five vertices");
    }
    const std::vector<Point> tetrahedron{{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
    const std::vector<Point> octahedron{{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    for (const auto& points : {tetrahedron, cube(), octahedron}) {
        const Graph graph = make_graph(points);
        const auto group = rotations(points, graph);
        require(group.size() == (points.size() == 4 ? 12 : 24), "small rotation group size");
        check_group(group);
        require(connected_subsets(graph) == brute_connected(graph), "small polyhedron connected subsets");
        const u64 answer = skeletons(graph, group);
        require(answer == brute_skeletons(graph, group), "direct small polyhedron orbit enumeration");
        if (points.size() == 4) require(answer == 10, "given tetrahedron example");
    }
    const auto points = dodecahedron();
    const Graph graph = make_graph(points);
    const auto group = rotations(points, graph);
    require(graph.vertices == 20 && graph.edges.size() == 30 && group.size() == 60, "dodecahedron dimensions");
    check_group(group);
    std::map<unsigned, unsigned> cycle_counts;
    u64 total = 0;
    for (const auto& permutation : group) {
        const auto partition = cycles(permutation);
        ++cycle_counts[static_cast<unsigned>(partition.size())];
        total += fixed_subsets(graph, partition);
    }
    require(cycle_counts == std::map<unsigned, unsigned>{{6, 24}, {10, 20}, {16, 15}, {30, 1}}, "rotation cycle types");
    require(total == group.size() * skeletons(graph, group), "all 60 rotations versus grouped enumeration");
    auto relabelled = graph.edges;
    for (auto& [u, v] : relabelled) {
        u = graph.vertices - 1 - u;
        v = graph.vertices - 1 - v;
    }
    std::reverse(relabelled.begin(), relabelled.end());
    require(connected_subsets(Graph(graph.vertices, std::move(relabelled))) == connected_subsets(graph),
            "independent frontier order");
    std::cout << "All checks passed: 1024 small graphs, tetrahedron example, direct tetrahedron/cube/octahedron "
                 "orbits, rotation groups, all 60 dodecahedron rotations, and frontier reordering.\n";
}

}

int main(int argc, char* argv[]) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            run_tests();
        } else {
            if (argc != 1) throw std::invalid_argument("Usage: Euler1019 [--self-test]");
            const auto points = dodecahedron();
            const Graph graph = make_graph(points);
            const auto group = rotations(points, graph);
            require(graph.vertices == 20 && graph.edges.size() == 30 && group.size() == 60, "dodecahedron dimensions");
            std::cout << skeletons(graph, group) << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
