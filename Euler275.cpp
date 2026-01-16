#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

using std::array;
using std::cout;
using std::cerr;
using std::max;
using std::size_t;
using std::string;
using std::uint64_t;
using std::vector;

struct Bitset {
    static constexpr int WORDS = 6;
    uint64_t w[WORDS]{};

    inline void set(int i) { w[i >> 6] |= (1ULL << (i & 63)); }
    inline void reset(int i) { w[i >> 6] &= ~(1ULL << (i & 63)); }
    inline bool test(int i) const { return (w[i >> 6] >> (i & 63)) & 1ULL; }
};

struct Cell {
    int x = 0;
    int y = 0;
    int deg = 0;
    int nbr[4]{};
    int mirror = -1;
};

struct Task {
    Bitset occ, forb, inU;
    array<int, 128> U{};
    int uCount = 0;

    array<int, 32> placed{};
    int k = 0;
    int sumx = 0;
};

class Solver {
  public:
    explicit Solver(int n) : n(n) {
        buildDomain();
        buildMaxCorr();
    }

    uint64_t solve(int threads) {
        Task st = initialState();

        vector<Task> tasks;
        int splitDepth = (n >= 18 ? 8 : max(2, n / 2));
        generateTasks(st, splitDepth, tasks);

        std::atomic<size_t> idx{0};
        vector<std::thread> pool;
        vector<uint64_t> total(threads, 0), sym(threads, 0);

        for (int t = 0; t < threads; ++t) {
            pool.emplace_back([&, t]() {
                while (true) {
                    size_t i = idx.fetch_add(1, std::memory_order_relaxed);
                    if (i >= tasks.size()) break;
                    Task local = tasks[i];
                    dfs(local, total[t], sym[t]);
                }
            });
        }
        for (auto& th : pool) th.join();

        uint64_t N = 0, S = 0;
        for (int t = 0; t < threads; ++t) {
            N += total[t];
            S += sym[t];
        }

        return (N + S) / 2;
    }

  private:
    int n;
    vector<Cell> dom;
    int root = -1;
    vector<int> maxCorr;

    inline bool isBlockCell(int idx) const { return dom[idx].y > 0; }

    void buildDomain() {
        vector<std::pair<int, int>> coords;
        coords.reserve((n + 1) * (n + 1));
        for (int y = 0; y <= n; ++y) {
            int lim = n - y;
            for (int x = -lim; x <= lim; ++x) coords.emplace_back(x, y);
        }

        std::unordered_map<long long, int> id;
        id.reserve(coords.size() * 2);

        auto key = [](int x, int y) -> long long {
            return (static_cast<long long>(static_cast<uint32_t>(x)) << 32) ^
                   static_cast<uint32_t>(y);
        };

        dom.resize(coords.size());
        for (int i = 0; i < static_cast<int>(coords.size()); ++i) {
            auto [x, y] = coords[i];
            dom[i].x = x;
            dom[i].y = y;
            id[key(x, y)] = i;
        }

        for (int i = 0; i < static_cast<int>(dom.size()); ++i) {
            int x = dom[i].x, y = dom[i].y;
            auto addNbr = [&](int nx, int ny) {
                auto it = id.find(key(nx, ny));
                if (it == id.end()) return;
                dom[i].nbr[dom[i].deg++] = it->second;
            };
            addNbr(x + 1, y);
            addNbr(x - 1, y);
            addNbr(x, y + 1);
            addNbr(x, y - 1);

            dom[i].mirror = id[key(-x, y)];
            if (x == 0 && y == 0) root = i;
        }
        if (root < 0) throw std::runtime_error("Root not found");
    }

    void buildMaxCorr() {
        vector<int> xs;
        xs.reserve(dom.size());
        for (const auto& c : dom) {
            if (c.y > 0) xs.push_back(c.x);
        }
        std::sort(xs.begin(), xs.end(), std::greater<int>());

        maxCorr.assign(n + 1, 0);
        long long pref = 0;
        for (int r = 1; r <= n; ++r) {
            pref += xs[r - 1];
            maxCorr[r] = static_cast<int>(pref);
        }
    }

    Task initialState() const {
        Task st;

        for (int i = 0; i < static_cast<int>(dom.size()); ++i) {
            if (dom[i].y == 0 && i != root) st.forb.set(i);
        }

        int firstBlock = -1;
        for (int ni = 0; ni < dom[root].deg; ++ni) {
            int v = dom[root].nbr[ni];
            if (dom[v].x == 0 && dom[v].y == 1) firstBlock = v;
        }
        if (firstBlock < 0) throw std::runtime_error("Missing (0,1)");

        st.occ.set(root);
        st.occ.set(firstBlock);
        st.placed[0] = root;
        st.placed[1] = firstBlock;
        st.k = 1;
        st.sumx = 0;

        for (int ni = 0; ni < dom[firstBlock].deg; ++ni) {
            int v = dom[firstBlock].nbr[ni];
            if (!isBlockCell(v)) continue;
            if (st.occ.test(v) || st.forb.test(v) || st.inU.test(v)) continue;
            st.U[st.uCount++] = v;
            st.inU.set(v);
        }
        return st;
    }

    inline bool isSymmetric(const Task& st) const {
        for (int i = 0; i <= st.k; ++i) {
            int idx = st.placed[i];
            if (!st.occ.test(dom[idx].mirror)) return false;
        }
        return true;
    }

    void dfs(Task& st, uint64_t& total, uint64_t& sym) {
        int remaining = n - st.k;
        if (std::abs(st.sumx) > maxCorr[remaining]) return;

        if (st.k == n) {
            if (st.sumx == 0) {
                total++;
                if (isSymmetric(st)) sym++;
            }
            return;
        }

        int processedCount = 0;
        int processed[128];

        while (st.uCount > 0) {
            int u = st.U[--st.uCount];
            st.inU.reset(u);

            st.occ.set(u);
            st.placed[st.k + 1] = u;

            int addedStart = st.uCount;
            for (int ni = 0; ni < dom[u].deg; ++ni) {
                int v = dom[u].nbr[ni];
                if (!isBlockCell(v)) continue;
                if (st.occ.test(v) || st.forb.test(v) || st.inU.test(v)) continue;
                st.U[st.uCount++] = v;
                st.inU.set(v);
            }

            st.k++;
            st.sumx += dom[u].x;
            dfs(st, total, sym);
            st.sumx -= dom[u].x;
            st.k--;

            while (st.uCount > addedStart) {
                int v = st.U[--st.uCount];
                st.inU.reset(v);
            }
            st.occ.reset(u);

            st.forb.set(u);
            processed[processedCount++] = u;
        }

        for (int i = processedCount - 1; i >= 0; --i) {
            int u = processed[i];
            st.forb.reset(u);
            st.U[st.uCount++] = u;
            st.inU.set(u);
        }
    }

    void generateTasks(Task& st, int splitDepth, vector<Task>& out) {
        int remaining = n - st.k;
        if (std::abs(st.sumx) > maxCorr[remaining]) return;

        if (st.k >= splitDepth) {
            out.push_back(st);
            return;
        }

        int processedCount = 0;
        int processed[128];

        while (st.uCount > 0) {
            int u = st.U[--st.uCount];
            st.inU.reset(u);

            st.occ.set(u);
            st.placed[st.k + 1] = u;

            int addedStart = st.uCount;
            for (int ni = 0; ni < dom[u].deg; ++ni) {
                int v = dom[u].nbr[ni];
                if (!isBlockCell(v)) continue;
                if (st.occ.test(v) || st.forb.test(v) || st.inU.test(v)) continue;
                st.U[st.uCount++] = v;
                st.inU.set(v);
            }

            st.k++;
            st.sumx += dom[u].x;
            generateTasks(st, splitDepth, out);
            st.sumx -= dom[u].x;
            st.k--;

            while (st.uCount > addedStart) {
                int v = st.U[--st.uCount];
                st.inU.reset(v);
            }
            st.occ.reset(u);

            st.forb.set(u);
            processed[processedCount++] = u;
        }

        for (int i = processedCount - 1; i >= 0; --i) {
            int u = processed[i];
            st.forb.reset(u);
            st.U[st.uCount++] = u;
            st.inU.set(u);
        }
    }
};

static void usage(const char* argv0) {
    cerr << "Usage: " << argv0 << " [-n N] [-t THREADS] [--validate]\n";
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int n = 18;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads <= 0) threads = 1;
    bool validate = false;

    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        if (a == "-n" && i + 1 < argc) {
            n = std::stoi(argv[++i]);
        } else if (a == "-t" && i + 1 < argc) {
            threads = max(1, std::stoi(argv[++i]));
        } else if (a == "--validate") {
            validate = true;
        } else {
            usage(argv[0]);
            return 1;
        }
    }

    auto run = [&](int nn) -> uint64_t {
        Solver s(nn);
        return s.solve(threads);
    };

    if (validate) {
        struct Check {
            int n;
            uint64_t expected;
        };
        vector<Check> checks = {
            {6, 18},
            {10, 964},
            {15, 360505},
            {18, 15030564},
        };
        for (const auto& c : checks) {
            uint64_t got = run(c.n);
            if (got != c.expected) {
                cerr << "Validation failed for n=" << c.n << " expected=" << c.expected
                     << " got=" << got << "\n";
                return 2;
            }
        }
    }

    cout << run(n) << "\n";
    return 0;
}
