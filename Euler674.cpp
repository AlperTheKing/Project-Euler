#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <future>
#include <atomic>
#include <iomanip>
#include <mutex>

// --- Configuration ---
const long long MOD = 1000000000;
const std::string FILENAME = "I-expressions.txt";

// --- Data Structures ---

enum NodeType { VARIABLE, I_OP };

struct Node {
    NodeType type;
    virtual ~Node() = default;
};

using NodePtr = std::shared_ptr<Node>;

struct VarNode : Node {
    std::string name;
    VarNode(std::string n) : name(n) { type = VARIABLE; }
};

struct INode : Node {
    NodePtr left;
    NodePtr right;
    INode(NodePtr l, NodePtr r) : left(l), right(r) { type = I_OP; }
};

// --- Serialization Helper ---
std::string serialize(NodePtr n) {
    if (n->type == VARIABLE) {
        return std::static_pointer_cast<VarNode>(n)->name;
    }
    auto i = std::static_pointer_cast<INode>(n);
    return "I(" + serialize(i->left) + "," + serialize(i->right) + ")";
}

// --- Robust Parser ---

class Parser {
    std::string input;
    size_t pos;

    void skipWhitespace() {
        while (pos < input.size() && isspace(input[pos])) pos++;
    }

    bool consume(char c) {
        skipWhitespace();
        if (pos < input.size() && input[pos] == c) {
            pos++;
            return true;
        }
        return false;
    }

    bool matchString(const std::string& s) {
        skipWhitespace();
        if (input.compare(pos, s.size(), s) == 0) {
            pos += s.size();
            return true;
        }
        return false;
    }

    NodePtr parseTerm() {
        if (matchString("I(") || matchString("J(")) { 
            auto left = parseTerm();
            if (!consume(',')) throw std::runtime_error("Expected ','");
            auto right = parseTerm();
            if (!consume(')')) throw std::runtime_error("Expected ')'");
            return std::make_shared<INode>(left, right);
        } else {
            skipWhitespace();
            std::string name;
            while (pos < input.size() && (isalnum(input[pos]) || input[pos] == '_')) {
                name += input[pos++];
            }
            if (name.empty()) throw std::runtime_error("Unexpected char");
            return std::make_shared<VarNode>(name);
        }
    }

public:
    Parser(const std::string& content) : input(content), pos(0) {}

    NodePtr parseNext() {
        return parseTerm();
    }

    bool hasMore() {
        skipWhitespace();
        return pos < input.size();
    }
};

// --- Unification Engine ---

using Bindings = std::map<std::string, NodePtr>;
struct NodePtrHash { size_t operator()(const NodePtr& p) const { return std::hash<Node*>()(p.get()); } };

NodePtr resolve(NodePtr node, const Bindings& bindings) {
    if (node->type == VARIABLE) {
        auto v = std::static_pointer_cast<VarNode>(node);
        auto it = bindings.find(v->name);
        if (it != bindings.end()) {
            return resolve(it->second, bindings);
        }
    }
    return node;
}

// Occurs check with visited set optimization
bool occurs(const std::string& varName, NodePtr term, const Bindings& bindings, std::unordered_set<NodePtr, NodePtrHash>& visited) {
    NodePtr r = resolve(term, bindings);
    if (visited.count(r)) return false;
    visited.insert(r);

    if (r->type == VARIABLE) {
        return std::static_pointer_cast<VarNode>(r)->name == varName;
    } else {
        auto iNode = std::static_pointer_cast<INode>(r);
        return occurs(varName, iNode->left, bindings, visited) || occurs(varName, iNode->right, bindings, visited);
    }
}

bool unify(NodePtr t1, NodePtr t2, Bindings& bindings) {
    NodePtr r1 = resolve(t1, bindings);
    NodePtr r2 = resolve(t2, bindings);

    if (r1 == r2) return true;

    if (r1->type == VARIABLE) {
        auto v1 = std::static_pointer_cast<VarNode>(r1);
        if (r2->type == VARIABLE && std::static_pointer_cast<VarNode>(r2)->name == v1->name) return true;
        std::unordered_set<NodePtr, NodePtrHash> visited;
        if (occurs(v1->name, r2, bindings, visited)) return false;
        bindings[v1->name] = r2;
        return true;
    }
    if (r2->type == VARIABLE) {
        auto v2 = std::static_pointer_cast<VarNode>(r2);
        std::unordered_set<NodePtr, NodePtrHash> visited;
        if (occurs(v2->name, r1, bindings, visited)) return false;
        bindings[v2->name] = r1;
        return true;
    }

    if (r1->type == I_OP && r2->type == I_OP) {
        auto i1 = std::static_pointer_cast<INode>(r1);
        auto i2 = std::static_pointer_cast<INode>(r2);
        return unify(i1->left, i2->left, bindings) && unify(i1->right, i2->right, bindings);
    }

    return false;
}

// --- Evaluation ---

long long evalJ(long long x, long long y) {
    long long sum = (1 + x + y) % MOD;
    long long term1 = (sum * sum) % MOD;
    long long diff = (y - x) % MOD;
    long long res = (term1 + diff) % MOD;
    if (res < 0) res += MOD;
    return res;
}

using EvalCache = std::map<NodePtr, long long>;

long long evaluate(NodePtr node, const Bindings& bindings, EvalCache& cache) {
    NodePtr r = resolve(node, bindings);
    auto it = cache.find(r);
    if (it != cache.end()) return it->second;

    long long result;
    if (r->type == VARIABLE) {
        result = 0; 
    } else {
        auto iNode = std::static_pointer_cast<INode>(r);
        long long lv = evaluate(iNode->left, bindings, cache);
        long long rv = evaluate(iNode->right, bindings, cache);
        result = evalJ(lv, rv);
    }
    cache[r] = result;
    return result;
}

long long solvePair(NodePtr e1, NodePtr e2) {
    Bindings bindings;
    if (unify(e1, e2, bindings)) {
        EvalCache cache;
        return evaluate(e1, bindings, cache);
    }
    return 0; 
}

// --- Main ---

int main() {
    try {
        std::cout << "--- Project Euler 674 Solver (Final Corrected) ---\n";
        
        std::ifstream file(FILENAME, std::ios::binary | std::ios::ate);
        if (!file.is_open()) return 1;
        
        std::streamsize size = file.tellg();
        file.seekg(0, std::ios::beg);
        std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        
        std::cout << "Parsing file (" << size << " bytes)..." << std::endl;
        Parser fileParser(content);
        
        // 1. Parse ALL expressions, preserving order and duplicates
        std::vector<NodePtr> fileExprs;
        while (fileParser.hasMore()) {
            fileExprs.push_back(fileParser.parseNext());
        }

        size_t totalExprs = fileExprs.size();
        std::cout << "Parsed " << totalExprs << " total expressions from file." << std::endl;

        // 2. Map to unique IDs to handle duplicates efficiently
        std::map<std::string, int> strToId;
        std::vector<NodePtr> uniqueExprs;
        std::vector<int> fileIds; // The file as a list of IDs

        for (auto& e : fileExprs) {
            std::string s = serialize(e);
            if (strToId.find(s) == strToId.end()) {
                strToId[s] = uniqueExprs.size();
                uniqueExprs.push_back(e);
            }
            fileIds.push_back(strToId[s]);
        }
        
        std::cout << "Found " << uniqueExprs.size() << " unique expressions." << std::endl;

        // 3. Pre-compute LSV for all unique pairs (i, j) where i <= j
        // We use a linear map index for pairs
        int nUnique = uniqueExprs.size();
        std::cout << "Pre-computing LSVs for unique pairs..." << std::endl;
        
        // Use a flat vector for the cache: cache[i * n + j]
        std::vector<long long> pairCache(nUnique * nUnique, -1);

        // Parallelize pre-computation
        unsigned int numThreads = std::thread::hardware_concurrency();
        if (numThreads == 0) numThreads = 4;
        
        std::atomic<int> progress(0);
        std::vector<std::future<void>> futures;
        size_t chunkSize = (nUnique + numThreads - 1) / numThreads;

        for (unsigned int t = 0; t < numThreads; ++t) {
            size_t start = t * chunkSize;
            size_t end = std::min(start + chunkSize, (size_t)nUnique);
            if (start >= end) break;

            futures.push_back(std::async(std::launch::async, [&, start, end]() {
                for (size_t i = start; i < end; ++i) {
                    for (size_t j = i; j < nUnique; ++j) { // Compute for j >= i
                         // If i == j, result implies unification of identical trees -> eval at 0.
                         long long val = solvePair(uniqueExprs[i], uniqueExprs[j]);
                         pairCache[i * nUnique + j] = val;
                         pairCache[j * nUnique + i] = val; // Symmetric
                    }
                    if (i % 5 == 0) progress.fetch_add(5);
                }
            }));
        }
        
        for (auto& f : futures) f.get();
        std::cout << "\nPre-computation done." << std::endl;

        // 4. Sum over file positions
        std::cout << "Summing over file pairs..." << std::endl;
        long long totalSum = 0;
        
        // Loop over file positions (i, j)
        for (size_t i = 0; i < totalExprs; ++i) {
            for (size_t j = i + 1; j < totalExprs; ++j) {
                int id1 = fileIds[i];
                int id2 = fileIds[j];
                
                // "pairs made of distinct expressions" -> Skip if same unique expression
                if (id1 == id2) continue;
                
                totalSum = (totalSum + pairCache[id1 * nUnique + id2]) % MOD;
            }
        }

        std::cout << "\nFinal Result (Last 9 digits): " 
                  << std::setfill('0') << std::setw(9) << totalSum << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "\nException: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}