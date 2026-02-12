#include <cassert>
#include <cctype>
#include <iostream>
#include <string>
#include <vector>

int main() {
    const std::vector<std::string> words = {
        "affine", "plane",
        "radically", "integral", "local", "field",
        "open", "oriented", "line", "section",
        "jacobian",
        "orthogonal", "kernel", "embedding"
    };

    std::string answer;
    answer.reserve(words.size());
    for (const std::string& w : words) {
        assert(!w.empty());
        answer.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(w.front()))));
    }

    std::cout << answer << '\n';
    return 0;
}
