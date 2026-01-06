#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

using namespace std;
using boost::multiprecision::cpp_int;

namespace {

cpp_int pow10_cppint(int exp) {
    cpp_int result = 1;
    cpp_int base = 10;
    int e = exp;
    while (e > 0) {
        if (e & 1) result *= base;
        e >>= 1;
        if (e) base *= base;
    }
    return result;
}

cpp_int isqrt_cppint(const cpp_int& n) {
    if (n <= 0) return 0;
    unsigned int bit = boost::multiprecision::msb(n);
    cpp_int x = cpp_int(1) << ((bit + 1) / 2);
    for (;;) {
        cpp_int y = (x + n / x) >> 1;
        if (y >= x) return x;
        x = y;
    }
}

long long sum_digits_str(const string& s) {
    long long sum = 0;
    for (char c : s) sum += c - '0';
    return sum;
}

long long sum_digits_sqrt_fraction(int n, int d) {
    if (d <= 0) return 0;
    cpp_int pow10_d = pow10_cppint(d);
    cpp_int A = cpp_int(n) * pow10_d * pow10_d;
    cpp_int M = isqrt_cppint(A);
    cpp_int frac = M % pow10_d;

    string s = frac.convert_to<string>();
    if (static_cast<int>(s.size()) < d) {
        s.insert(0, d - static_cast<int>(s.size()), '0');
    }
    return sum_digits_str(s);
}

void check(const string& name, long long got, long long expected) {
    if (got != expected) {
        cerr << "Validation failed: " << name << " got " << got
             << " expected " << expected << "\n";
        exit(1);
    }
}

}  // namespace

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    struct Job {
        int n;
        int d;
        long long expected;
        string name;
    };

    vector<Job> checks = {
        {2, 10, 31, "S(2,10)"},
        {2, 100, 481, "S(2,100)"},
    };

    vector<long long> results(checks.size(), 0);
    long long final_result = 0;

    vector<thread> workers;
    workers.reserve(checks.size() + 1);

    for (size_t i = 0; i < checks.size(); ++i) {
        workers.emplace_back([&, i]() {
            results[i] = sum_digits_sqrt_fraction(checks[i].n, checks[i].d);
        });
    }
    workers.emplace_back([&]() {
        final_result = sum_digits_sqrt_fraction(13, 1000);
    });

    for (auto& t : workers) t.join();

    for (size_t i = 0; i < checks.size(); ++i) {
        check(checks[i].name, results[i], checks[i].expected);
    }

    cout << final_result << "\n";
    return 0;
}
