#include <bits/stdc++.h>
using namespace std;

static constexpr uint64_t MOD = 888888883ULL;
static constexpr uint64_t MUL = 8888ULL;
static constexpr uint64_t A0  = 88888888ULL;

// Compute F(N) efficiently.
static uint64_t compute_F(uint32_t N) {
    // cnt[mask][invParity]
    // mask = (px<<2)|(py<<1)|pz, where px/py/pz are parities of counts of x,y,z.
    // invParity is inversion parity of the 50-letter block in order x<y<z.
    uint64_t cnt[8][2] = {};

    uint64_t a = A0; // a_0

    for (uint32_t i = 0; i < N; ++i) {
        uint8_t px = 0, py = 0, pz = 0;
        uint8_t inv = 0;

        // Scan 50 symbols b_{50i}..b_{50i+49} (note b_n = a_n mod 3).
        for (int k = 0; k < 50; ++k) {
            uint8_t b = static_cast<uint8_t>(a % 3ULL);
            if (b == 0) {
                // letter x: contributes inversions with previous y and z
                inv ^= (py ^ pz);
                px ^= 1;
            } else if (b == 1) {
                // letter y: contributes inversions with previous z
                inv ^= pz;
                py ^= 1;
            } else {
                // letter z: contributes no inversions
                pz ^= 1;
            }
            a = (a * MUL) % MOD; // advance to a_{n+1}
        }

        uint8_t mask = static_cast<uint8_t>((px << 2) | (py << 1) | pz);
        cnt[mask][inv]++;
    }

    // Count ordered pairs (i,j).
    uint64_t total = 0;
    for (int mask = 0; mask < 8; ++mask) {
        uint64_t c0 = cnt[mask][0];
        uint64_t c1 = cnt[mask][1];
        if (mask == 0) {
            // parity vector (0,0,0): need inv(A)=inv(B)
            total += c0 * c0 + c1 * c1;
        } else {
            // parity vector has two 1s: need inv(A) != inv(B)
            total += 2ULL * c0 * c1;
        }
    }
    return total;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Validation checkpoints from the statement.
    {
        uint64_t f10 = compute_F(10);
        uint64_t f100 = compute_F(100);
        if (f10 != 13ULL || f100 != 1224ULL) {
            cerr << "Validation failed: F(10)=" << f10
                 << ", F(100)=" << f100 << "\n";
            return 1;
        }
    }

    cout << compute_F(1'000'000) << "\n";
    return 0;
}
