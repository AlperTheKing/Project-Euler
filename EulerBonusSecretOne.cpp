#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "lodepng.h"

using namespace std;

namespace {

struct Image {
    int w = 0;
    int h = 0;
    vector<uint8_t> gray;
};

string lower_ext(const string& path) {
    size_t pos = path.find_last_of('.');
    if (pos == string::npos) return "";
    string ext = path.substr(pos + 1);
    for (char& c : ext) c = static_cast<char>(tolower(static_cast<unsigned char>(c)));
    return ext;
}

bool read_token(istream& in, string& out) {
    string tok;
    while (in >> tok) {
        if (!tok.empty() && tok[0] == '#') {
            string line;
            getline(in, line);
            continue;
        }
        out = tok;
        return true;
    }
    return false;
}

bool load_pgm(const string& path, Image& img, string& err) {
    ifstream in(path, ios::binary);
    if (!in) {
        err = "failed to open PGM file";
        return false;
    }
    string magic;
    if (!read_token(in, magic)) {
        err = "missing PGM header";
        return false;
    }
    if (magic != "P5" && magic != "P2") {
        err = "unsupported PGM magic (expected P5 or P2)";
        return false;
    }
    string tok;
    int w = 0, h = 0, maxv = 0;
    auto read_int = [&](int& v) -> bool {
        if (!read_token(in, tok)) return false;
        try {
            v = stoi(tok);
        } catch (...) {
            return false;
        }
        return true;
    };
    if (!read_int(w) || !read_int(h) || !read_int(maxv)) {
        err = "invalid PGM header";
        return false;
    }
    if (w <= 0 || h <= 0 || maxv <= 0 || maxv > 255) {
        err = "invalid PGM dimensions or maxval";
        return false;
    }
    img.w = w;
    img.h = h;
    img.gray.assign(static_cast<size_t>(w) * h, 0);

    if (magic == "P5") {
        int c = in.get();
        if (c != EOF && !isspace(c)) {
            in.unget();
        }
        in.read(reinterpret_cast<char*>(img.gray.data()), img.gray.size());
        if (!in) {
            err = "unexpected EOF in PGM data";
            return false;
        }
        if (maxv != 255) {
            for (auto& v : img.gray) {
                v = static_cast<uint8_t>(v * 255 / maxv);
            }
        }
    } else {
        for (size_t i = 0; i < img.gray.size(); ++i) {
            if (!read_token(in, tok)) {
                err = "unexpected EOF in ASCII PGM data";
                return false;
            }
            int v = 0;
            try {
                v = stoi(tok);
            } catch (...) {
                err = "invalid ASCII PGM value";
                return false;
            }
            v = max(0, min(maxv, v));
            img.gray[i] = static_cast<uint8_t>(v * 255 / maxv);
        }
    }
    return true;
}

bool load_png(const string& path, Image& img, string& err) {
    vector<unsigned char> rgba;
    unsigned w = 0, h = 0;
    unsigned error = lodepng::decode(rgba, w, h, path);
    if (error) {
        err = lodepng_error_text(error);
        return false;
    }
    if (w == 0 || h == 0) {
        err = "empty PNG";
        return false;
    }
    img.w = static_cast<int>(w);
    img.h = static_cast<int>(h);
    img.gray.assign(static_cast<size_t>(w) * h, 0);
    bool non_gray = false;
    for (size_t i = 0; i < img.gray.size(); ++i) {
        unsigned char r = rgba[4 * i + 0];
        unsigned char g = rgba[4 * i + 1];
        unsigned char b = rgba[4 * i + 2];
        unsigned char a = rgba[4 * i + 3];
        if (a != 255) {
            r = static_cast<unsigned char>((r * a + 255 * (255 - a) + 127) / 255);
            g = static_cast<unsigned char>((g * a + 255 * (255 - a) + 127) / 255);
            b = static_cast<unsigned char>((b * a + 255 * (255 - a) + 127) / 255);
        }
        if (r == g && g == b) {
            img.gray[i] = r;
        } else {
            non_gray = true;
            int lum = (77 * r + 150 * g + 29 * b + 128) >> 8;
            img.gray[i] = static_cast<uint8_t>(lum);
        }
    }
    if (non_gray) {
        cerr << "Warning: PNG is not grayscale, using luma conversion.\n";
    }
    return true;
}

bool load_image(const string& path, Image& img, string& err) {
    string ext = lower_ext(path);
    if (ext == "png") {
        return load_png(path, img, err);
    }
    if (ext == "pgm") {
        return load_pgm(path, img, err);
    }
    string err_png;
    if (load_png(path, img, err_png)) {
        return true;
    }
    string err_pgm;
    if (load_pgm(path, img, err_pgm)) {
        return true;
    }
    err = "failed to decode image: " + err_png + "; " + err_pgm;
    return false;
}

inline uint8_t sum4_mod7(uint8_t a, uint8_t b, uint8_t c, uint8_t d) {
    int s = a + b + c + d;
    if (s >= 7) s -= 7;
    if (s >= 7) s -= 7;
    if (s >= 7) s -= 7;
    return static_cast<uint8_t>(s);
}

void apply_shift_sum(const vector<uint8_t>& src, vector<uint8_t>& dst,
                     int W, int H, int sx, int sy, int threads) {
    if (W <= 0 || H <= 0) return;
    vector<int> xp(W), xm(W), yp(H), ym(H);
    for (int x = 0; x < W; ++x) {
        int p = x + sx;
        if (p >= W) p -= W;
        int m = x - sx;
        if (m < 0) m += W;
        xp[x] = p;
        xm[x] = m;
    }
    for (int y = 0; y < H; ++y) {
        int p = y + sy;
        if (p >= H) p -= H;
        int m = y - sy;
        if (m < 0) m += H;
        yp[y] = p;
        ym[y] = m;
    }

    size_t total = static_cast<size_t>(W) * H;
    if (threads <= 1 || total < 100000) {
        for (int y = 0; y < H; ++y) {
            size_t row = static_cast<size_t>(y) * W;
            size_t rowp = static_cast<size_t>(yp[y]) * W;
            size_t rowm = static_cast<size_t>(ym[y]) * W;
            for (int x = 0; x < W; ++x) {
                dst[row + x] = sum4_mod7(src[rowp + x], src[rowm + x],
                                         src[row + xp[x]], src[row + xm[x]]);
            }
        }
        return;
    }

    threads = min(threads, H);
    vector<thread> pool;
    pool.reserve(threads);
    for (int t = 0; t < threads; ++t) {
        int y0 = (H * t) / threads;
        int y1 = (H * (t + 1)) / threads;
        pool.emplace_back([&, y0, y1]() {
            for (int y = y0; y < y1; ++y) {
                size_t row = static_cast<size_t>(y) * W;
                size_t rowp = static_cast<size_t>(yp[y]) * W;
                size_t rowm = static_cast<size_t>(ym[y]) * W;
                for (int x = 0; x < W; ++x) {
                    dst[row + x] = sum4_mod7(src[rowp + x], src[rowm + x],
                                             src[row + xp[x]], src[row + xm[x]]);
                }
            }
        });
    }
    for (auto& th : pool) th.join();
}

vector<int> base7_digits(uint64_t steps) {
    vector<int> digits;
    if (steps == 0) {
        digits.push_back(0);
        return digits;
    }
    while (steps > 0) {
        digits.push_back(static_cast<int>(steps % 7));
        steps /= 7;
    }
    return digits;
}

void apply_steps(vector<uint8_t>& img, int W, int H, uint64_t steps, int threads) {
    if (steps == 0) return;
    vector<uint8_t> tmp(img.size());
    vector<int> digits = base7_digits(steps);
    int sx = (W == 0) ? 0 : 1 % W;
    int sy = (H == 0) ? 0 : 1 % H;
    for (size_t k = 0; k < digits.size(); ++k) {
        int d = digits[k];
        for (int i = 0; i < d; ++i) {
            apply_shift_sum(img, tmp, W, H, sx, sy, threads);
            img.swap(tmp);
        }
        if (k + 1 < digits.size()) {
            if (W > 0) sx = static_cast<int>((static_cast<long long>(sx) * 7) % W);
            if (H > 0) sy = static_cast<int>((static_cast<long long>(sy) * 7) % H);
        }
    }
}

vector<uint8_t> simulate_naive(vector<uint8_t> img, int W, int H, int steps) {
    vector<uint8_t> tmp(img.size());
    for (int s = 0; s < steps; ++s) {
        apply_shift_sum(img, tmp, W, H, 1 % W, 1 % H, 1);
        img.swap(tmp);
    }
    return img;
}

uint8_t powmod4(uint64_t e) {
    int base = 4 % 7;
    int r = 1;
    while (e > 0) {
        if (e & 1ULL) r = (r * base) % 7;
        base = (base * base) % 7;
        e >>= 1ULL;
    }
    return static_cast<uint8_t>(r);
}

bool run_validation() {
    const int W = 9;
    const int H = 7;
    vector<uint8_t> img(static_cast<size_t>(W) * H);
    uint32_t seed = 1;
    for (size_t i = 0; i < img.size(); ++i) {
        seed = seed * 1664525u + 1013904223u;
        img[i] = static_cast<uint8_t>((seed >> 16) % 7);
    }

    const int steps = 20;
    vector<uint8_t> naive = simulate_naive(img, W, H, steps);
    vector<uint8_t> fast = img;
    apply_steps(fast, W, H, steps, 1);
    if (naive != fast) {
        cerr << "Validation failed: base-7 steps mismatch.\n";
        return false;
    }

    vector<uint8_t> naive7 = simulate_naive(img, W, H, 7);
    vector<uint8_t> tmp(img.size());
    apply_shift_sum(img, tmp, W, H, 7 % W, 7 % H, 1);
    if (naive7 != tmp) {
        cerr << "Validation failed: L^7 shift identity mismatch.\n";
        return false;
    }

    int sum0 = 0;
    for (uint8_t v : img) sum0 = (sum0 + v) % 7;
    apply_shift_sum(img, tmp, W, H, 1 % W, 1 % H, 1);
    int sum1 = 0;
    for (uint8_t v : tmp) sum1 = (sum1 + v) % 7;
    if (sum1 != (sum0 * 4) % 7) {
        cerr << "Validation failed: sum scaling mismatch.\n";
        return false;
    }

    vector<uint8_t> const_img(static_cast<size_t>(W) * H, 3);
    vector<uint8_t> const_fast = const_img;
    apply_steps(const_fast, W, H, 123, 1);
    uint8_t expected = static_cast<uint8_t>((3 * powmod4(123)) % 7);
    for (uint8_t v : const_fast) {
        if (v != expected) {
            cerr << "Validation failed: constant field check.\n";
            return false;
        }
    }

    cerr << "Validation checkpoints passed.\n";
    return true;
}

void write_pgm(const string& path, const vector<uint8_t>& img, int W, int H) {
    ofstream out(path, ios::binary);
    if (!out) {
        cerr << "Failed to write output file: " << path << "\n";
        return;
    }
    out << "P5\n" << W << " " << H << "\n255\n";
    vector<uint8_t> scaled(img.size());
    for (size_t i = 0; i < img.size(); ++i) {
        scaled[i] = static_cast<uint8_t>(img[i] * 255 / 6);
    }
    out.write(reinterpret_cast<const char*>(scaled.data()), scaled.size());
}

void write_color_png(const string& path, const vector<uint8_t>& img, int W, int H, int scale) {
    if (scale < 1) scale = 1;
    if (W <= 0 || H <= 0) return;
    unsigned W2 = static_cast<unsigned>(W * scale);
    unsigned H2 = static_cast<unsigned>(H * scale);
    vector<unsigned char> rgba(static_cast<size_t>(W2) * H2 * 4, 255);

    static const unsigned char palette[7][3] = {
        {255, 255, 255},  // 0
        {0, 0, 0},        // 1
        {220, 20, 60},    // 2
        {30, 144, 255},   // 3
        {34, 139, 34},    // 4
        {255, 165, 0},    // 5
        {148, 0, 211}     // 6
    };

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            uint8_t v = img[static_cast<size_t>(y) * W + x];
            const unsigned char* c = palette[v % 7];
            int y0 = y * scale;
            int x0 = x * scale;
            for (int dy = 0; dy < scale; ++dy) {
                size_t row = static_cast<size_t>(y0 + dy) * W2;
                for (int dx = 0; dx < scale; ++dx) {
                    size_t idx = (row + (x0 + dx)) * 4;
                    rgba[idx + 0] = c[0];
                    rgba[idx + 1] = c[1];
                    rgba[idx + 2] = c[2];
                    rgba[idx + 3] = 255;
                }
            }
        }
    }

    unsigned error = lodepng::encode(path, rgba, W2, H2);
    if (error) {
        cerr << "Failed to write PNG: " << lodepng_error_text(error) << "\n";
    }
}

void print_ascii(const vector<uint8_t>& img, int W, int H) {
    int hist[7] = {};
    for (uint8_t v : img) hist[v]++;
    int bg = 0;
    for (int i = 1; i < 7; ++i) {
        if (hist[i] > hist[bg]) bg = i;
    }

    int minx = W, miny = H, maxx = -1, maxy = -1;
    for (int y = 0; y < H; ++y) {
        size_t row = static_cast<size_t>(y) * W;
        for (int x = 0; x < W; ++x) {
            if (img[row + x] != bg) {
                minx = min(minx, x);
                miny = min(miny, y);
                maxx = max(maxx, x);
                maxy = max(maxy, y);
            }
        }
    }
    if (maxx < 0) {
        cerr << "ASCII preview skipped: all pixels match background value.\n";
        return;
    }

    int width = maxx - minx + 1;
    int height = maxy - miny + 1;
    int max_w = 160;
    int max_h = 60;
    int scale = 1;
    if (width > max_w) scale = max(scale, (width + max_w - 1) / max_w);
    if (height > max_h) scale = max(scale, (height + max_h - 1) / max_h);

    cerr << "ASCII preview (background=" << bg << ", scale=" << scale << "):\n";
    for (int y = miny; y <= maxy; y += scale) {
        string line;
        line.reserve(width / scale + 1);
        size_t row = static_cast<size_t>(y) * W;
        for (int x = minx; x <= maxx; x += scale) {
            line.push_back(img[row + x] == bg ? ' ' : '#');
        }
        while (!line.empty() && line.back() == ' ') line.pop_back();
        cout << line << "\n";
    }
}

void print_usage(const char* argv0) {
    cerr << "Usage: " << argv0 << " <image.(png|pgm)> [--threads N] [--steps N]\n"
         << "       [--out output.pgm] [--color-out output.png] [--scale N]\n"
         << "       [--no-ascii] [--no-validate]\n";
}

}  // namespace

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string input_path;
    string output_path = "output_mod7.pgm";
    string color_output_path;
    uint64_t steps = 1000000000000ULL;
    int threads = static_cast<int>(thread::hardware_concurrency());
    if (threads <= 0) threads = 1;
    bool ascii = true;
    bool validate = true;
    int scale = 6;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--threads" && i + 1 < argc) {
            try {
                threads = stoi(argv[++i]);
            } catch (...) {
                cerr << "Invalid --threads value.\n";
                return 1;
            }
        } else if (arg == "--steps" && i + 1 < argc) {
            try {
                steps = stoull(argv[++i]);
            } catch (...) {
                cerr << "Invalid --steps value.\n";
                return 1;
            }
        } else if (arg == "--out" && i + 1 < argc) {
            output_path = argv[++i];
        } else if (arg == "--color-out" && i + 1 < argc) {
            color_output_path = argv[++i];
        } else if (arg == "--scale" && i + 1 < argc) {
            try {
                scale = stoi(argv[++i]);
            } catch (...) {
                cerr << "Invalid --scale value.\n";
                return 1;
            }
        } else if (arg == "--no-ascii") {
            ascii = false;
        } else if (arg == "--no-validate") {
            validate = false;
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if (input_path.empty()) {
            input_path = arg;
        } else {
            cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (input_path.empty()) {
        print_usage(argv[0]);
        return 1;
    }

    if (validate && !run_validation()) {
        return 1;
    }

    Image img;
    string err;
    if (!load_image(input_path, img, err)) {
        cerr << "Failed to load image: " << err << "\n";
        return 1;
    }

    for (auto& v : img.gray) v %= 7;

    apply_steps(img.gray, img.w, img.h, steps, threads);

    write_pgm(output_path, img.gray, img.w, img.h);
    if (!color_output_path.empty()) {
        write_color_png(color_output_path, img.gray, img.w, img.h, scale);
    }
    if (ascii) {
        print_ascii(img.gray, img.w, img.h);
    } else {
        cerr << "ASCII preview disabled. Output written to " << output_path << "\n";
    }

    return 0;
}
