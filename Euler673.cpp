#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;

constexpr i64 MOD = 999'999'937LL;

struct Comp {
    std::vector<int> b;
    std::vector<int> d;
};

i64 mod_pow(i64 a, int e) {
    i64 r = 1;
    i64 x = a % MOD;
    int p = e;
    while (p > 0) {
        if ((p & 1) != 0) {
            r = static_cast<i64>((static_cast<__int128>(r) * x) % MOD);
        }
        p >>= 1;
        if (p > 0) {
            x = static_cast<i64>((static_cast<__int128>(x) * x) % MOD);
        }
    }
    return r;
}

std::vector<std::pair<int, int>> parse_pairs(const std::string& data) {
    std::vector<std::pair<int, int>> out;
    std::istringstream in(data);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        const std::size_t pos = line.find(',');
        assert(pos != std::string::npos);
        const int a = std::stoi(line.substr(0, pos));
        const int b = std::stoi(line.substr(pos + 1));
        out.push_back({a, b});
    }
    return out;
}

std::vector<Comp> build_components(int n, const std::vector<std::pair<int, int>>& beds,
                                   const std::vector<std::pair<int, int>>& desks) {
    std::vector<int> b(static_cast<std::size_t>(n + 1));
    std::vector<int> d(static_cast<std::size_t>(n + 1));
    for (int i = 1; i <= n; ++i) {
        b[static_cast<std::size_t>(i)] = i;
        d[static_cast<std::size_t>(i)] = i;
    }

    for (const auto& p : beds) {
        b[static_cast<std::size_t>(p.first)] = p.second;
        b[static_cast<std::size_t>(p.second)] = p.first;
    }
    for (const auto& p : desks) {
        d[static_cast<std::size_t>(p.first)] = p.second;
        d[static_cast<std::size_t>(p.second)] = p.first;
    }

    std::vector<char> seen(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> loc(static_cast<std::size_t>(n + 1), -1);
    std::vector<Comp> comps;

    for (int s = 1; s <= n; ++s) {
        if (seen[static_cast<std::size_t>(s)] != 0) {
            continue;
        }

        std::vector<int> stack{s};
        seen[static_cast<std::size_t>(s)] = 1;
        std::vector<int> verts;

        while (!stack.empty()) {
            const int v = stack.back();
            stack.pop_back();
            verts.push_back(v);

            const int nb = b[static_cast<std::size_t>(v)];
            const int nd = d[static_cast<std::size_t>(v)];
            if (seen[static_cast<std::size_t>(nb)] == 0) {
                seen[static_cast<std::size_t>(nb)] = 1;
                stack.push_back(nb);
            }
            if (seen[static_cast<std::size_t>(nd)] == 0) {
                seen[static_cast<std::size_t>(nd)] = 1;
                stack.push_back(nd);
            }
        }

        for (int i = 0; i < static_cast<int>(verts.size()); ++i) {
            loc[static_cast<std::size_t>(verts[static_cast<std::size_t>(i)])] = i;
        }

        Comp c;
        c.b.resize(verts.size());
        c.d.resize(verts.size());
        for (int i = 0; i < static_cast<int>(verts.size()); ++i) {
            const int v = verts[static_cast<std::size_t>(i)];
            c.b[static_cast<std::size_t>(i)] =
                loc[static_cast<std::size_t>(b[static_cast<std::size_t>(v)])];
            c.d[static_cast<std::size_t>(i)] =
                loc[static_cast<std::size_t>(d[static_cast<std::size_t>(v)])];
        }

        for (int v : verts) {
            loc[static_cast<std::size_t>(v)] = -1;
        }

        comps.push_back(std::move(c));
    }

    return comps;
}

bool try_map_root(const Comp& a, const Comp& b, int start_b) {
    const int m = static_cast<int>(a.b.size());
    std::vector<int> map_a(static_cast<std::size_t>(m), -1);
    std::vector<int> map_b(static_cast<std::size_t>(m), -1);
    std::vector<int> stack;
    stack.reserve(static_cast<std::size_t>(m));

    map_a[0] = start_b;
    map_b[static_cast<std::size_t>(start_b)] = 0;
    stack.push_back(0);
    int assigned = 1;

    auto step = [&](int nx, int ny) -> bool {
        const int cur = map_a[static_cast<std::size_t>(nx)];
        if (cur == -1) {
            if (map_b[static_cast<std::size_t>(ny)] != -1) {
                return false;
            }
            map_a[static_cast<std::size_t>(nx)] = ny;
            map_b[static_cast<std::size_t>(ny)] = nx;
            stack.push_back(nx);
            ++assigned;
            return true;
        }
        return cur == ny;
    };

    while (!stack.empty()) {
        const int x = stack.back();
        stack.pop_back();
        const int y = map_a[static_cast<std::size_t>(x)];

        const int nx_b = a.b[static_cast<std::size_t>(x)];
        const int ny_b = b.b[static_cast<std::size_t>(y)];
        if (!step(nx_b, ny_b)) {
            return false;
        }

        const int nx_d = a.d[static_cast<std::size_t>(x)];
        const int ny_d = b.d[static_cast<std::size_t>(y)];
        if (!step(nx_d, ny_d)) {
            return false;
        }
    }

    return assigned == m;
}

int isomorphism_count(const Comp& a, const Comp& b) {
    if (a.b.size() != b.b.size()) {
        return 0;
    }

    const int m = static_cast<int>(b.b.size());
    int count = 0;
    for (int s = 0; s < m; ++s) {
        if ((a.b[0] == 0) != (b.b[static_cast<std::size_t>(s)] == s)) {
            continue;
        }
        if ((a.d[0] == 0) != (b.d[static_cast<std::size_t>(s)] == s)) {
            continue;
        }
        if (try_map_root(a, b, s)) {
            ++count;
        }
    }
    return count;
}

i64 count_valid_permutations(int n, const std::vector<std::pair<int, int>>& beds,
                             const std::vector<std::pair<int, int>>& desks) {
    const std::vector<Comp> comps = build_components(n, beds, desks);
    const int c = static_cast<int>(comps.size());

    std::vector<char> used(static_cast<std::size_t>(c), 0);
    i64 ans = 1;

    for (int i = 0; i < c; ++i) {
        if (used[static_cast<std::size_t>(i)] != 0) {
            continue;
        }

        used[static_cast<std::size_t>(i)] = 1;
        const int aut = isomorphism_count(comps[static_cast<std::size_t>(i)],
                                          comps[static_cast<std::size_t>(i)]);
        assert(aut > 0);

        int mult = 1;
        for (int j = i + 1; j < c; ++j) {
            if (used[static_cast<std::size_t>(j)] != 0) {
                continue;
            }
            if (isomorphism_count(comps[static_cast<std::size_t>(i)],
                                  comps[static_cast<std::size_t>(j)]) > 0) {
                used[static_cast<std::size_t>(j)] = 1;
                ++mult;
            }
        }

        i64 fac = 1;
        for (int t = 2; t <= mult; ++t) {
            fac = static_cast<i64>((static_cast<__int128>(fac) * t) % MOD);
        }

        ans = static_cast<i64>((static_cast<__int128>(ans) * mod_pow(aut, mult)) % MOD);
        ans = static_cast<i64>((static_cast<__int128>(ans) * fac) % MOD);
    }

    return ans;
}

const std::string kBeds500 = R"BED(1,352
2,251
3,243
4,311
5,312
6,231
7,416
8,166
9,384
10,330
11,89
12,316
13,394
14,57
15,177
16,341
17,149
18,199
19,310
20,421
21,459
22,189
23,308
25,46
26,460
27,63
28,34
29,295
30,412
31,49
32,247
33,326
35,429
36,211
37,104
38,376
39,263
40,187
41,500
42,161
43,407
44,60
45,118
47,85
48,56
50,403
51,239
52,415
53,444
54,246
55,288
58,410
59,432
61,160
62,438
65,225
67,392
68,333
69,478
70,365
71,213
72,294
73,360
74,499
75,267
76,493
77,108
78,428
79,440
80,116
81,117
82,257
83,306
84,195
86,216
87,226
88,151
90,172
91,355
93,237
94,156
95,147
96,343
97,424
98,112
100,240
101,163
102,371
103,361
106,315
107,230
109,227
110,336
111,377
113,253
114,367
115,139
119,290
120,121
122,174
123,159
124,282
125,142
126,435
127,452
128,338
129,445
130,448
131,183
132,179
133,196
134,268
135,170
136,357
137,368
138,242
141,400
143,278
144,208
145,204
146,378
148,252
152,238
153,380
154,269
155,270
157,373
158,339
162,272
164,332
167,466
168,327
169,266
173,471
175,494
176,463
178,433
180,198
181,258
182,464
184,467
186,409
188,191
190,346
192,446
193,418
197,474
200,354
201,479
202,436
205,404
206,475
207,349
210,476
212,372
214,379
217,232
218,320
219,366
220,356
221,299
222,426
223,484
224,264
228,359
229,405
233,483
234,391
235,381
236,262
241,420
244,393
245,337
248,296
249,370
250,325
254,351
255,302
256,364
260,431
261,447
265,318
271,397
273,375
274,469
275,453
276,328
279,313
280,292
283,434
284,495
285,386
286,454
287,385
289,363
291,304
293,382
297,389
298,383
301,461
303,347
305,369
309,439
314,344
317,491
319,443
321,470
322,437
323,411
324,406
331,340
334,358
335,456
342,414
345,396
350,390
374,408
387,422
395,487
398,423
399,485
401,441
402,473
413,457
417,449
425,451
427,477
430,496
442,488
450,480
462,490
465,482
472,492
481,486
497,498
)BED";

const std::string kDesks500 = R"DESK(1,106
2,429
4,94
6,349
7,362
8,100
9,368
10,477
11,132
12,328
13,319
14,185
15,225
16,279
17,72
18,51
19,482
20,270
21,50
22,449
23,398
24,492
25,478
26,360
27,428
28,240
29,103
30,32
31,484
33,296
34,371
35,251
36,444
37,460
38,491
39,180
40,264
41,81
42,67
43,472
44,173
45,209
46,135
47,218
48,69
49,137
52,324
54,210
55,330
56,316
57,436
58,212
59,193
61,233
62,211
63,158
64,439
65,419
66,455
70,153
71,198
73,283
74,76
75,250
77,418
78,191
79,144
80,352
82,195
85,192
86,329
87,340
88,306
89,181
90,299
91,383
92,358
93,276
95,119
96,176
97,129
98,235
99,487
101,336
102,332
104,183
105,281
107,432
109,273
110,416
111,304
112,465
113,385
114,167
115,497
116,399
117,404
118,384
121,450
122,247
124,486
125,402
127,344
130,175
131,310
133,149
134,147
136,289
138,453
139,184
140,351
141,256
142,422
143,440
145,481
146,201
148,293
150,231
151,170
152,223
154,314
155,395
156,249
157,414
159,341
160,446
161,410
162,286
163,166
164,307
165,242
168,282
172,234
174,224
177,320
178,258
179,433
182,461
186,257
187,312
188,339
189,388
190,373
194,325
196,322
197,427
200,483
202,357
204,268
205,206
207,285
208,244
213,425
214,294
215,354
216,420
217,382
219,297
220,378
222,400
226,323
229,337
230,499
237,437
238,408
239,489
241,405
243,389
245,392
246,361
248,364
252,424
253,454
254,374
255,301
259,292
260,409
261,347
262,280
263,451
265,381
266,379
267,355
271,456
272,369
274,298
275,431
277,412
278,473
284,490
287,494
288,413
290,343
291,396
295,348
300,315
302,464
303,345
305,468
308,480
309,448
311,331
313,463
317,495
318,367
321,438
327,366
334,426
338,466
342,346
350,445
353,471
356,447
363,476
365,415
370,442
372,458
375,434
376,462
377,423
380,406
386,469
387,459
390,498
391,479
394,496
397,485
401,457
403,474
407,470
411,452
417,443
421,430
435,475
467,493
488,500
)DESK";

}  // namespace

int main() {
    assert(count_valid_permutations(4, {{2, 3}}, {{1, 3}, {2, 4}}) == 2LL);

    assert(count_valid_permutations(6, {{1, 2}, {3, 4}, {5, 6}},
                                    {{3, 6}, {4, 5}}) == 8LL);

    const std::vector<std::pair<int, int>> beds36{
        {2, 13}, {4, 30}, {5, 27}, {6, 16}, {10, 18}, {12, 35}, {14, 19},
        {15, 20}, {17, 26}, {21, 32}, {22, 33}, {24, 34}, {25, 28},
    };
    const std::vector<std::pair<int, int>> desks36{
        {1, 35}, {2, 22}, {3, 36}, {4, 28}, {5, 25}, {7, 18}, {9, 23},
        {13, 19}, {14, 33}, {15, 34}, {20, 24}, {26, 29}, {27, 30},
    };
    assert(count_valid_permutations(36, beds36, desks36) == 663'552LL);

    const std::vector<std::pair<int, int>> beds500 = parse_pairs(kBeds500);
    const std::vector<std::pair<int, int>> desks500 = parse_pairs(kDesks500);
    assert(beds500.size() == 235);
    assert(desks500.size() == 236);

    std::cout << count_valid_permutations(500, beds500, desks500) << "\n";
    return 0;
}
