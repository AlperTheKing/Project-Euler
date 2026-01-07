use std::cmp::min;
use std::sync::Arc;
use std::thread;

const MOD: u64 = 1_000_000_007;

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if a == 0 {
        (b, 0, 1)
    } else {
        let (g, x1, y1) = extended_gcd(b % a, a);
        let x = y1 - (b / a) * x1;
        let y = x1;
        (g, x, y)
    }
}

fn mod_inverse(a: i64, m: i64) -> Option<i64> {
    let (g, x, _) = extended_gcd(a, m);
    if g != 1 {
        None
    } else {
        let mut res = x % m;
        if res < 0 {
            res += m;
        }
        Some(res)
    }
}

fn mod_pow(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    if modulus == 1 {
        return 0;
    }
    let mut result: u64 = 1;
    base %= modulus;
    while exp > 0 {
        if exp & 1 == 1 {
            result = ((result as u128 * base as u128) % modulus as u128) as u64;
        }
        base = ((base as u128 * base as u128) % modulus as u128) as u64;
        exp >>= 1;
    }
    result
}

struct PrimeFactorizer {
    spf: Vec<u32>,
}

impl PrimeFactorizer {
    fn new(max_n: usize) -> Self {
        let mut spf = vec![0u32; max_n + 1];
        for i in 0..=max_n {
            spf[i] = i as u32;
        }
        if max_n >= 1 {
            spf[1] = 1;
        }
        let limit = (max_n as f64).sqrt() as usize;
        for i in 2..=limit {
            if spf[i] == i as u32 {
                for j in (i * i..=max_n).step_by(i) {
                    if spf[j] == j as u32 {
                        spf[j] = i as u32;
                    }
                }
            }
        }
        Self { spf }
    }

    fn get_factors(&self, mut n: u32) -> Vec<u32> {
        let mut factors = Vec::new();
        while n > 1 {
            let p = self.spf[n as usize];
            factors.push(p);
            n /= p;
        }
        factors
    }

    fn euler_phi(&self, n: u32) -> u32 {
        let mut result = n;
        let mut temp = n;
        while temp > 1 {
            let p = self.spf[temp as usize];
            while temp % p == 0 {
                temp /= p;
            }
            result -= result / p;
        }
        result
    }
}

fn multiplicative_order(a: u64, m: u64, factorizer: &PrimeFactorizer) -> Option<u64> {
    if gcd(a, m) != 1 {
        return None;
    }
    if m == 1 {
        return Some(1);
    }
    let phi = factorizer.euler_phi(m as u32) as u64;
    let mut factors = factorizer.get_factors(phi as u32);
    factors.sort_unstable();
    factors.dedup();

    let mut order = phi;
    for &p in &factors {
        let p_u64 = p as u64;
        while order % p_u64 == 0 && mod_pow(a, order / p_u64, m) == 1 {
            order /= p_u64;
        }
    }
    Some(order)
}

fn min_l_magnitude(p: u64, q: u64) -> u64 {
    let limit = (q + p - 1) / p;
    let mut pow10 = 1u64;
    let mut l = 1u64;
    while pow10 < limit {
        pow10 *= 10;
        l += 1;
    }
    l
}

fn solve_n(u: u64, v: u64, factorizer: &PrimeFactorizer) -> u64 {
    let p = u * u * u;
    let q = v * v * v;
    if p >= 10 * q {
        return 0;
    }

    let k_val = 10 * q - p;
    let min_l = min_l_magnitude(p, q);

    let mut best_l = u64::MAX;
    let mut best_d = 10u64;

    for d in 1..=9u64 {
        if p * (d + 1) >= 10 * q {
            continue;
        }

        let g = gcd(k_val, d);
        let m_prime = k_val / g;
        if gcd(10, m_prime) != 1 {
            continue;
        }

        let l_base = match multiplicative_order(10, m_prime, factorizer) {
            Some(ord) => ord,
            None => continue,
        };

        let l_cand = ((min_l + l_base - 1) / l_base) * l_base;
        if l_cand < best_l || (l_cand == best_l && d < best_d) {
            best_l = l_cand;
            best_d = d;
        }
    }

    if best_l == u64::MAX {
        return 0;
    }

    let ten_pow_l = mod_pow(10, best_l, MOD);
    let term1 = ((best_d % MOD) as u128 * (q % MOD) as u128 % MOD as u128) as u64;
    let term2 = (ten_pow_l + MOD - 1) % MOD;
    let k_inv = mod_inverse(k_val as i64, MOD as i64).unwrap() as u64;

    let mut n_mod = (term1 as u128 * term2 as u128 % MOD as u128) as u64;
    n_mod = (n_mod as u128 * k_inv as u128 % MOD as u128) as u64;
    n_mod
}

fn compute_total_serial(m: u64, factorizer: &PrimeFactorizer) -> u64 {
    let mut total = 0u64;
    for u in 1..=m {
        for v in 1..=m {
            if gcd(u, v) == 1 {
                total += solve_n(u, v, factorizer);
                if total >= MOD {
                    total %= MOD;
                }
            }
        }
    }
    total % MOD
}

fn compute_total(m: u64, factorizer: &PrimeFactorizer, threads: usize) -> u64 {
    if threads <= 1 {
        return compute_total_serial(m, factorizer);
    }

    let m_usize = m as usize;
    let threads = min(threads, m_usize.max(1));
    let chunk = (m_usize + threads - 1) / threads;
    let factorizer = Arc::new(factorizer);

    let mut handles = Vec::new();
    for t in 0..threads {
        let start = t * chunk + 1;
        if start > m_usize {
            break;
        }
        let end = min(m_usize, (t + 1) * chunk);
        let factorizer = Arc::clone(&factorizer);
        handles.push(thread::spawn(move || {
            let mut local = 0u64;
            for u in start..=end {
                let u64u = u as u64;
                for v in 1..=m {
                    if gcd(u64u, v) == 1 {
                        local += solve_n(u64u, v, factorizer.as_ref());
                        if local >= MOD {
                            local %= MOD;
                        }
                    }
                }
            }
            local % MOD
        }));
    }

    let mut total = 0u64;
    for handle in handles {
        total = (total + handle.join().unwrap()) % MOD;
    }
    total
}

fn run_validations(factorizer: &PrimeFactorizer) -> bool {
    let mut ok = true;

    let t1 = compute_total_serial(1, factorizer);
    if t1 != 1 {
        eprintln!("Validation failed: T(1) expected 1, got {}", t1);
        ok = false;
    } else {
        eprintln!("Validation checkpoint: T(1) ok.");
    }

    let expected_t3 = 262429173u64;
    let t3 = compute_total_serial(3, factorizer);
    if t3 != expected_t3 {
        eprintln!(
            "Validation failed: T(3) expected {}, got {}",
            expected_t3, t3
        );
        ok = false;
    } else {
        eprintln!("Validation checkpoint: T(3) ok.");
    }

    if ok {
        eprintln!("All validation checkpoints passed.");
    }
    ok
}

fn main() {
    let m_limit: u64 = 200;
    let max_n = 80_000_000usize; // 10 * 200^3

    eprintln!("Building SPF up to {}...", max_n);
    let factorizer = PrimeFactorizer::new(max_n);

    eprintln!("Running validation checkpoints...");
    if !run_validations(&factorizer) {
        return;
    }

    let threads = thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .min(16);

    let total = compute_total(m_limit, &factorizer, threads);
    println!("T({}) modulo {}: {}", m_limit, MOD, total);
}
