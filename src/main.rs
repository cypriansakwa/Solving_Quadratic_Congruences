fn main() {
    let b = 297;
    let m = 313;

    println!("Determine solutions for x^2 \u{2261} {} (mod {})", b, m);

    match solve_quadratic_congruence(b, m) {
        Some(solutions) => {
            if solutions.is_empty() {
                println!("x^2 \u{2261} {} (mod {}) has no solutions", b, m);
            } else {
                println!("Solutions found: {:?}", solutions);
            }
        }
        None => {
            println!("x^2 \u{2261} {} (mod {}) has no solutions", b, m);
        }
    }
}

fn solve_quadratic_congruence(b: i64, m: i64) -> Option<Vec<i64>> {
    if m % 2 == 0 {
        return solve_mod_power_of_two(b, m);
    }

    if legendre_symbol(b, m) != 1 {
        return None; // b is not a quadratic residue mod m then no solutions
    }

    // Tonelli-Shanks algorithm for prime moduli
    let mut q = m - 1;
    let mut s = 0;
    while q % 2 == 0 {
        q /= 2;
        s += 1;
    }

    let mut z = 2;
    while legendre_symbol(z, m) != -1 {
        z += 1;
    }

    let mut c = mod_exp(z, q, m);
    let mut r = mod_exp(b, (q + 1) / 2, m);
    let mut t = mod_exp(b, q, m);
    let mut m_s = s;

    while t != 1 {
        let mut t2i = t;
        let mut i = 0;
        for j in 1..m_s {
            t2i = (t2i * t2i) % m;
            if t2i == 1 {
                i = j;
                break;
            }
        }

        let b2i = mod_exp(c, 2_i64.pow((m_s - i - 1) as u32), m);
        r = (r * b2i) % m;
        c = (b2i * b2i) % m;
        t = (t * c) % m;
        m_s = i;
    }

    Some(vec![r, m - r])
}

fn solve_mod_power_of_two(b: i64, m: i64) -> Option<Vec<i64>> {
    if m == 2 {
        return Some(vec![b % 2]);
    }

    let k = (m as f64).log2() as i64;

    // Solution to x^2 \u{2261} 1 (mod 2^k)
    if b % 2 == 1 {
        let mut solutions = vec![1, m - 1];
        for i in 1..k {
            let step = 1 << i;
            let new_solutions: Vec<_> = solutions.iter().map(|&x| (x + step) % m).collect();
            for sol in new_solutions {
                if !solutions.contains(&sol) {
                    solutions.push(sol);
                }
            }
        }
        return Some(solutions);
    }

    // For b even
    None
}

fn legendre_symbol(a: i64, p: i64) -> i64 {
    let ls = mod_exp(a, (p - 1) / 2, p);
    if ls == p - 1 {
        -1
    } else {
        ls
    }
}

fn mod_exp(base: i64, exp: i64, modulus: i64) -> i64 {
    let mut result = 1;
    let mut base = base % modulus;
    let mut exp = exp;

    while exp > 0 {
        if exp % 2 == 1 {
            result = (result * base) % modulus;
        }
        exp >>= 1;
        base = (base * base) % modulus;
    }
    result
}


