use std::cmp;
use std::convert::TryInto;
use std::convert::TryFrom;

const GCD_THRES: i64 = 4_294_967_295;


pub fn gcd(x: i64, y: i64) -> i64
{
    if cmp::max(x, y) < GCD_THRES {
        gcd_rec(
            x.try_into().unwrap(),
            y.try_into().unwrap()
        ).try_into().unwrap()
    } else {
        gcd_bin(
            x.try_into().unwrap(),
            y.try_into().unwrap()
        ).try_into().unwrap()
    }
}


fn gcd_bin(mut x: u64, mut y: u64) -> u64
{
    if x == 0 {return y;}
    if y == 0 {return x;}

    let shift: u32 = (x|y).trailing_zeros();
    x >>= x.trailing_zeros();

    let gcd: u64 = loop {
        y >>= y.trailing_zeros();
        if x > y {
            let t = y;
            y = x;
            x = t;
        }
        y -= x;
        if y == 0 {
            break x << shift;
        }
    };
    gcd
}


fn gcd_rec(x: u64, y: u64) -> u64
{
    if y == 0 {return x;}
    gcd_rec(y, x % y)
}


pub fn mod_exp_u64(mut b: u64, mut e: u64, m: u64) -> u64
{
    let mut r: u64 = 1;
    b %= m;

    if b == 0 {return 0;}

    while e > 0 {
        if e & 1 != 0 {
            r = mod_mult_u64(r, b, m);
        }
        e >>= 1;
        b = mod_mult_u64(b, b, m);
    }
    r
}


pub fn mod_mult_u64(mut x: u64, mut y: u64, m: u64) -> u64
{
    if y == 0 || x < m / y {
        return (x * y) % m;
    }
    let mut s = 0;

    while y > 0 {
        if y & 1 != 0 {
            s = (s + x) % m;
        }
        y >>= 1;
        x = (2 * x) % m;
    }
    s
}


pub fn mod_mult_i64(x: i64, y: i64, m: i64) -> i64
{
    // for non-negative but i64 type arguments
    mod_mult_u64(
        x.try_into().unwrap(),
        y.try_into().unwrap(),
        m.try_into().unwrap()
    ).try_into().unwrap()
}


pub fn mod_sum_i64(x: i64, y: i64, m: i64) -> i64
{
    // args x == (x % m), y == (y % m)
    if x <= m - y {
        (x + y) % m
    } else {
        (y - (m - x)) % m
    }
}


pub fn multip_inverse(x: i64, n: i64) -> i64
{
    let mut r_p = n;
    let mut r_n = x;
    let mut t_p = 0;
    let mut t_n = 1;

    while r_n > 0 {
        let q = r_p / r_n;

        let r_t = r_n;
        r_n = r_p - q * r_n;
        r_p = r_t;

        let t_t = t_n;
        t_n = t_p - q * t_n;
        t_p = t_t;
    }

    if r_p > 1 {return 0;}

    if t_p < 0 {
        t_p + n
    } else {
        t_p
    }
}


fn find_start_indices_for_diff_modulos(x: &Vec<(i64, i64)>) -> (Vec<u32>, Vec<u32>)
{
    let mut start_indices: Vec<u32> = Vec::new();
    let mut diff_mod_counts: Vec<u32> = Vec::new();

    let mut curr_mod = x[0].1;
    let mut curr_mod_count = 1;
    start_indices.push(0);

    for j in 1..x.len() {
        if x[j].1 != curr_mod {
            curr_mod = x[j].1;
            start_indices.push(j.try_into().unwrap());

            diff_mod_counts.push(curr_mod_count);
            curr_mod_count = 0;
            continue;
        }
        curr_mod_count += 1;
    }
    let result = (start_indices, diff_mod_counts);

    result
}


fn make_index_combinations(diff_mod_counts: Vec<u32>) -> Vec<Vec<u32>>
{
    let numbers = diff_mod_counts.len();

    let mut n_combi = 1;
    for j in 0..numbers {
        n_combi *= diff_mod_counts[j];
    }

    let mut combi_counter = 0;

    loop {
        break;
    }

    let mut x: Vec<Vec<u32>> = Vec::new();
    let temp = vec![0];
    x.push(temp);

    x
}


pub fn crt(x: Vec<(i64, i64)>, n: i64) -> Vec<i64>
{
    let mut sols: Vec<i64> = Vec::new();
    if x.len() == 0 {return sols;}

    let result = find_start_indices_for_diff_modulos(&x);
    let (start_indices, diff_mod_counts) = result;
    let combinations = make_index_combinations(diff_mod_counts);

    for combi in combinations {
        let mut s = 0;

        for (i, c) in combi.iter().enumerate() {
            let idx = usize::try_from(*c + start_indices[i]).unwrap();

            let z = n / x[idx].1;
            let inv = multip_inverse(z, x[idx].1);
            let r = mod_mult_i64(mod_mult_i64(x[idx].0, z, n), inv, n);

            s = mod_sum_i64(s, r, n);
        }
        sols.push(s);
    }

    sols
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gcd_small()
    {
        assert_eq!(gcd_rec(2,3), 1);
        assert_eq!(gcd_rec(3,2), 1);
        assert_eq!(gcd_rec(34,85), 17);
        assert_eq!(gcd_rec(224, 412), 4);
        assert_eq!(gcd_rec(526, 17_210), 2);
    }

    #[test]
    fn test_gcd_mid()
    {
        assert_eq!(gcd_rec(10_500, 975), 75);
        assert_eq!(gcd_rec(110_010, 750), 30);
        assert_eq!(gcd_rec(100_000, 15_888), 16);
        assert_eq!(gcd_rec(2147483647, 7483647), 1);
    }

    #[test]
    fn test_gcd_large()
    {
        assert_eq!(gcd_bin(1001116321, 1001118301), 1);
        assert_eq!(gcd_bin(9223372036854775807, 24141901), 7);
        assert_eq!(gcd_bin(9223372036854775807, 23523434234), 1);
        assert_eq!(gcd_bin(9223372036854775807, 9933434335423), 73);
        assert_eq!(gcd_bin(9223372036854775807, 3), 1);
    }

    #[test]
    fn test_mod_sum_i64()
    {
        assert_eq!(mod_sum_i64(9223372036854775781, 4, 9223372036854775783), 2);
        assert_eq!(mod_sum_i64(9223372036854775781, 9223372036854775781, 9223372036854775783), 9223372036854775779);
        assert_eq!(mod_sum_i64(9223372036854775782, 9223372036854775782, 9223372036854775783), 9223372036854775781);
        assert_eq!(mod_sum_i64(9223372036854775783, 9223372036854775783, 9223372036854775807), 9223372036854775759);
    }

    #[test]
    fn test_mod_exp_u64()
    {
        assert_eq!(mod_exp_u64(1, 1, 2), 1);
        assert_eq!(mod_exp_u64(3, 4, 3), 0);
        assert_eq!(mod_exp_u64(2, 17563959203, 35127918407), 29505221767);
        assert_eq!(mod_exp_u64(99876541124, 998899, 9223214), 7615604);
        assert_eq!(mod_exp_u64(2, 9999999, 9223372036854775807), 512);
        assert_eq!(mod_exp_u64(9987654, 999999901010111, 9223372036854775807), 2940910929841963431);
    }

}
