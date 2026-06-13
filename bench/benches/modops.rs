#[macro_use]
extern crate criterion;
use criterion::Criterion;
use num_modular::{
    FixedMersenne64, FixedMontgomery64, FixedProth64, FixedTrinomialSolinas64, ModularCoreOps,
    Montgomery, PreMulInv2by1, Reducer,
};
use rand::random;

// ── Benchmark helper macros ─────────────────────────────────────────────────

/// Register a benchmark for a single-arg operation (pow, sqr).
macro_rules! bench_fn1 {
    ($group:expr, $name:expr, $reducer:expr, $bases:expr, $op:ident $(, $arg:expr)?) => {
        $group.bench_function($name, |b| {
            b.iter(|| {
                $bases
                    .iter()
                    .map(|&v| $reducer.$op(v $(, $arg)?))
                    .reduce(|a, b| a.wrapping_add(b))
            })
        });
    };
}

/// Register a benchmark for `inv` (returns `Option`).
macro_rules! bench_fn_inv {
    ($group:expr, $name:expr, $reducer:expr, $bases:expr) => {
        $group.bench_function($name, |b| {
            b.iter(|| {
                $bases
                    .iter()
                    .map(|&v| $reducer.inv(v))
                    .reduce(|a, b| Some(a?.wrapping_add(b?)))
            })
        });
    };
}

/// Register a benchmark for two-arg operations (add, mul) using zip.
macro_rules! bench_fn2 {
    ($group:expr, $name:expr, $reducer:expr, $lhs:expr, $rhs:expr, $op:ident) => {
        $group.bench_function($name, |b| {
            b.iter(|| {
                $lhs.iter()
                    .zip($rhs.iter())
                    .map(|(a, b)| $reducer.$op(a, b))
                    .reduce(|a, b| a.wrapping_add(b))
            })
        });
    };
}

// ── u128 benchmarks ─────────────────────────────────────────────────────────

pub fn bench_u128(c: &mut Criterion) {
    const N: usize = 256;
    let mut cases: [(u128, u128, u128); N] = [(0, 0, 0); N];
    for i in 0..N {
        let (a, b, mut m) = (random(), random(), random());
        if m == 0 {
            m = 1;
        }
        cases[i] = (a, b, m);
    }

    let mut group = c.benchmark_group("u128 modular ops");
    group.bench_function("addm", |b| {
        b.iter(|| {
            cases
                .iter()
                .map(|&(a, b, m)| a.addm(b, &m))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("mulm", |b| {
        b.iter(|| {
            cases
                .iter()
                .map(|&(a, b, m)| a.mulm(b, &m))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });

    group.finish();
}

// ── mod_inv ─────────────────────────────────────────────────────────────────

pub fn bench_mod_inv_goldilocks(c: &mut Criterion) {
    const N: usize = 256;

    type Sol = FixedTrinomialSolinas64<64, 32, 1>;
    type Mer = FixedMersenne64<64, { (1u64 << 32) - 1 }>;
    const MOD: u64 = <Sol>::MODULUS;

    let mer = Mer::new(&MOD);
    let sol = Sol::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let fm = FixedMontgomery64::<MOD>::new(&MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut bases_mer = [0u64; N];
    let mut bases_sol = [0u64; N];
    let mut bases_fm = [0u64; N];
    for i in 0..N {
        let v = random::<u64>() % MOD;
        bases_mer[i] = mer.transform(if v == 0 { 1 } else { v });
        bases_sol[i] = sol.transform(if v == 0 { 1 } else { v });
        bases_fm[i] = fm.transform(if v == 0 { 1 } else { v });
    }

    let mut group = c.benchmark_group("mod_inv (2^64 - 2^32 + 1)");
    bench_fn_inv!(group, "Mersenne", mer, bases_mer);
    bench_fn_inv!(group, "Solinas", sol, bases_sol);
    bench_fn_inv!(group, "Montgomery", monty, bases_sol);
    bench_fn_inv!(group, "FixedMontgomery64", fm, bases_fm);
    bench_fn_inv!(group, "PreMulInv2by1", premul, bases_sol);
    group.finish();
}

pub fn bench_mod_inv_proth(c: &mut Criterion) {
    const N: usize = 256;

    type Pro = FixedProth64<32, 3>;
    const MOD: u64 = <Pro>::MODULUS;

    let pro = Pro::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let fm = FixedMontgomery64::<MOD>::new(&MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut bases = [0u64; N];
    let mut bases_fm = [0u64; N];
    for i in 0..N {
        let v = random::<u64>() % MOD;
        bases[i] = pro.transform(if v == 0 { 1 } else { v });
        bases_fm[i] = fm.transform(if v == 0 { 1 } else { v });
    }

    let mut group = c.benchmark_group("mod_inv proth 3*2^32+1");
    bench_fn_inv!(group, "FixedProth64", pro, bases);
    bench_fn_inv!(group, "Montgomery", monty, bases);
    bench_fn_inv!(group, "FixedMontgomery64", fm, bases_fm);
    bench_fn_inv!(group, "PreMulInv2by1", premul, bases);
    group.finish();
}

// ── mod_pow ─────────────────────────────────────────────────────────────────

pub fn bench_mod_pow_goldilocks(c: &mut Criterion) {
    const N: usize = 256;

    type Sol = FixedTrinomialSolinas64<64, 32, 1>;
    type Mer = FixedMersenne64<64, { (1u64 << 32) - 1 }>;
    const MOD: u64 = <Sol>::MODULUS;

    let mer = Mer::new(&MOD);
    let sol = Sol::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let fm = FixedMontgomery64::<MOD>::new(&MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut bases_mer = [0u64; N];
    let mut bases_sol = [0u64; N];
    let mut bases_fm = [0u64; N];
    for i in 0..N {
        let v = random::<u64>() % MOD;
        bases_mer[i] = mer.transform(v);
        bases_sol[i] = sol.transform(v);
        bases_fm[i] = fm.transform(v);
    }

    let exp = MOD - 2;

    let mut group = c.benchmark_group("mod_pow");
    bench_fn1!(
        group,
        "Mersenne (2^64 - 2^32 + 1)",
        mer,
        bases_mer,
        pow,
        &exp
    );
    bench_fn1!(
        group,
        "Solinas (2^64 - 2^32 + 1)",
        sol,
        bases_sol,
        pow,
        &exp
    );
    bench_fn1!(
        group,
        "Montgomery (2^64 - 2^32 + 1)",
        monty,
        bases_sol,
        pow,
        &exp
    );
    bench_fn1!(
        group,
        "FixedMontgomery64 (2^64 - 2^32 + 1)",
        fm,
        bases_fm,
        pow,
        &exp
    );
    bench_fn1!(
        group,
        "PreMulInv2by1 (2^64 - 2^32 + 1)",
        premul,
        bases_sol,
        pow,
        &exp
    );
    group.finish();
}

pub fn bench_mod_pow_proth(c: &mut Criterion) {
    const N: usize = 256;

    type Pro = FixedProth64<32, 3>;
    const MOD: u64 = <Pro>::MODULUS;

    let pro = Pro::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let fm = FixedMontgomery64::<MOD>::new(&MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut bases_pro = [0u64; N];
    let mut bases_fm = [0u64; N];
    for i in 0..N {
        let v = random::<u64>() % MOD;
        bases_pro[i] = pro.transform(v);
        bases_fm[i] = fm.transform(v);
    }

    let exp = MOD - 2;

    let mut group = c.benchmark_group("mod_pow proth 3*2^32+1");
    bench_fn1!(group, "FixedProth64", pro, bases_pro, pow, &exp);
    bench_fn1!(group, "Montgomery", monty, bases_pro, pow, &exp);
    bench_fn1!(group, "FixedMontgomery64", fm, bases_fm, pow, &exp);
    bench_fn1!(group, "PreMulInv2by1", premul, bases_pro, pow, &exp);
    group.finish();
}

// ── mod_add ─────────────────────────────────────────────────────────────────

pub fn bench_mod_add_goldilocks(c: &mut Criterion) {
    const N: usize = 256;

    type Sol = FixedTrinomialSolinas64<64, 32, 1>;
    type Mer = FixedMersenne64<64, { (1u64 << 32) - 1 }>;
    const MOD: u64 = <Sol>::MODULUS;

    let mer = Mer::new(&MOD);
    let sol = Sol::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let fm = FixedMontgomery64::<MOD>::new(&MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut lhs_mer = [0u64; N];
    let mut rhs_mer = [0u64; N];
    let mut lhs_sol = [0u64; N];
    let mut rhs_sol = [0u64; N];
    let mut lhs_fm = [0u64; N];
    let mut rhs_fm = [0u64; N];
    for i in 0..N {
        lhs_mer[i] = mer.transform(random::<u64>() % MOD);
        rhs_mer[i] = mer.transform(random::<u64>() % MOD);
        lhs_sol[i] = sol.transform(random::<u64>() % MOD);
        rhs_sol[i] = sol.transform(random::<u64>() % MOD);
        lhs_fm[i] = fm.transform(random::<u64>() % MOD);
        rhs_fm[i] = fm.transform(random::<u64>() % MOD);
    }

    let mut group = c.benchmark_group("mod_add (2^64 - 2^32 + 1)");
    bench_fn2!(group, "Mersenne", mer, lhs_mer, rhs_mer, add);
    bench_fn2!(group, "Solinas", sol, lhs_sol, rhs_sol, add);
    bench_fn2!(group, "Montgomery", monty, lhs_sol, rhs_sol, add);
    bench_fn2!(group, "FixedMontgomery64", fm, lhs_fm, rhs_fm, add);
    bench_fn2!(group, "PreMulInv2by1", premul, lhs_sol, rhs_sol, add);
    group.finish();
}

pub fn bench_mod_add_proth(c: &mut Criterion) {
    const N: usize = 256;

    type Pro = FixedProth64<32, 3>;
    const MOD: u64 = <Pro>::MODULUS;

    let pro = Pro::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let fm = FixedMontgomery64::<MOD>::new(&MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut lhs = [0u64; N];
    let mut rhs = [0u64; N];
    let mut lhs_fm = [0u64; N];
    let mut rhs_fm = [0u64; N];
    for i in 0..N {
        lhs[i] = pro.transform(random::<u64>() % MOD);
        rhs[i] = pro.transform(random::<u64>() % MOD);
        lhs_fm[i] = fm.transform(random::<u64>() % MOD);
        rhs_fm[i] = fm.transform(random::<u64>() % MOD);
    }

    let mut group = c.benchmark_group("mod_add proth 3*2^32+1");
    bench_fn2!(group, "FixedProth64", pro, lhs, rhs, add);
    bench_fn2!(group, "Montgomery", monty, lhs, rhs, add);
    bench_fn2!(group, "FixedMontgomery64", fm, lhs_fm, rhs_fm, add);
    bench_fn2!(group, "PreMulInv2by1", premul, lhs, rhs, add);
    group.finish();
}

// ── mod_mul ─────────────────────────────────────────────────────────────────

pub fn bench_mod_mul_goldilocks(c: &mut Criterion) {
    const N: usize = 256;

    type Sol = FixedTrinomialSolinas64<64, 32, 1>;
    type Mer = FixedMersenne64<64, { (1u64 << 32) - 1 }>;
    const MOD: u64 = <Sol>::MODULUS;

    let mer = Mer::new(&MOD);
    let sol = Sol::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let fm = FixedMontgomery64::<MOD>::new(&MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut lhs_mer = [0u64; N];
    let mut rhs_mer = [0u64; N];
    let mut lhs_sol = [0u64; N];
    let mut rhs_sol = [0u64; N];
    let mut lhs_fm = [0u64; N];
    let mut rhs_fm = [0u64; N];
    for i in 0..N {
        lhs_mer[i] = mer.transform(random::<u64>() % MOD);
        rhs_mer[i] = mer.transform(random::<u64>() % MOD);
        lhs_sol[i] = sol.transform(random::<u64>() % MOD);
        rhs_sol[i] = sol.transform(random::<u64>() % MOD);
        lhs_fm[i] = fm.transform(random::<u64>() % MOD);
        rhs_fm[i] = fm.transform(random::<u64>() % MOD);
    }

    let mut group = c.benchmark_group("mod_mul (2^64 - 2^32 + 1)");
    bench_fn2!(group, "Mersenne", mer, lhs_mer, rhs_mer, mul);
    bench_fn2!(group, "Solinas", sol, lhs_sol, rhs_sol, mul);
    bench_fn2!(group, "Montgomery", monty, lhs_sol, rhs_sol, mul);
    bench_fn2!(group, "FixedMontgomery64", fm, lhs_fm, rhs_fm, mul);
    bench_fn2!(group, "PreMulInv2by1", premul, lhs_sol, rhs_sol, mul);
    group.finish();
}

pub fn bench_mod_mul_proth(c: &mut Criterion) {
    const N: usize = 256;

    type Pro = FixedProth64<32, 3>;
    const MOD: u64 = <Pro>::MODULUS;

    let pro = Pro::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let fm = FixedMontgomery64::<MOD>::new(&MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut lhs = [0u64; N];
    let mut rhs = [0u64; N];
    let mut lhs_fm = [0u64; N];
    let mut rhs_fm = [0u64; N];
    for i in 0..N {
        lhs[i] = pro.transform(random::<u64>() % MOD);
        rhs[i] = pro.transform(random::<u64>() % MOD);
        lhs_fm[i] = fm.transform(random::<u64>() % MOD);
        rhs_fm[i] = fm.transform(random::<u64>() % MOD);
    }

    let mut group = c.benchmark_group("mod_mul proth 3*2^32+1");
    bench_fn2!(group, "FixedProth64", pro, lhs, rhs, mul);
    bench_fn2!(group, "Montgomery", monty, lhs, rhs, mul);
    bench_fn2!(group, "FixedMontgomery64", fm, lhs_fm, rhs_fm, mul);
    bench_fn2!(group, "PreMulInv2by1", premul, lhs, rhs, mul);
    group.finish();
}

criterion_group!(
    benches,
    bench_mod_inv_goldilocks,
    bench_mod_inv_proth,
    bench_u128,
    bench_mod_add_goldilocks,
    bench_mod_add_proth,
    bench_mod_mul_goldilocks,
    bench_mod_mul_proth,
    bench_mod_pow_goldilocks,
    bench_mod_pow_proth
);
criterion_main!(benches);
