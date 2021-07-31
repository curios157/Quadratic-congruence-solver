# Quadratic congruence equation solver

[![main](https://github.com/elmomoilanen/Quadratic-congruence-solver/actions/workflows/main.yml/badge.svg)](https://github.com/elmomoilanen/Quadratic-congruence-solver/actions/workflows/main.yml)

A solver for quadratic congruence equations `ax^2 + bx + c = d (mod n)`. Coefficients of the equation must be passed as command line arguments in the illustrated order. Allowed value range for modulus `n` is `[2, 2^63-2]` and for the others `[-2^63+2, 2^63-2]`. Notice that when setting the coefficient `a` to zero, the solver finds solutions to a linear equation instead of quadratic. Other interesting subcase, namely an equation to determine whether `d` is a quadratic residue, occurs if `a` is set equal to one and coefs `b` and `c` to zero. Solutions, if exist, are given as members of ring `Z/nZ`. In most cases, the solver should be able to find all solutions for a solvable equation but there can be few cases where only a subset of the complete set of solutions is returned.

## Use ##
Following examples use Rust's package manager Cargo to build and then run the binary file. Building takes place only when running for the first time or if the code has been modified since the last usage. Optionally, one may separately run the build command with Cargo or rustc and subsequently execute directly the binary without Cargo.

When running the following command in a shell, please use numerical values instead of character placeholders a and b etc. Notice that two dashes separate Cargo and binary's arguments.

```Bash
cargo run --release -- a b c d n
```

To shortly mention the other usage option, after the first release build (done e.g. with Cargo) the binary file can be run directly in the following manner

```bash
./target/release/quadratic_congruence_solver a b c d n
```

To solve actual quadratic or linear congruence equations, consider the following examples.

For example, in order to solve an equation of the form `3x^2 + 7x - 9 = 1 (mod 1729)`, run the following command
```Bash
cargo run --release -- 3 7 -9 1 1729
```
and as a result the solver will return all four solutions 1, 573, 820 and 1483 that belong to ring Z/1729Z.

Another, perhaps slightly more serious problem `-x^2 + x - 1 = 0 (mod 2147483647^2)` can be solved similarly
```Bash
cargo run --release -- -1 1 -1 0 $((2147483647 ** 2))
```
and it would return two different solutions.

## Test ##

Run all, approx. one hundred, unit tests by
```Bash
cargo test
```
