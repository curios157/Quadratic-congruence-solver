# Quadratic congruence equation solver

A solver for quadratic congruence equations `ax^2 + bx + c = d (mod n)`. Coefficients must be passed as command line arguments in the illustrated order and for other coefs than the modulus n allowed value range is `[-2^63+2, 2^63-2]` and for n `[2, 2^63-2]`. Notice that when setting coef a to 0, the solver finds solutions to a linear equation. Solutions, if exist, are given as members of the ring `Z/nZ`. In most cases, the solver should be able to find all solutions for a solvable equation but there can be few cases where only a subset of the complete set of solutions is returned.

## Use ##
Following examples use Rust's package manager Cargo to build and then run the binary file. If needed, one can separate these steps.

When running the following command, please set numerical values instead of placeholders a etc. Two dashes separate Cargo's and binary arguments.
```Rust
cargo run --release -- a b c d n
```
For example, to solve an equation `x^2 - x - 3 = 0 (mod 9)`, run the following
```Rust
cargo run --release -- 1 -1 -3 0 9
```
and the solver will return two solutions 4 and 6 that belong to the ring Z/9Z. Another, perhaps slightly more serious problem `-x^2 + x - 1 = 0 (mod 2147483647^2)` can be solved similarly
```Rust
cargo run --release -- -1 1 -1 0 $((2147483647 ** 2))
```
and it would return again two different solutions.

## Test ##

Run all, approx. one hundred, unit tests by
```rust
cargo test
```
