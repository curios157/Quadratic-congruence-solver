//! Solver for quadratic congruence equations `ax^2 + bx + c = d (mod n)`.
//!
//! Coefficients of the equation must be passed as command line arguments.
//! Allowed value range for modulus n is `[2, 2^63-2]` and for the others
//! `[-2^63+2, 2^63-2]`. Invalid coef values make the program to prompt
//! a specific error message to the shell.
//!
//! The congruence equation will be linear if `a = 0`. Other interesting
//! subcase occurs if `a = 1` and coefs `b` and `c` are set to zero,
//! meaning that the equation to be solved is of the form `x^2 = d (mod n)`.
//!
//! Solutions, if exist, are given as members of ring `Z/nZ`. In most
//! cases the solver should be able to find all solutions for a solvable
//! equation but there can be cases where only a subset of the complete
//! set of solutions is returned.
//!
//! There are two recommended ways to call this program. First of these
//! uses Rust's package manager Cargo to do both building and running
//! of the binary. Second way is to build with Cargo but run the binary
//! directly without Cargo.
//!
//! E.g., the first case happens as follows
//! ```bash
//! cargo run --release -- a b c d n
//! ```
//!
//! and the second as follows
//! ```bash
//! cargo build --release
//! ./target/release/quadratic_congruence_solver a b c d n
//! ```
//!
//! It's up to user which way to use.
//!
use std::{env, process};

mod equation;
mod euclid;
mod prime;
mod residue;
mod solver;

pub use equation::{Coefs, Equation};

fn main() {
    let args: Vec<String> = env::args().collect();

    let equation = Equation::new(&args).unwrap_or_else(|err| {
        if err == "help" {
            process::exit(0);
        }
        eprintln!("Error when parsing command line arguments: {}", err);
        process::exit(1);
    });

    solver::solve_equation(equation);
}
