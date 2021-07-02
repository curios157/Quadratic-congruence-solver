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
        println!("Error when parsing command line arguments: {}", err);
        process::exit(1);
    });

    solver::solve_equation(equation);
}
