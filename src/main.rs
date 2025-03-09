/*
cool monte carlo stock simulation. uses stochastic differential equation to simulate stock prices over time.
code by @thewhitelisted
*/

// imports
use rand::prelude::*;
use rand_distr::{Normal, Distribution};
use std::{f64::consts::E, vec};
use plotters::prelude::*;
use tokio::task;

// async monte_carlo(s0, mu, sigma, dt, steps, simulations) returns 
// a vector of vectors of f64s representing the simulated stock paths
async fn monte_carlo(
    s0: f64,
    mu: f64,
    sigma: f64,
    dt: f64,
    steps: usize,
    simulations: usize,
) -> Vec<Vec<f64>> {
    let normal: Normal<f64> = Normal::new(0.0, 1.0).unwrap();

    // task handling
    let mut tasks: Vec<task::JoinHandle<Vec<f64>>> = Vec::new();

    for _ in 0..simulations {
        tasks.push(task::spawn(async move {
            let mut rng: ThreadRng = rand::thread_rng();
            let mut path: Vec<f64> = vec![s0; steps + 1];
            
            for t in 1..=steps {
                let z: f64 = normal.sample(&mut rng);
                // cool epic stochastic differential equation
                path[t] = path[t - 1] * E.powf((mu - 0.5 * sigma.powi(2)) * dt + sigma * dt.sqrt() * z);
            }
            return path;
        }));
    }

    // tie it up in a nice bow and return
    let results = futures::future::join_all(tasks).await;
    return results.into_iter().map(|r| r.unwrap()).collect();
}

// plot_sims(paths, steps) returns a Result<(), Box<dyn std::error::Error>>
// plots the simulated stock paths and saves it as an image
fn plot_sims(paths: &Vec<Vec<f64>>, steps: usize) -> Result<(), Box<dyn std::error::Error>> {
    // create a bitmap backend
    let root = BitMapBackend::new("monte_carlo.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    // create a chart
    let mut chart = ChartBuilder::on(&root)
        .caption("Monte Carlo Stock Simulation", ("Arial", 30).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(0..steps, 0.0..paths.iter().flatten().cloned().fold(f64::MIN, f64::max))?;

    chart.configure_mesh().draw()?;

    // plotting stuff
    for path in paths.iter().take(10) {
        chart.draw_series(LineSeries::new(
            (0..steps).map(|i| (i, path[i])),
            &BLACK,
        ))?;
    }

    // tie it up in a nice bow and return
    root.present()?;
    Ok(())
}

#[tokio::main]
async fn main() {
    // parameters
    let s0: f64 = 100.0;
    let mu: f64 = 0.05;
    let sigma: f64 = 0.2;
    let dt: f64 = 1.0 / 252.0;
    let steps: usize = 252;
    let simulations: usize = 10_000;

    let simulated_paths = monte_carlo(s0, mu, sigma, dt, steps, simulations).await;

    // print out head of results
    for i in 0..simulations {
        println!("simulation {}: {:?}", i + 1, &simulated_paths[i][..5]);
    }

    // plots and saves results
    if let Err(e) = plot_sims(&simulated_paths, steps) {
        eprintln!("An error occurred: {}", e);
    }
}
