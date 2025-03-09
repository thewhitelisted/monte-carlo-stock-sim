/*
monte carlo and heston model simulations for stocks and volatility
coded by @thewhitelisted
*/

// imports
use rand::prelude::*;
use rand_distr::{Normal, Distribution};
// use std::f64::consts::E;
use std::vec;
use plotters::prelude::*;
use tokio::task;

// heston params struct to hold data
#[derive(Clone, Copy)]
struct HestonParams {
    mu: f64,
    kappa: f64,
    theta: f64,
    xi: f64,
    rho: f64,
}

// async monte_carlo(s0, mu, sigma, dt, steps, simulations) returns 
// a vector of vectors of f64s representing the simulated stock paths
// async fn monte_carlo(
//     s0: f64,
//     mu: f64,
//     sigma: f64,
//     dt: f64,
//     steps: usize,
//     simulations: usize,
// ) -> Vec<Vec<f64>> {
//     let normal: Normal<f64> = Normal::new(0.0, 1.0).unwrap();

//     // task handling
//     let mut tasks: Vec<task::JoinHandle<Vec<f64>>> = Vec::new();

//     for _ in 0..simulations {
//         tasks.push(task::spawn(async move {
//             let mut rng: ThreadRng = rand::thread_rng();
//             let mut path: Vec<f64> = vec![s0; steps + 1];
            
//             for t in 1..=steps {
//                 let z: f64 = normal.sample(&mut rng);
//                 // cool epic stochastic differential equation
//                 path[t] = path[t - 1] * E.powf((mu - 0.5 * sigma.powi(2)) * dt + sigma * dt.sqrt() * z);
//             }
//             return path;
//         }));
//     }

//     // tie it up in a nice bow and return
//     let results = futures::future::join_all(tasks).await;
//     return results.into_iter().map(|r| r.unwrap()).collect();
// }

// heston_monte_carlo(s0, v0, params, dt, steps, simulations) returns
// a tuple of vectors of vectors of f64s representing the simulated stock and volatility paths
async fn heston_monte_carlo(
    s0: f64,
    v0: f64,
    params: HestonParams,
    dt: f64,
    steps: usize,
    simulations: usize,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    // task handling
    let mut tasks: Vec<task::JoinHandle<(Vec<f64>, Vec<f64>)>> = Vec::new();

    // go through each simulation
    for _ in 0..simulations {
        let params = params.clone();
        tasks.push(tokio::task::spawn(async move {
            let mut rng: ThreadRng = rand::thread_rng();
            let mut price_path: Vec<f64> = vec![s0; steps + 1];
            let mut vol_path: Vec<f64> = vec![v0; steps + 1];
            
            // Create the normal distribution outside the loop
            let normal: Normal<f64> = Normal::new(0.0, 1.0).unwrap();

            for t in 1..=steps {
                // brownian motion
                let z1: f64 = normal.sample(&mut rng);
                let z2: f64 = params.rho * z1 + 
                    (1.0 - params.rho.powi(2)).sqrt() * 
                    normal.sample(&mut rng);

                // volatility
                let prev_vol: f64 = vol_path[t-1].max(0.0);
                vol_path[t] = prev_vol + params.kappa * (params.theta - prev_vol) * dt
                    + params.xi * prev_vol.sqrt() * dt.sqrt() * z2;

                // log transforming wowie
                let drift: f64 = (params.mu - 0.5 * prev_vol) * dt;
                let diffusion: f64 = prev_vol.sqrt() * dt.sqrt() * z1;
                price_path[t] = price_path[t-1] * (drift + diffusion).exp();
            }

            return (price_path, vol_path);
        }));
    }
    let results: Vec<Result<(Vec<f64>, Vec<f64>), task::JoinError>> = futures::future::join_all(tasks).await;
    let (price, vol): (Vec<_>, Vec<_>) = results.into_iter().map(|r| {
        let res: (Vec<f64>, Vec<f64>) = r.unwrap();
        (res.0, res.1)
    }).unzip();

    (price, vol)
}

// plot_sims(paths, steps) returns a Result<(), Box<dyn std::error::Error>>
// plots the simulated stock paths and saves it as an image
// fn plot_sims(paths: &Vec<Vec<f64>>, steps: usize) -> Result<(), Box<dyn std::error::Error>> {
//     // create a bitmap backend
//     let root = BitMapBackend::new("monte_carlo.png", (800, 600)).into_drawing_area();
//     root.fill(&WHITE)?;

//     // create a chart
//     let mut chart = ChartBuilder::on(&root)
//         .caption("Monte Carlo Stock Simulation", ("Arial", 30).into_font())
//         .margin(10)
//         .x_label_area_size(40)
//         .y_label_area_size(50)
//         .build_cartesian_2d(0..steps, 0.0..paths.iter().flatten().cloned().fold(f64::MIN, f64::max))?;

//     chart.configure_mesh().draw()?;

//     // plotting stuff
//     for path in paths.iter().take(10) {
//         chart.draw_series(LineSeries::new(
//             (0..steps).map(|i| (i, path[i])),
//             &BLACK,
//         ))?;
//     }

//     // tie it up in a nice bow and return
//     root.present()?;
//     Ok(())
// }

// plot_heston(price_paths, vol_paths, steps) returns a Result<(), Box<dyn std::error::Error>>
// plots the simulated stock and volatility paths and saves it as an image
// i used gpt for this one cuz idk how to make things seeable
fn plot_heston(
    price_paths: &Vec<Vec<f64>>,
    vol_paths: &Vec<Vec<f64>>,
    steps: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    // two thingys
    let root = BitMapBackend::new("heston.png", (800, 1200)).into_drawing_area();
    root.fill(&WHITE)?;
    let (top, bottom) = root.split_vertically(600);

    // plotting top
    let max_price = price_paths.iter().flatten().cloned().fold(f64::MIN, f64::max);
    let mut chart = ChartBuilder::on(&top)
        .caption("Heston Price Paths", ("Arial", 20).into_font())
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(0..steps, 0.0..max_price * 1.1)?;

    chart.configure_mesh().draw()?;
    for path in price_paths.iter().take(10) {
        chart.draw_series(LineSeries::new(
            (0..=steps).map(|i| (i, path[i])),
            &BLUE,
        ))?;
    }

    // plotting bottom
    let max_vol = vol_paths.iter().flatten().cloned().fold(f64::MIN, f64::max);
    let mut chart = ChartBuilder::on(&bottom)
        .caption("Heston Volatility Paths", ("Arial", 20).into_font())
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(0..steps, 0.0..max_vol * 1.1)?;

    chart.configure_mesh().draw()?;
    for path in vol_paths.iter().take(10) {
        chart.draw_series(LineSeries::new(
            (0..=steps).map(|i| (i, path[i])),
            &RED,
        ))?;
    }

    root.present()?;
    Ok(())
}

#[tokio::main]
async fn main() {
    // parameters
    let params = HestonParams {
        mu: 0.05,
        kappa: 2.0,
        theta: 0.2,
        xi: 0.3,
        rho: -0.7,
    };
    let simulations: usize = 10_000;

    let (price_paths, vol_paths) = heston_monte_carlo(
        100.0,
        0.2,
        params,
        1.0 / 252.0,
        252,
        simulations,
    ).await;

    // let simulated_paths = monte_carlo(s0, mu, sigma, dt, steps, simulations).await;

    // Print first 5 steps of first simulation
    println!("Price Path 1: {:?}", &price_paths[0][..5]);
    println!("Vol Path 1: {:?}", &vol_paths[0][..5]);

    // Plot results
    if let Err(e) = plot_heston(&price_paths, &vol_paths, 252) {
        eprintln!("Plotting error: {}", e);
    }
}
