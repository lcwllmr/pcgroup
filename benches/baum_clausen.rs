use std::hint::black_box;
use std::time::{Duration, Instant};

fn generate_input_from_param(k: u32) -> pcgroup::Presentation {
    // NOTE: 2-groups are usually worst case examples for the Baum-Clausen algorithm
    pcgroup::zoo::dihedral(1 << (k - 1))
}

fn my_algorithm(data: &mut pcgroup::Presentation) -> Vec<pcgroup::Representation> {
    pcgroup::irreducible_representations(data).unwrap()
}

fn run_for_param(k: u32, min_duration: Duration, min_iters: usize) -> (u128, Duration) {
    let n = (2_u128) * (1 << (k - 1) as u128);

    // Warm-up pass
    {
        let mut data = generate_input_from_param(k);
        assert_eq!(data.order(), n);
        black_box(my_algorithm(&mut data));
    }

    let mut total_iters = 0;
    let mut total_time = Duration::ZERO;
    let start_all = Instant::now();

    while total_time < min_duration || total_iters < min_iters {
        let mut data = generate_input_from_param(k);
        let start = Instant::now();
        black_box(my_algorithm(&mut data));
        total_time += start.elapsed();
        total_iters += 1;

        if start_all.elapsed() > Duration::from_secs(10) {
            break;
        }
    }

    (n, total_time / total_iters as u32)
}

fn main() {
    // Choose p values (geometric growth roughly doubling N each step)
    let k_values: Vec<u32> = (10..=18).collect();
    let mut results: Vec<(u32, u128, Duration)> = Vec::new(); // (p, N, Duration)

    println!("{:-^80}", " Arbitrary-N Scaling Analysis: O(N log N) ");
    println!(
        "{:>6} | {:>10} | {:>10} | {:>10} | {:>10} | {:>14}",
        "k", "N = 2*2^k", "Time", "Measured", "Expected", "c = T / (N lg N)"
    );
    println!("{:-<80}", "");

    for (i, &k) in k_values.iter().enumerate() {
        let (n, t) = run_for_param(k, Duration::from_millis(200), 10);
        results.push((k, n, t));

        let t_secs = t.as_secs_f64();
        let n_f = n as f64;
        let c_factor = (t_secs / (n_f * n_f.log2())) * 1e9; // ns per operation

        if i == 0 {
            println!(
                "{:>6} | {:>10} | {:>10.2?} | {:>10} | {:>10} | {:>11.2} ns",
                k, n, t, "-", "-", c_factor
            );
        } else {
            let (_, prev_n, prev_t) = results[i - 1];
            let prev_n_f = prev_n as f64;
            let prev_t_secs = prev_t.as_secs_f64();

            // Generalized theoretical ratio: f(N_i) / f(N_{i-1})
            let measured_ratio = t_secs / prev_t_secs;
            let expected_ratio = (n_f * n_f.log2()) / (prev_n_f * prev_n_f.log2());

            println!(
                "{:>6} | {:>10} | {:>10.2?} | {:>10.3} | {:>10.3} | {:>11.2} ns",
                k, n, t, measured_ratio, expected_ratio, c_factor
            );
        }
    }
    println!("{:-<80}", "");
}
