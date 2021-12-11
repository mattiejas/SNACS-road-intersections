mod anf;

fn main() {
    let now = std::time::Instant::now();

    let edges = split_line(
        &std::fs::read_to_string("../data/croatia/croatia-highways.csv").unwrap(),
        " ",
    );
    let mut model = anf::ANF::new(edges, 7, 128);
    let (neighbourhood_sizes, mean_sizes) = model.compute(5);
    println!("neighbourhood_sizes = {:#?}", neighbourhood_sizes);
    println!("mean = {:#?}", mean_sizes);

    println!(
        "ran for {}ms",
        (now.elapsed().as_nanos()) as f32 / 1000000f32
    );
}

fn split_line(line: &str, delim: &str) -> Vec<(usize, usize)> {
    line.split('\n')
        .filter(|s| !s.is_empty())
        .map(|s| {
            let (x, y) = s.trim().split_once(delim).unwrap();
            (
                x.trim().parse::<usize>().unwrap(),
                y.trim().parse::<usize>().unwrap(),
            )
        })
        .collect()
}
