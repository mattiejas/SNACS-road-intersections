mod anf;

fn main() {
    let edges = split_line(&std::fs::read_to_string("edgelist.csv").unwrap());
    // let edges = vec![(0, 1), (0, 2), (1, 2), (2, 0), (2, 3), (3, 3)];
    let mut model = anf::ANF::new(edges, 7, 64);
    let N = model.compute(5);
    println!("{:?}", N);
}

fn split_line(line: &str) -> Vec<(usize, usize)> {
    line.split('\n')
        .filter(|s| !s.is_empty())
        .map(|s| {
            let (x, y) = s.trim().split_once(',').unwrap();
            (
                x.trim().parse::<usize>().unwrap(),
                y.trim().parse::<usize>().unwrap(),
            )
        })
        .collect()
}
