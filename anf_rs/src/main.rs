use std::env;
mod anf;

fn main() {
    // Prints each argument on a separate line
    if env::args().len() < 2 {
        println!("Usage: anf <filename> [distance] [r] [k]");
        return;
    }

    let distance = env::args()
        .nth(2)
        .unwrap_or("5".to_string())
        .parse::<usize>()
        .unwrap();

    let r = env::args()
        .nth(3)
        .unwrap_or("7".to_string())
        .parse::<usize>()
        .unwrap();

    let k = env::args()
        .nth(4)
        .unwrap_or("128".to_string())
        .parse::<usize>()
        .unwrap();

    let filename = env::args().nth(1).unwrap();
    let edges = split_line(&std::fs::read_to_string(filename).unwrap(), " ");
    let mut model = anf::ANF::new(edges, r, k);
    let result = model.compute(distance);

    println!("{}", serde_json::to_string(&result).unwrap());
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
