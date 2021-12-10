fn main() {
    let bit_one = 0b01;
    let bit_two = 0b10;
    let bit_three = bit_one | bit_two;
    println!("{:b}", bit_three);
}

fn get_unique_nodes<T: Ord, Clone>(edges: &[(T, T)]) -> Vec<T> {
    let mut nodes = vec![];
    for &(from, to) in edges.clone().iter() {
        nodes.push(from);
        nodes.push(to);
    }
    nodes.sort();
    nodes.dedup();
    nodes
}

struct ANF<T> {
    r: u8,
    k: u8,
    edges: Vec<(T, T)>,
    nodes: Vec<T>,
}

impl<T> ANF<T> {
    fn new(edges: Vec<(T, T)>, nodes: Vec<T>, r: u8, k: u8) -> ANF<T> {
        ANF { r, k, edges, nodes }
    }
}
