use std::{collections::HashMap, fmt::Debug, hash::Hash};

use rand::Rng;

fn get_unique_nodes<T: Ord + Clone>(edges: &[(T, T)]) -> Vec<T> {
    let mut nodes: Vec<T> = vec![];
    for (from, to) in edges.iter().cloned() {
        nodes.push(from);
        nodes.push(to);
    }
    nodes.sort();
    nodes.dedup();
    nodes
}

pub struct ANF {
    r: usize,
    k: usize,
    edges: Vec<(usize, usize)>,
    nodes: Vec<usize>,
}

impl ANF {
    pub fn new(edges: Vec<(usize, usize)>, r: usize, k: usize) -> ANF {
        let nodes = get_unique_nodes(&edges);
        ANF { r, k, edges, nodes }
    }

    pub fn compute(&mut self, max_distance: usize) -> Vec<f64> {
        // let M = vec![HashMap::new<T, usize>(); max_distance];
        let bitmask_len = self.bitmask_length();
        println!(
            "bitmask_len: {}, nodes: {}, edges: {}",
            bitmask_len,
            self.nodes.len(),
            self.edges.len()
        );

        let mut matrices: Vec<HashMap<usize, Vec<usize>>> = vec![];
        matrices.insert(0, HashMap::new());

        let mut individual_neigh: Vec<HashMap<&usize, f64>> = vec![];
        individual_neigh.insert(0, HashMap::new());

        for n in &self.nodes {
            (&mut matrices[0]).insert(*n, self.generate_concat_bitmasks(bitmask_len));
        }

        for h in 1..max_distance {
            // Set the value of M[h] to previous values
            let mut new_matrix = matrices[h - 1].clone();

            for (x, y) in &self.edges {
                let y_vals = matrices[h - 1].get(y).unwrap();
                for i in 0..new_matrix[x].len() {
                    new_matrix.get_mut(x).unwrap()[i] |= y_vals[i];
                }
            }

            matrices.insert(h, new_matrix);

            // approximate distrance based on the k masks
            individual_neigh.push(HashMap::new());
            for node in &self.nodes {
                individual_neigh[h].insert(node, self.approx_dist(matrices[h].get(&node).unwrap()));
            }
        }

        // return mean neighbourhood size for every distance
        let mut neigh_size: Vec<f64> = vec![self.nodes.len() as f64];
        for h in 1..max_distance {
            let mut sum = 0.0;
            // println!("{:?}", individual_neigh[h]);
            for node in &self.nodes {
                sum += individual_neigh[h].get(node).unwrap();
            }
            // mean_neigh_size.push(sum / self.nodes.len() as f64);
            neigh_size.push(sum);
        }

        return neigh_size;
    }

    fn bitmask_length(&self) -> usize {
        (self.nodes.len() as f64).log2().ceil() as usize + self.r as usize
    }

    fn generate_bitmask(&self, len: usize) -> usize {
        // generate bitmask of size len
        let mut rng = rand::thread_rng();
        let mut bitmask = 0;
        for i in 0..len {
            if rng.gen_bool(0.5_f64.powi((i + 1 as usize) as i32)) {
                bitmask |= 1 << i;
            }
        }
        bitmask
    }

    fn generate_concat_bitmasks(&self, len: usize) -> Vec<usize> {
        return (0..self.k).map(|_| self.generate_bitmask(len)).collect();
    }

    fn get_mean_least_zero_bit(&self, bitmask: &[usize]) -> f64 {
        let mut mean: f64 = 0_f64;
        for b in bitmask {
            mean += b.trailing_ones() as f64;
        }
        return mean / self.k as f64;
    }

    fn approx_dist(&self, bitmask: &[usize]) -> f64 {
        return 2_i64.pow(self.get_mean_least_zero_bit(bitmask) as u32) as f64 / 0.77351_f64;
    }
}
