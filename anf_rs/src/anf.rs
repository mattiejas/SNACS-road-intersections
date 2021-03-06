use serde::Serialize;
use std::collections::HashMap;

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

#[derive(Serialize)]
pub struct Result {
    r: usize,
    k: usize,
    max_distance: usize,
    neighbourhoods: Vec<f64>,
    mean_individual_neighbourhoods: Vec<f64>,
    time_in_ms: f32,
    bitmask_length: usize,
    nodes: usize,
    edges: usize,
    individual_neighbourhoods: Option<Vec<HashMap<usize, f64>>>,
}

impl ANF {
    pub fn new(edges: Vec<(usize, usize)>, r: usize, k: usize) -> ANF {
        let nodes = get_unique_nodes(&edges);
        ANF { r, k, edges, nodes }
    }

    pub fn compute(&mut self, max_distance: usize, show_in: bool) -> Result {
        let now = std::time::Instant::now();

        let bitmask_len = self.bitmask_length();
        let mut matrices: Vec<HashMap<usize, Vec<usize>>> = vec![];
        matrices.insert(0, HashMap::new());

        let mut individual_neigh: Vec<HashMap<usize, f64>> = vec![];
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
                individual_neigh[h].insert(
                    node.clone(),
                    self.approx_dist(matrices[h].get(&node).unwrap()),
                );
            }
        }

        // return mean neighbourhood size for every distance
        let mut neigh_size: Vec<f64> = vec![self.nodes.len() as f64];
        let mut mean_neigh_size: Vec<f64> = vec![1.0];

        for h in 1..max_distance {
            let mut sum = 0.0;
            // println!("{:?}", individual_neigh[h]);
            for node in &self.nodes {
                sum += individual_neigh[h].get(node).unwrap();
            }
            mean_neigh_size.push(sum / self.nodes.len() as f64);
            neigh_size.push(sum);
        }

        let time = (now.elapsed().as_nanos()) as f32 / 1000000f32;

        return Result {
            r: self.r,
            k: self.k,
            max_distance,
            neighbourhoods: neigh_size,
            mean_individual_neighbourhoods: mean_neigh_size,
            time_in_ms: time,
            bitmask_length: bitmask_len,
            nodes: self.nodes.len(),
            edges: self.edges.len(),
            individual_neighbourhoods: if show_in {
                Some(individual_neigh)
            } else {
                None
            },
        };
    }

    fn bitmask_length(&self) -> usize {
        (self.nodes.len() as f64).log2().ceil() as usize + self.r as usize
    }

    fn generate_bitmask(&self, len: usize) -> usize {
        // generate bitmask of size len
        let mut rng = rand::thread_rng();
        let mut sum: f64 = 0.;
        let p: f64 = rng.gen();

        for i in 0..len {
            sum += 0.5_f64.powi((i + 1) as i32);

            if p < sum {
                return 1 << i;
            }
        }
        return 0;
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
        return 2_f64.powf(self.get_mean_least_zero_bit(bitmask)) / 0.77351_f64;
    }
}
