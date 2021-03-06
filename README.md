# SNACS-road-intersections

This project contains the work done by Matthias Aarnoutse and Niels Witte regarding the final project for the Social Network Analysis for Computer Scientists (SNACS) course. This project is based on the [ANF: A Fast and Scalable Tool for Data Mining in Massive Graphs](http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf) paper and provides a more modern look at a fairly old paper.

The main goal of our work is to find out more about the road networks of countries in Europe. Data was gathered from openstreetmaps and filtered to only contain highways. The ANF algorithm is implemented in [Rust](https://www.rust-lang.org/) with the data processing being done in Python.

`main.py` can be used to analyse the road network of a single country. Run `python main.py malta` to get an analysis of Malta. See [here](https://download.geofabrik.de/europe.html) for available options

Rust build can be created by navigating to `anf_rs` and using `cargo build --release`. Replace the generated artifact with the ANF executable in the root of the repo.
