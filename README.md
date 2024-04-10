# Fuzzy Kidney Solver with Patient and Donor Choice
Code modified from [this repository](https://github.com/jamestrimble/kidney_solver) 

This program is the proof of concept of using fuzzy graphs in kidney exchanges that incorporates patient and donor choice. Below are the specific descriptions of the files and what is needed to run the proof of concept.

## Prerequisites
 - Python 3
 - Gurobi

## Use
The program uses file inputs, as seen in the data folder, to input the number of patient-donor pairs, edges, and non-directed donors. Below is an introduction to the files in the repository:
 - `kidney_solver.py`: driver of the code, required arguments are the cycle-cap, chain-cap, and specific formulation
 - `kidney_ndd.py`: class for chains and finding all possible optimal chains
 - `kidney_ip.py`: specific formulation of Integer Linear Programming problem
 - `kidney_digraph.py`: inputs the kidney exchange graph and finds the score of all the cycles and chains
