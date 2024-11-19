# Fuzzy Kidney Solver with Patient and Donor Choice
Code modified from [this repository](https://github.com/jamestrimble/kidney_solver) 

This program is the proof of concept of using fuzzy graphs in kidney exchanges that incorporates patient and donor choice. Below are the specific descriptions of the files and what is needed to run the proof of concept. This project was accepted for presentation at the American Association for the Advancement of Science (AAAS) Annual Meeting 2025.

## Prerequisites
 - Python 3
 - Gurobi

## Use
In each folder, there is a file called `kidney_solver.py` which is used to run the optimization algorithm. Below is a general description for each of the files in each of the folders.
 - `kidney_solver.py`: driver of the code, required arguments are the cycle-cap, chain-cap, and specific formulation
 - `kidney_ndd.py`: class for chains and finding all possible optimal chains
 - `kidney_ip.py`: specific formulation of Integer Linear Programming problem
 - `kidney_digraph.py`: inputs the kidney exchange graph and finds the score of all the cycles and chains

The folder that is simply named `kidney_solver` contains the implementations of the fuzzy graph representation, failure aware representation, and maximum cardinality. The folder that is named `choice_kidney_solver` contains weighted choice and hybrid choice. 

## Commands
Going into each of the folders, call the below command to run the program.
```
python kidney_solver.py 3 5 max
```
which runs `kidney_solver.py` with a cycle cap of 3 and a chain cap of 5, running maximum cardinality.

