"""
IP formulations for kidney exchange, including PICEF
"""

import argparse
import time
import sys

import kidney_digraph
import kidney_ip
import kidney_ndds
import kidney_utils

def solve_kep(cfg, formulation, use_relabelled=True):

    formulations = {
        "weighted":   ("Maximum Cardinality",
                  kidney_ip.optimise_weighted_choice),
        "hybrid": ("Failure Aware", kidney_ip.optimise_hybrid_choice)
    }

    if formulation in formulations:
        formulation_name, formulation_fun = formulations[formulation]
        opt_result = formulation_fun(cfg)
        kidney_utils.check_validity(opt_result, cfg.digraph, cfg.ndds, cfg.max_cycle, cfg.max_chain)
        opt_result.formulation_name = formulation_name
        return opt_result
    else:
        raise ValueError("Unrecognised IP formulation name")

def start():
    parser = argparse.ArgumentParser("Solve a kidney-exchange instance")
    parser.add_argument("cycle_cap", type=int,
            help="The maximum permitted cycle length")
    parser.add_argument("chain_cap", type=int,
            help="The maximum permitted number of edges in a chain")
    parser.add_argument("formulation",
            help="The IP formulation with choice")
    parser.add_argument("--use-relabelled", "-r", required=False,
            action="store_true",
            help="Relabel vertices in descending order of in-deg + out-deg")
    parser.add_argument("--timelimit", "-t", required=False, default=None,
            type=float,
            help="IP solver time limit in seconds (default: no time limit)")
    parser.add_argument("--verbose", "-v", required=False,
            action="store_true",
            help="Log Gurobi output to screen and log file")
    parser.add_argument("--edge-success-prob", "-q", required=False,
            type=float, default=0.7,
            help="Edge success probability")
    parser.add_argument("--vertex-success-prob", "-p", required=False,
                        type=float, default=0.7,
                        help="Edge success probability")
    parser.add_argument("--donor", required=False,
                        type=float, default=0.1,
                        help="delta for weighted and hybrid choice")
    parser.add_argument("--patient", required=False,
                        type=float, default=0.1,
                        help="lambda for weighted and hybrid choice")
    parser.add_argument("--threshold", required=False,
                        type=float, default=0.05,
                        help="threshold for choice")
    parser.add_argument("--lp-file", "-l", required=False, default=None,
            metavar='FILE',
            help="Write the IP model to FILE, then exit.")
    parser.add_argument("--relax", "-x", required=False,
            action='store_true',
            help="Solve the LP relaxation.")

    args = parser.parse_args()
    args.formulation = args.formulation.lower()

    # ftime = open("time_weighted.txt", "w")
    # fwc = open("weighted_utility.txt", "w")

    for i in range(23):
        f = open("choice_data/choice_{}.txt".format(i), "r")
        lines = f.readlines()
        input_lines = [line for line in lines if len(line.strip()) > 0]
        n_digraph_edges = int(input_lines[0].split()[1])
        digraph_lines = input_lines[:n_digraph_edges + 2]

        d = kidney_digraph.read_digraph(digraph_lines)

        if len(input_lines) > len(digraph_lines):
            ndd_lines = input_lines[n_digraph_edges + 2:]
            altruists = kidney_ndds.read_ndds(ndd_lines, d)
        else:
            altruists = []
        start_time = time.time()
        cfg = kidney_ip.OptConfig(d, altruists, args.cycle_cap, args.chain_cap, args.verbose,
                                  args.timelimit, args.edge_success_prob, args.vertex_success_prob, args.donor,
                                  args.patient, args.threshold, args.lp_file, args.relax)
        opt_weighted_solution = solve_kep(cfg, "weighted", args.use_relabelled)
        opt_hybrid_solution = solve_kep(cfg, "hybrid", args.use_relabelled)
        print(str(opt_weighted_solution.total_score) + ", " + str(opt_hybrid_solution.total_score))
        time_taken = time.time() - start_time

if __name__=="__main__":
    start()
