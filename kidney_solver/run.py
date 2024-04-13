import kidney_digraph
import kidney_ndds
import kidney_ip
import kidney_utils
import argparse


def solve_kep(cfg, formulation, use_relabelled=True):
    formulations = {
        "max": ("Maximum Cardinality", kidney_ip.optimise_max_cardinality),
        "failure": ("Failure Aware", kidney_ip.optimise_failure_aware),
        "fuzzy": ("Fuzzy Graph", kidney_ip.optimise_fuzzy)
    }

    if formulation in formulations:
        formulation_name, formulation_fun = formulations[formulation]
        if use_relabelled:
            opt_result = kidney_ip.optimise_relabelled(formulation_fun, cfg)
        else:
            opt_result = formulation_fun(cfg)
        kidney_utils.check_validity(opt_result, cfg.digraph, cfg.ndds, cfg.max_cycle, cfg.max_chain)
        opt_result.formulation_name = formulation_name
        return opt_result
    else:
        raise ValueError("Unrecognised IP formulation name")


def start():
    parser = argparse.ArgumentParser("Experiments")
    parser.add_argument("cycle_cap", type=int,
                        help="The maximum permitted cycle length")
    parser.add_argument("chain_cap", type=int,
                        help="The maximum permitted number of edges in a chain")
    parser.add_argument("formulation",
                        help="The IP formulation (max, failure, fuzzy, WC, HC)")
    parser.add_argument("--use-relabelled", "-r", required=False,
                        action="store_true",
                        help="Relabel vertices in descending order of in-deg + out-deg")
    parser.add_argument("--timelimit", "-t", required=False, default=None,
                        type=float,
                        help="IP solver time limit in seconds (default: no time limit)")
    parser.add_argument("--edge-success-prob", "-q", required=False,
                        type=float, default=1.0,
                        help="Edge success probability")
    parser.add_argument("--vertex-success-prob", "-p", required=False,
                        type=float, default=1.0,
                        help="Vertex success probability")
    parser.add_argument("--lp-file", "-l", required=False, default=None,
                        metavar='FILE',
                        help="Write the IP model to FILE, then exit.")
    parser.add_argument("--relax", "-x", required=False,
                        action='store_true',
                        help="Solve the LP relaxation.")

    args = parser.parse_args()
    args.formulation = args.formulation.lower()

    f = open(("{}.txt".format(args.formulation)), "w")

    num_runs = 200

    for run in range(num_runs):
        file1 = open(("data/data{}.txt".format(run)), 'r')
        lines = file1.readlines()
        input_lines = [line for line in lines if len(line.strip()) > 0]
        n_digraph_edges = int(input_lines[0].split()[1])
        digraph_lines = input_lines[:n_digraph_edges + 2]

        d = kidney_digraph.read_digraph(digraph_lines)

        if len(input_lines) > len(digraph_lines):
            ndd_lines = input_lines[n_digraph_edges + 2:]
            altruists = kidney_ndds.read_ndds(ndd_lines, d)
        else:
            altruists = []

        cfg = kidney_ip.OptConfig(d, altruists, args.cycle_cap, args.chain_cap, args.timelimit,
                                  args.edge_success_prob, args.vertex_success_prob,
                                  args.lp_file, args.relax)
        opt_solution = solve_kep(cfg, args.formulation, args.use_relabelled)

        f.write(str(opt_solution.total_score) + '\n')

    f.close()


if __name__ == "__main__":
    start()
