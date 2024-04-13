from kidney_digraph import *
from kidney_ndds import *

from gurobipy import *


class OptConfig(object):
    """The inputs (problem instance and parameters) for an optimisation run

    Data members:
        digraph
        ndds
        max_cycle
        max_chain
        verbose: True if and only if Gurobi output should be writtent to screen and log file
        timelimit
        edge_success_prob
        vertex_success_prob
        donor
        patient
        threshold
        lp_file: The name of a .lp file to write, or None if the file should not be written
        relax: True if and only if the LP relaxation should be solved also
    """

    def __init__(self, digraph, ndds, max_cycle, max_chain, verbose=False,
                 timelimit=None, edge_success_prob=1, donor=0.1, patient=0.1,
                 vertex_success_prob=1, threshold=0.1, lp_file=None, relax=False):
        self.digraph = digraph
        self.ndds = ndds
        self.max_cycle = max_cycle
        self.max_chain = max_chain
        self.verbose = verbose
        self.timelimit = timelimit
        self.edge_success_prob = edge_success_prob
        self.vertex_success_prob = vertex_success_prob
        self.patient = patient
        self.donor = donor
        self.threshold = threshold
        self.lp_file = lp_file
        self.relax = relax


class OptSolution(object):
    """An optimal solution for a kidney-exchange problem instance.

    Data members:
        ip_model: The Gurobi Model object
        cycles: A list of cycles in the optimal solution, each represented
            as a list of vertices
        chains: A list of chains in the optimal solution, each represented
            as a Chain object
        total_score: The total score of the solution
    """

    def __init__(self, ip_model, cycles, chains, digraph, edge_success_prob=0.7, vertex_success_prob=0.7, donor=0.1,
                 patient=0.1, threshold=0.05):
        self.ip_model = ip_model
        self.cycles = cycles
        self.chains = chains
        self.digraph = digraph
        self.total_score = (sum(c.expected_score for c in chains) +
                            sum(expected_utility(c, digraph, edge_success_prob, vertex_success_prob) for c in cycles))
        self.edge_success_prob = edge_success_prob
        self.vertex_success_prob = vertex_success_prob
        self.donor = donor
        self.patient = patient
        self.threshold = threshold

    def display(self):
        """Print the optimal cycles and chains to standard output."""

        print(("cycle_count: {}".format(len(self.cycles))))
        print(("chain_count: {}".format(len(self.chains))))
        print("cycles:")
        # cs is a list of cycles, with each cycle represented as a list of vertex IDs
        cs = [[v.id for v in c] for c in self.cycles]
        # Put the lowest-indexed vertex at the start of each cycle
        for i in range(len(cs)):
            min_index_pos = cs[i].index(min(cs[i]))
            cs[i] = cs[i][min_index_pos:] + cs[i][:min_index_pos]
        # Sort the cycles
        cs.sort()
        for c in cs:
            print(("\t".join(str(v_id) for v_id in c)))
        print("chains:")
        for c in self.chains:
            print((str(c.ndd_index) + "\t" + "\t".join(str(v) for v in c.vtx_indices)))


def optimise(model, cfg):
    if cfg.lp_file:
        model.update()
        model.write(cfg.lp_file)
        sys.exit(0)
    elif cfg.relax:
        model.update()
        r = model.relax()
        r.optimize()
        print(("lp_relax_obj_val:", r.obj_val))
        print(("lp_relax_solver_status:", r.status))
        sys.exit(0)
    else:
        model.optimize()


def create_ip_model(time_limit, verbose):
    """Create a Gurobi Model."""

    m = Model("kidney-mip")
    if not verbose:
        m.params.outputflag = 0
    m.params.mipGap = 0
    if time_limit is not None:
        m.params.timelimit = time_limit
    return m


###################################################################################################
#                                                                                                 #
#                                        Weighted Choice                                          #
#                                                                                                 #
###################################################################################################

def optimise_weighted_choice(cfg):
    """Optimise using the cycle formulation (with one var per cycle and one var per chain).

    Args:
        cfg: an OptConfig object

    Returns:
        an OptSolution object
    """

    cycles = cfg.digraph.find_cycles(cfg.max_cycle)
    chains = weighted_choice_find_chains(cfg.digraph, cfg.ndds, cfg.max_chain, cfg.edge_success_prob,
                                         cfg.vertex_success_prob, cfg.donor, cfg.patient)

    m = create_ip_model(cfg.timelimit, cfg.verbose)
    m.params.method = 2

    cycle_vars = [m.addVar(vtype=GRB.BINARY) for __ in cycles]
    chain_vars = [m.addVar(vtype=GRB.BINARY) for __ in chains]
    m.update()

    ndd_to_vars = [[] for __ in cfg.ndds]
    vtx_to_vars = [[] for __ in cfg.digraph.vs]

    for var, c in zip(cycle_vars, cycles):
        for v in c:
            vtx_to_vars[v.id].append(var)

    for var, c in zip(chain_vars, chains):
        ndd_to_vars[c.ndd_index].append(var)
        for v in c.vtx_indices:
            vtx_to_vars[v].append(var)

    # Each donor-patient pair and each each NDD is in at most one chosen cycle or chain
    for l in vtx_to_vars + ndd_to_vars:
        if len(l) > 0:
            m.addConstr(quicksum(l) <= 1)

    obj_expr = (quicksum(weighted_choice_utility(c, cfg.digraph) * var
                         for (c, var) in zip(cycles, cycle_vars)) +
                quicksum(c.score * var for (c, var) in zip(chains, chain_vars)))

    m.setObjective(obj_expr, GRB.MAXIMIZE)
    optimise(m, cfg)

    return OptSolution(ip_model=m,
                       cycles=[c for c, v in zip(cycles, cycle_vars) if v.x > 0.5],
                       chains=[c for c, v in zip(chains, chain_vars) if v.x > 0.5],
                       digraph=cfg.digraph,
                       edge_success_prob=cfg.edge_success_prob,
                       vertex_success_prob=cfg.vertex_success_prob,
                       donor=cfg.donor,
                       patient=cfg.patient)

###################################################################################################
#                                                                                                 #
#                                        Hybrid Choice                                            #
#                                                                                                 #
###################################################################################################

def optimise_hybrid_choice(cfg):
    """Optimise using the cycle formulation (with one var per cycle and one var per chain).

    Args:
        cfg: an OptConfig object

    Returns:
        an OptSolution object
    """

    cycles = cfg.digraph.find_cycles(cfg.max_cycle)
    chains = hybrid_choice_find_chains(cfg.digraph, cfg.ndds, cfg.max_chain, cfg.edge_success_prob,
                                         cfg.vertex_success_prob, cfg.donor, cfg.patient, cfg.threshold)

    cycles = [cycle for cycle in cycles if choice_utility(cycle, cfg.digraph) >= cfg.threshold]

    m = create_ip_model(cfg.timelimit, cfg.verbose)
    m.params.method = 2

    cycle_vars = [m.addVar(vtype=GRB.BINARY) for __ in cycles]
    chain_vars = [m.addVar(vtype=GRB.BINARY) for __ in chains]
    m.update()

    ndd_to_vars = [[] for __ in cfg.ndds]
    vtx_to_vars = [[] for __ in cfg.digraph.vs]

    for var, c in zip(cycle_vars, cycles):
        for v in c:
            vtx_to_vars[v.id].append(var)

    for var, c in zip(chain_vars, chains):
        ndd_to_vars[c.ndd_index].append(var)
        for v in c.vtx_indices:
            vtx_to_vars[v].append(var)

    # Each donor-patient pair and each each NDD is in at most one chosen cycle or chain
    for l in vtx_to_vars + ndd_to_vars:
        if len(l) > 0:
            m.addConstr(quicksum(l) <= 1)

    obj_expr = (quicksum(expected_utility(c, cfg.digraph) * var
                         for (c, var) in zip(cycles, cycle_vars)) +
                quicksum(c.expected_score * var for (c, var) in zip(chains, chain_vars)))

    m.setObjective(obj_expr, GRB.MAXIMIZE)
    optimise(m, cfg)

    return OptSolution(ip_model=m,
                       cycles=[c for c, v in zip(cycles, cycle_vars) if v.x > 0.5],
                       chains=[c for c, v in zip(chains, chain_vars) if v.x > 0.5],
                       digraph=cfg.digraph,
                       edge_success_prob=cfg.edge_success_prob,
                       vertex_success_prob=cfg.vertex_success_prob,
                       donor=cfg.donor,
                       patient=cfg.patient,
                       threshold=cfg.threshold)
