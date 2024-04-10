"""Solving the kidney-exchange problem using the Gurobi IP solver."""

import copy
import sys

from kidney_digraph import *
from kidney_ndds import *
import kidney_utils

from gurobipy import *


###################################################################################################
#                                                                                                 #
#                                  Code used by all formulations                                  #
#                                                                                                 #
###################################################################################################

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
        eef_alt_constraints: True if and only if alternative EEF constraints should be used
        lp_file: The name of a .lp file to write, or None if the file should not be written
        relax: True if and only if the LP relaxation should be solved also
    """

    def __init__(self, digraph, ndds, max_cycle, max_chain, verbose=False,
                 timelimit=None, edge_success_prob=1, vertex_success_prob=1,
                 lp_file=None, relax=False):
        self.digraph = digraph
        self.ndds = ndds
        self.max_cycle = max_cycle
        self.max_chain = max_chain
        self.verbose = verbose
        self.timelimit = timelimit
        self.edge_success_prob = edge_success_prob
        self.vertex_success_prob = vertex_success_prob
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

    def __init__(self, ip_model, cycles, chains, digraph, edge_success_prob=0.7, vertex_success_prob=0.7):
        self.ip_model = ip_model
        self.cycles = cycles
        self.chains = chains
        self.digraph = digraph
        self.total_score = (sum(c.score for c in chains) +
                            sum(fuzzy_cycle_score(c, digraph, edge_success_prob, vertex_success_prob) for c in cycles))
        self.edge_success_prob = edge_success_prob
        self.vertex_success_prob = vertex_success_prob

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

    # def relabelled_copy(self, old_to_new_vertices, new_digraph):
    #     """Create a copy of the solution with vertices relabelled.
    #
    #     If the solution was found on a relabelled copy of the instance digraph, this
    #     method can be used to transform the solution back to the original digraph. Each
    #     Vertex v in the OptSolution on which this method is called is replaced in the
    #     returned copy by old_to_new_vertices[v.id].
    #     """
    #
    #     relabelled_cycles = [[old_to_new_vertices[v.id] for v in c] for c in self.cycles]
    #     relabelled_chains = [Chain(c.ndd_index,
    #                                [old_to_new_vertices[i].id for i in c.vtx_indices],
    #                                c.score)
    #                          for c in self.chains]
    #     return OptSolution(self.ip_model, relabelled_cycles, relabelled_chains,
    #                        new_digraph, self.edge_success_prob)


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


def optimise_relabelled(formulation_fun, cfg):
    """Optimise on a relabelled graph such that vertices are sorted in descending
        order of (indegree + outdegree)"""

    in_degs = [0] * cfg.digraph.n
    for e in cfg.digraph.es:
        in_degs[e.tgt.id] += 1

    sorted_vertices = sorted(cfg.digraph.vs,
                             key=lambda v: len(v.edges) + in_degs[v.id],
                             reverse=True)

    relabelled_digraph = cfg.digraph.induced_subgraph(sorted_vertices)

    # old_to_new_vtx[i] is the vertex in the new graph corresponding to vertex
    # i in the original digraph
    old_to_new_vtx = [None] * cfg.digraph.n
    for i, v in enumerate(sorted_vertices):
        old_to_new_vtx[v.id] = relabelled_digraph.vs[i]

    relabelled_ndds = create_relabelled_ndds(cfg.ndds, old_to_new_vtx)
    relabelled_cfg = copy.copy(cfg)
    relabelled_cfg.digraph = relabelled_digraph
    relabelled_cfg.ndds = relabelled_ndds

    opt_result = formulation_fun(relabelled_cfg)
    return opt_result.relabelled_copy(sorted_vertices, cfg.digraph)


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
#                                       Maximum Cardinality                                       #
#                                                                                                 #
###################################################################################################

def optimise_maximum_cardinality(cfg):
    """Optimise using the cycle formulation (with one var per cycle and one var per chain).

    Args:
        cfg: an OptConfig object

    Returns:
        an OptSolution object
    """

    cycles = cfg.digraph.find_cycles(cfg.max_cycle)
    chains = find_chains(cfg.digraph, cfg.ndds, cfg.max_chain, 1, 1)

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

    obj_expr = (quicksum(cycle_score(c, cfg.digraph) * var
                         for (c, var) in zip(cycles, cycle_vars)) +
                quicksum(c.score * var for (c, var) in zip(chains, chain_vars)))

    m.setObjective(obj_expr, GRB.MAXIMIZE)
    optimise(m, cfg)

    return OptSolution(ip_model=m,
                       cycles=[c for c, v in zip(cycles, cycle_vars) if v.x > 0.5],
                       chains=[c for c, v in zip(chains, chain_vars) if v.x > 0.5],
                       digraph=cfg.digraph,
                       edge_success_prob=cfg.edge_success_prob,
                       vertex_success_prob=cfg.vertex_success_prob)



###################################################################################################
#                                                                                                 #
#                                   Failure-Aware Representation                                  #
#                                                                                                 #
###################################################################################################

def optimise_failure_aware(cfg):
    """Optimise using the cycle formulation (with one var per cycle and one var per chain).

    Args:
        cfg: an OptConfig object

    Returns:
        an OptSolution object
    """

    cycles = cfg.digraph.find_cycles(cfg.max_cycle)
    chains = find_chains(cfg.digraph, cfg.ndds, cfg.max_chain, cfg.edge_success_prob, 1)

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

    obj_expr = (quicksum(failure_aware_cycle_score(c, cfg.digraph, cfg.edge_success_prob) * var
                         for (c, var) in zip(cycles, cycle_vars)) +
                quicksum(c.score * var for (c, var) in zip(chains, chain_vars)))

    m.setObjective(obj_expr, GRB.MAXIMIZE)
    optimise(m, cfg)

    return OptSolution(ip_model=m,
                       cycles=[c for c, v in zip(cycles, cycle_vars) if v.x > 0.5],
                       chains=[c for c, v in zip(chains, chain_vars) if v.x > 0.5],
                       digraph=cfg.digraph,
                       edge_success_prob=cfg.edge_success_prob)


###################################################################################################
#                                                                                                 #
#                                   Fuzzy Graph Representation                                    #
#                                                                                                 #
###################################################################################################

def optimise_fuzzy_graph(cfg):
    """Optimise using the cycle formulation (with one var per cycle and one var per chain).

    Args:
        cfg: an OptConfig object

    Returns:
        an OptSolution object
    """

    cycles = cfg.digraph.find_cycles(cfg.max_cycle)
    chains = find_chains(cfg.digraph, cfg.ndds, cfg.max_chain, cfg.edge_success_prob, cfg.vertex_success_prob)

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

    obj_expr = (quicksum(fuzzy_cycle_score(c, cfg.digraph, cfg.edge_success_prob, cfg.vertex_success_prob) * var
                         for (c, var) in zip(cycles, cycle_vars)) +
                quicksum(c.score * var for (c, var) in zip(chains, chain_vars)))

    m.setObjective(obj_expr, GRB.MAXIMIZE)
    optimise(m, cfg)

    return OptSolution(ip_model=m,
                       cycles=[c for c, v in zip(cycles, cycle_vars) if v.x > 0.5],
                       chains=[c for c, v in zip(chains, chain_vars) if v.x > 0.5],
                       digraph=cfg.digraph,
                       edge_success_prob=cfg.edge_success_prob,
                       vertex_success_prob=cfg.vertex_success_prob)
