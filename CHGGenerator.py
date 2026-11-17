from constrainthg.hypergraph import Hypergraph, Node
import constrainthg.relations as R
import numpy as np
import random

# random.seed(9 + i)
# np.random.seed(1 + i)

def prop_viability(prop: float):
    def Rprop_via(*args, **kwargs):
        return random.random() < prop
    return Rprop_via

def gen_labels(num_nodes: int):
    """Generates a list of alphabetic labels."""
    alph = range(ord('A'), ord('Z')+1)
    if num_nodes <= 26:
        labels = [chr(i) for i in alph]
    elif num_nodes <= 26**2:
        labels = [chr(i) + chr(j) for i in alph for j in alph]
    else:
        labels = [chr(i) + chr(j) + chr(k) for i in alph for j in alph for k in alph]
    labels = labels[:num_nodes]
    return labels

def max_geq(val: int):
    """Returns a method that returns True if the max argument is greater
    than or equal to `val`."""
    def Rcyclecounter(*args, **kwargs):
        args = R.extend(args, kwargs)
        return max(args) >= val
    return Rcyclecounter

def target_makes_cycle(target: str, source: str, sources_for_targets: dict, checked_nodes: list=None) -> bool:
    '''Returns true if the given target would create a cycle to the source.'''
    if checked_nodes is None:
        checked_nodes = []

    if source in checked_nodes:
        return False

    if target in sources_for_targets[source] or source == target:
        return True
    checked_nodes.append(source)

    for next_source in sources_for_targets[source]:
        if target_makes_cycle(target, next_source, sources_for_targets, checked_nodes):
            return True
    return False

def get_new_target(nodes: list, untargeted_nodes: list, targeted_nodes: list) -> Node:
    if len(untargeted_nodes) > 0:
        target = random.choice(untargeted_nodes)
        untargeted_nodes.remove(target)
        targeted_nodes.append(target)
    else:
        target = random.choice(nodes)
    return target

def add_cycle(cycle_length: int, hg: Hypergraph, nodes: list, untargeted_nodes: list, targeted_nodes: list, prop_partial_edges: float, sources_for_targets: dict) -> list:
    cycle_targets = [get_new_target(nodes, untargeted_nodes, targeted_nodes) for _ in range(cycle_length)]
    for i in range(cycle_length):
        num_sources = max(1, int(random.gammavariate(1,1)))
        other_sources = random.sample(untargeted_nodes, num_sources - 1) if len(untargeted_nodes) >= num_sources - 1 else []
        sources = [cycle_targets[i]] + other_sources
        rel = R.Rincrement if i == 0 else R.Rmax
        via = prop_viability(0.5) if random.random() < prop_partial_edges else None
        target = cycle_targets[(i + 1) % cycle_length]
        hg.add_edge(sources, target, rel, via)
        for s in sources:
            sources_for_targets[target.label].add(s.label)

        # Break out of the cycle
        req_cycle_iterations = np.random.poisson(5)
        via_out = max_geq(req_cycle_iterations)
        out_target = get_new_target(nodes, untargeted_nodes, targeted_nodes)
        hg.add_edge([cycle_targets[0]], out_target, R.Rmax, via_out)
        for s in sources:
            sources_for_targets[target.label].add(s.label)

def chg_generator(num_nodes: int, connectedness: float=1.0, prop_partial_edges: float=0.0, 
                  prop_cycles: float=0.0, prop_multiedges: float=0.0):
    """Generates a random CHG with a specified structure"""
    num_edges = int(num_nodes * connectedness)
    num_edges_remaining = num_edges
        
    hg = Hypergraph(no_weights=True)
    labels = gen_labels(num_nodes)
    nodes = [hg.add_node(l) for l in labels]

    sources_for_targets = {n.label : set() for n in nodes} ####FIXME
    targeted_nodes = []
    untargeted_nodes = [n for n in nodes] #deep copy

    # Make Cycles
    num_cycle_edges_remaining = int(prop_cycles * num_edges_remaining)
    while num_cycle_edges_remaining > 0:
        cycle_length = random.randint(1, num_cycle_edges_remaining)
        add_cycle(cycle_length, hg, nodes, untargeted_nodes, targeted_nodes, prop_partial_edges, sources_for_targets)

        num_edges_remaining -=  cycle_length + 1
        num_cycle_edges_remaining -= cycle_length

    # Make Edges
    num_acyclic_edges = num_edges_remaining
    for i in range(num_acyclic_edges):
        is_partial = random.random() < prop_partial_edges
        via = prop_viability(0.5) if is_partial else None

        while(True):
            target = get_new_target(nodes, untargeted_nodes, targeted_nodes)
            possible_sources = [n for n in nodes if not target_makes_cycle(target.label, n.label, sources_for_targets)]
            if len(possible_sources) > 0:
                break
            
        num_sources = max(1, min(int(random.gammavariate(2,1)), len(possible_sources)))
        source_nodes = random.sample(possible_sources, num_sources)

        for n in source_nodes:
            sources_for_targets[target.label].add(n.label)       
        hg.add_edge(source_nodes, target, R.Rsum, via=via)

        num_edges_remaining -= 1

    num_multiedges_remaining = max(0, int(num_nodes * prop_multiedges) - num_edges_remaining)

    # Make Multi-Edges
    for i in range(num_multiedges_remaining):
        target = random.choice(targeted_nodes)
        num_sources = max(1, int(random.gammavariate(2,1)))
        source_nodes = random.sample(nodes, num_sources)
        via = prop_viability(0.5) if random.random() < prop_partial_edges else None
        hg.add_edge(source_nodes, target, R.Rsum, via=via)

    # Add Results
    hg.targeted_nodes = targeted_nodes

    return hg

def solve_n_chgs(n: int,
                 num_nodes: int,
                 connectedness: float=0.0,
                 prop_partial_edges: float=0.0,
                 prop_cycles: float=0.0,
                 prop_multiedges: float=0.0,
                 to_print: bool=False,
        ):
    inputs = {'AA': 5, 'AB': 7}

    solved_chgs = 0
    formed_chgs = 0
    tests = 0
    num_edges= int(num_nodes * connectedness)
    
    
    
    for i in range(n):
        random.seed(9 + i)
        np.random.seed(1 + i)
        try:
            hg = chg_generator(num_nodes, connectedness, prop_partial_edges, prop_cycles, prop_multiedges)
            formed_chgs += 1
        except Exception as e:
            if to_print:
                print(e)
            continue

        try:
            target = list(hg.nodes.keys())[-1]
            if to_print:
                print(hg)
                print(hg.summary(target))
            t = hg.solve(target, inputs, to_print=to_print, search_depth=num_edges*3, logging_level=9)
            tests += 1
            if t is not None:
                solved_chgs += 1
                print("Target Value:", t.value)
            else:
                print("Target Value: None")       
        except Exception as e:
            tests += 1
            if to_print:
                print(e)
            continue

    print(f'Created {formed_chgs} CHGS.\nSolved {solved_chgs} CHGs out of {tests} tests and {num_edges} edges')

def main():
    solve_n_chgs(30, prop_multiedges=0.5, num_nodes=60, connectedness=0.5, prop_partial_edges=0.0, prop_cycles=0.0, to_print=True)


if __name__ == '__main__':
    main()

    