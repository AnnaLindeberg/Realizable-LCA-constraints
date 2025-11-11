
import networkx as nx
from lca_types import ptwo, ptwo_bin_rel
from helper_functions import get_transitive_closure




def check_realization(G: nx.DiGraph, R: ptwo_bin_rel) -> bool:
    tc_R = get_transitive_closure(R)
    for ab, cds in tc_R.items():
        for cd in cds:
            a, b = ab
            c, d = cd
            lca_ab = nx.lowest_common_ancestor(G, a, b)
            lca_cd = nx.lowest_common_ancestor(G, c, d)

            if cd not in tc_R or ab not in tc_R[cd]:    # (I1)
                if lca_ab not in nx.descendants(G, lca_cd):
                    return False
            
            else:                                       # (I2)
                if lca_ab != lca_cd:
                    return False

    return True

def main():
    from main import Algorithm_1
    from parse_input import read_constraints_csv

    X,R = read_constraints_csv("examples/figure_5.csv")
    G_r, N_r = Algorithm_1(X, R) # type: ignore
    print(check_realization(N_r, R))
    g = nx.DiGraph()
    g.add_edge("x", "xy")
    g.add_edge("y", "z")
    g.add_edge("x", "z")
    print(check_realization(g, R))


if __name__ == "__main__":
    main()