import networkx as nx
import collections.abc # just for typechecking that we have things that can be used as nodes

# For printing, testing and debugging
import matplotlib.pyplot as plt
import pprint

type node = collections.abc.Hashable
type ptwo = tuple[node,node] # TODO: These are ordered, is that fine? (Currently working on the assumpution that we don't have both ab and ba in original relation)
type ptwo_bin_rel = dict[ptwo, set[ptwo]] # TODO: Should they be dictionaries or would it be preferable to have for example a graph represent the relation?


def get_extended_support(r: ptwo_bin_rel) -> set[ptwo]: # Should this take X as parameter aswell? in this case change the second for loop. YES
    """
    Input: A dict representing a binaryrelation r (on P2(X))
    Output: The extended support of r
    """

    unique_tuples = set(r.keys())
    for vs in r.values():
        for v in vs:
            unique_tuples.add(v)
            
    out = unique_tuples.copy()
    for (u,v) in unique_tuples:
        out.add((u,u))
        out.add((v,v))
    return out


def get_r_plus(r: ptwo_bin_rel, supp_plus: set[ptwo]) -> ptwo_bin_rel:
    s: ptwo_bin_rel = r
    
    #R1
    for p in supp_plus:
        if not p in s:
            s[p] = {p}
        else:
            s[p].add(p)

    while True:
        change_made = False

        #R2
        change_made = R2(s)
        
        #R3
        change_made = R3(s, supp_plus)
        
        if not change_made:
            break           

    return s


def R2(s: ptwo_bin_rel) -> bool: # TODO: Is there a better transitive closer algorithm or similar to use?
    """
    For all p S q, 
        for all q S r, add p S r 

    repeating this until no change is made gives a transitive closure.
    """
    change_made = False
    for p, qs in s.items():
            p_size = len(qs)

            to_add: set[ptwo] = set()
            for q in qs:
                rs = s[q]
                to_add.update(rs)
            if len(to_add) > 0: # Needed to not add empty sets as elements
                s[p].update(to_add)

            if p_size != len(s[p]):
                change_made = True
    return change_made

def R3(s: ptwo_bin_rel, supp_plus: set[ptwo]) -> bool:
    change_made = False
    for (a, b) in s:
        for (c, d) in s:
            # Note that the nested loop structure makes it so we also check (c, a), (d, a) etc later.
            # Could add all 8 combinations now and do a "tringular nesting" instead of a square, but dont think it would help with anything
            supported_combinations = []
            if (a, c) in supp_plus:
                supported_combinations.append((a, c))
            if (a, d) in supp_plus:
                supported_combinations.append((a, d))
            if (b, c) in supp_plus:
                supported_combinations.append((b, c))
            if (b, d) in supp_plus:
                supported_combinations.append((b, d))
            if len(supported_combinations) == 0:
                continue

            overlap = s[(a,b)].intersection(s[(c,d)])
            if len(overlap) == 0:
                continue

            for p in supported_combinations:
                pre_len = len(s[p])
                s[p].update(overlap)
                if pre_len != len(s[p]):
                    change_made = True
    return change_made



def X1(r_plus: ptwo_bin_rel) -> bool:
    for (a, b), pqs in r_plus.items():
        for (p, q) in pqs:
            if p == q and (p != a or p != b):
                return False
    return True

def X2(r: ptwo_bin_rel, r_plus: ptwo_bin_rel) -> bool:

    r_tc = r # TODO: Has this been changed to just look at r instead of the transitive closure?
    while True:
        if not R2(r_tc):
            break
    
    for ab, xys in r_tc.items():
        for xy in xys:
            if (not ab in r_tc[xy]) and (ab in r_plus[xy]):
                return False
            
    return True


def get_equiv_r_plus(r_plus: ptwo_bin_rel) -> ptwo_bin_rel:

    equiv_rel: ptwo_bin_rel = {p: {p} for p in r_plus.keys()}

    for p, qs in r_plus.items():
        for q in qs:
            if p in r_plus[q]:
                equiv_rel[p].add(q)

    return equiv_rel


def get_q_set(equiv_r_plus: ptwo_bin_rel) -> set[ptwo]:

    to_keep = set()
    to_remove = set()
    for a in equiv_r_plus.keys():
        if a in to_remove:
            continue
        to_keep.add(a)
        bs = equiv_r_plus[a]
        to_remove.update(bs)
        
    return to_keep


def get_order_r_plus(q_set: set[ptwo], r_plus: ptwo_bin_rel) -> ptwo_bin_rel: 
    order_r_plus = {p: {p} for p in q_set} # These will all be removed again later

    for p in q_set:
        for q in q_set:
            if q in r_plus[p]:
                order_r_plus[p].add(q)

    return order_r_plus

def get_canoncial_dag(order_r_plus: ptwo_bin_rel) -> nx.DiGraph:

    # for p, qs in order_r_plus.items():
    #     qs.remove(p)

    # pprint.pp(order_r_plus)

    g_r = nx.DiGraph(order_r_plus).reverse() # reverse it so we get edges q -> p instead of p -> q
    g_r.remove_edges_from(nx.selfloop_edges(g_r)) 
    
    leaf_list = []
    for node in g_r.nodes:
        if node[0] == node[1]: #type: ignore
            leaf_list.append(node)
    
    leaf_dict = {node: node[0] for node in leaf_list}
    g_r = nx.relabel_nodes(g_r, leaf_dict)

    return g_r

def get_canoncial_network(g_r) -> nx.DiGraph:

    n_r = nx.transitive_reduction(g_r)

    roots = find_roots(n_r)

    if len(roots) != 1:
        for node in roots:
            n_r.add_edge("rho", node)

    return n_r

def find_roots(g: nx.DiGraph) -> set[ptwo]:
    roots = set()
    print(g.in_degree)
    for node, degree in g.in_degree:
        if degree == 0:
            roots.add(node)
    return roots
        

def Algorithm_1(r: ptwo_bin_rel) -> bool | tuple[nx.DiGraph, nx.DiGraph]:
    r = sort_relation(r)
    supp_plus = get_extended_support(r)                 # 1
    r_plus = get_r_plus(r, supp_plus)                   # 2
    if X1(r_plus) and X2(r, r_plus):                    # 3
        equiv_r_plus = get_equiv_r_plus(r_plus)         # 4
        q_set = get_q_set(equiv_r_plus)                 # 5
        order_r_plus = get_order_r_plus(q_set, r_plus)  # 6
        g_r = get_canoncial_dag(order_r_plus)           # 7
        n_r = get_canoncial_network(g_r)                # 8
        return g_r, n_r                                 # 9
    return False                                        # 10



def sort_relation(r: ptwo_bin_rel) -> ptwo_bin_rel:
    # Makes it so that we never have both (a, b) and (b, a) in the original relation.
    # Note, currently destructive on r
    # Could probably make an easier implementation if assuming the node elements in ptwo are comparable, which should be fine

    canonical_pairs = set()

    for (p, q) in r:
        if (q, p) in canonical_pairs:
            r[(q,p)].update(r.pop((p, q)))
        else:
            canonical_pairs.add((p, q))
    
    to_change = []
    for pq, xys in r.items():
        for (x, y) in xys:
            if (y, x) in canonical_pairs:
                to_change.append((pq, (x, y)))
            else:
                canonical_pairs.add((x, y))

    for (pq, (x, y)) in to_change:
        r[pq].remove((x, y))
        r[pq].add((y, x))

    return r




def main():
    r: ptwo_bin_rel = {
        ("x","y"): {("a","b"), ("q","t")},
        ("u","v"): {("b","a")}
    }
    r2: ptwo_bin_rel = {
        
    }
    r = sort_relation(r)
    rp = get_r_plus(r, get_extended_support(r))
    pprint.pp(rp)
    print(X1(rp))
    print(X2(r, rp))

    res = Algorithm_1(r)
    if type(res) == tuple:
        g_r, n_r = res
    else:
        return
    pos = nx.spiral_layout(g_r)
    nx.draw(g_r, pos)
    nx.draw_networkx_labels(g_r, pos)
    plt.show()

    pos = nx.spiral_layout(n_r)
    nx.draw(n_r, pos)
    nx.draw_networkx_labels(n_r, pos)
    plt.show()
    # G = nx.DiGraph(r)
    # nx.draw(G)
    # plt.show()
    
    # print(G.nodes)
    # print(type(G.nodes))

if __name__ == "__main__":
    main()