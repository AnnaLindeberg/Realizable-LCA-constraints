import networkx as nx
import collections.abc # just for typechecking that we have things that can be used as nodes

# For printing, testing and debugging
import matplotlib.pyplot as plt
import pprint


type leaf = collections.abc.Hashable
"""
Elements of the leaf set X can be of any Hashable type. 
This restriction is placed so that they will work as nodes in networkx.
"""

type ptwo = tuple[leaf,leaf] 
"""
The type ptwo represents 1 and 2 element subsets of the leaf set X.  
Since sets are not ordered, {a, b} = {b, a} and {a,a} = {a}. 
They are implemented here as ordered pairs (tuples). 
For tuples (a,b) != (b,a). So we make sure that for each set {a,b}={b,a}, we consistently represent it either only as (a,b) or only as (b,a).   
Similarily since (a,a) != (a,), we always represent (a,) as (a,a).
"""

type ptwo_bin_rel = dict[ptwo, set[ptwo]] 
"""
The type ptwo_bin_rel represent a binary relation R over ptwo times ptwo.
A relation R is represented by a dictionary D such that {a,b}R{x,y} if and only if D[(a,b)] contains (x,y).  
"""


# TODO: add a step/function to verify that the relation R actually is defined over the squared powerset of X



def get_extended_support(X: set[leaf], R: ptwo_bin_rel) -> set[ptwo]:
    """
    Input: 
        X: A set of leaf nodes.
        R: A dict representing a binaryrelation on P_2(X) x P_2(X).
    Output: 
        The extended support of R

    Given a binary relation R, get the extended support supp_plus_R.

    supp_R = { (p,q) in P_2(X) | there is some (a,b) in P_2(X) with (p,q)R(a,b) or (a,b)R(p,q)}

    supp_plus_R = supp_R union {(x,x) | x in X}
    """

    out = set(R.keys())
    for qs in R.values():
        for q in qs:
            out.add(q)

    # At this point out = supp_R

    for x in X:
        out.add((x,x))

    # At this point out = supp_plus_R

    return out


def get_r_plus(R: ptwo_bin_rel, supp_plus: set[ptwo]) -> ptwo_bin_rel:
    s: ptwo_bin_rel = R
    
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

    r_tc = r 
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

def find_roots(G: nx.DiGraph) -> set[ptwo]:
    """
    Input:
        G: A networkx DiGraph
    Output:
        The set of root nodes in G 
    """
    roots = set()
    print(G.in_degree)
    for node, degree in G.in_degree:
        if degree == 0:
            roots.add(node)
    return roots
        

def Algorithm_1(X: set[leaf], R: ptwo_bin_rel) -> bool | tuple[nx.DiGraph, nx.DiGraph]:
    R = unify_representation(R)
    supp_plus_R = get_extended_support(X, R)            # 1
    R_plus = get_r_plus(R, supp_plus_R)                 # 2
    if X1(R_plus) and X2(R, R_plus):                    # 3
        equiv_R_plus = get_equiv_r_plus(R_plus)         # 4
        Q_set = get_q_set(equiv_R_plus)                 # 5
        order_R_plus = get_order_r_plus(Q_set, R_plus)  # 6
        G_r = get_canoncial_dag(order_R_plus)           # 7
        N_r = get_canoncial_network(G_r)                # 8
        return G_r, N_r                                 # 9
    return False                                        # 10



def unify_representation(R: ptwo_bin_rel) -> ptwo_bin_rel:
    """
    Input:
        R: a dictionary representing a binary relation over P_2(X) x P_2(X), with some sets {a,b} possibly represented as both (a,b) and (b,a).
    Output:
        Modifies R and returns it as a new representation of the relation where any element {a,b} of P_2(X) now has a consistent representation.
    """

    # When a representation (a,b) or (b,a) for the each set {a,b} is first seen, save it as canoncial.
    canonical_pairs = set()

    # Goes over the keys (left-hand sides of the relation) 
    # and merges any pairs of keys (a,b) and (b,a) into the canonical key.
    for (p, q) in R:
        if (q, p) in canonical_pairs and p != q:    # If the opposite ordering is already canoncial
            R[(q,p)].update(R.pop((p, q)))          # adds all elements of R[(p,q)] to R[(q, p)] and removes R[(p,q)]
        else:
            canonical_pairs.add((p, q))
    

    # Goes over all values (right-hand sides of the relation) 
    # and saves any entries (p,q)R(x,y) where (y,x) has already been decided as canoncial.
    to_change = []
    for pq, xys in R.items():
        for (x, y) in xys:
            if (y, x) in canonical_pairs and x != y:
                to_change.append((pq, (x, y)))
            else:
                canonical_pairs.add((x, y))

    # for each saved entry (p,q)R(x,y), remove it and add the canonical representation (p,q)R(y,x)
    for (pq, (x, y)) in to_change:
        R[pq].remove((x, y))
        R[pq].add((y, x))

    return R




def main():
    X: set[leaf] = {"x", "y", "a", "b", "q", "t", "u", "v"}
    r: ptwo_bin_rel = {
        ("x","y"): {("a","b"), ("q","t")},
        ("u","v"): {("b","a")}
    }
    r2: ptwo_bin_rel = {
        
    }
    r = unify_representation(r)
    rp = get_r_plus(r, get_extended_support(X, r))
    pprint.pp(rp)
    print(X1(rp))
    print(X2(r, rp))

    res = Algorithm_1(X, r)
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