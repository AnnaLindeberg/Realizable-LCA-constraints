# -*- coding: utf-8 -*-

import networkx as nx
import collections.abc # just for typechecking that we have things that can be used as nodes
import copy

# For printing, testing and debugging
import matplotlib.pyplot as plt
import pprint

from parse_input import read_constraints_csv


type leaf = collections.abc.Hashable
"""
Elements of the leaf set X can be of any Hashable type. 
This restriction is placed so that they will work as nodes in networkx.
"""

type ptwo = tuple[leaf,leaf] 
"""
The type ptwo represents elements in 𝒫₂(X), 1 and 2 element subsets of the leaf set X.  
Since sets are not ordered, {a, b} = {b, a} and {a,a} = {a}. 
They are implemented here as ordered pairs (tuples). 
For tuples (a,b) != (b,a). So we make sure that for each set {a,b}={b,a}, we consistently represent it either only as (a,b) or only as (b,a).   
Similarily since (a,a) != (a,), we always represent both {a} and {a,a} as (a,a).
"""

type ptwo_bin_rel = dict[ptwo, set[ptwo]] 
"""
The type ptwo_bin_rel represent a binary relation R over ptwo times ptwo.
A relation R is represented by a dictionary D such that {a,b}R{x,y} if and only if D[(a,b)] contains (x,y).  
"""


# TODO: add a step/function to verify that the relation R actually is defined over 𝒫_2(X) x P_2(X)

# TODO: Double check that no changes are made by sideeffects to R or R_plus


def get_extended_support(X: set[leaf], R: ptwo_bin_rel) -> set[ptwo]:
    """
    Input: 
        X: A set of leaf nodes.
        R: A dict representing a binaryrelation on 𝒫₂(X) ⨉ 𝒫₂(X).
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


def get_R_plus(R: ptwo_bin_rel, supp_plus_R: set[ptwo]) -> ptwo_bin_rel:
    """
    Input:
        R: A binary relation on 𝒫₂(X) ⨉ 𝒫₂(X)
        supp_plus_R: The extended support of R
    Output:
        The relexive-, transitive-, 2-consistent closure of R, called R_plus.

    This is done as described in Theorem 4.5:
        Set S = R
        R1: Add pSp for each p in supp_plus_R
        Repeatedly apply the following rules until they can no longer be applied:
            R2: if pSq and qSr, add pSr
            R3: if ab in supp_plus_R and acSxy and bdSxy for some c,d in X, add abSxy
    """
    S: ptwo_bin_rel = copy.deepcopy(R)
    
    #R1
    for p in supp_plus_R:
        if not p in S:
            S[p] = {p}
        else:
            S[p].add(p)

    while True:
        change_made = False

        #R2
        change_made = R2(S)
        
        #R3
        change_made = R3(S, supp_plus_R)
        
        if not change_made:
            break           

    return S


def R2(S: ptwo_bin_rel) -> bool: # TODO: Is there a better transitive closer algorithm or similar to use?
    """
    Input:
        S: A binary relation on 𝒫₂(X) ⨉ 𝒫₂(X).
    Output:
        True if one or more application of the rule R2 was performed, False otherwise.

    For each pSq in S, S[p] contains q. 
    For each qSr in S, S[q] contains r.
    Loops over all pair pSq in the dictionary, and for each such pair adds all r's in S[q] to S[p].

    Repeated application until False is returned gives the transitive closure of S.
    """
    change_made = False
    for p, qs in S.items():
            p_size = len(qs)

            to_add: set[ptwo] = set()
            for q in qs:
                if q not in S:
                    continue
                rs = S[q]
                to_add.update(rs)
            if len(to_add) > 0: # Needed to not add empty sets as elements
                S[p].update(to_add) 

            if p_size != len(S[p]):
                change_made = True
    return change_made

def R3(S: ptwo_bin_rel, supp_plus: set[ptwo]) -> bool:
    """
    Input:
        S: A binary relation on 𝒫₂(X) ⨉ 𝒫₂(X).
        supp_plus: The extended support of S
    Output:
        True if one or more applications of the rule R3 was performed, False otherwise.

    Looks at all pairs of keys (a,b) and (c,d) in S, 
    for each combination (m,n) among (a,c), (a,d), (b,c), (b,d) in the extended support, 
    Takes all (x,y)'s in the intersecion of S[(a,b)] and S[(c,d)] and adds (m,n)S(x,y).

    Note that this preserves the propery that any element {a,b} of 𝒫₂(X) in S has a consistent representation 
    since we only add values to keys already in the extended support, and the extended support uses consistent represenation.
    """
    change_made = False
    for (a, b) in S:
        for (c, d) in S:
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

            overlap = S[(a,b)].intersection(S[(c,d)])
            if len(overlap) == 0:
                continue

            for p in supported_combinations:
                pre_len = len(S[p])
                S[p].update(overlap)
                if pre_len != len(S[p]):
                    change_made = True
    return change_made



# TODO: Double check X1 and X2

def X1(R: ptwo_bin_rel) -> bool:
    """
    Input:
        R_plus: The reflexive-, transitive-, 2-consistent closure of a binary relation R on 𝒫₂(X) ⨉ 𝒫₂(X).
    Output:
        True if the relation R satisifies condition X1 given in definition 5.1

    The condition says that for all a,b,x in X: {a,b} != {x,x} implies ({a,b},{x,x}) not in R_plus.
    Equivalently, for all a,b,x in X: if {a,b}R_plus{x,x}, then {a,b} == {x,x}. 
    """
    for (a, b), pqs in R.items():
        for (p, q) in pqs:
            if p == q and (p != a or p != b):
                return False
    return True

def X2(R: ptwo_bin_rel, R_plus: ptwo_bin_rel) -> bool:
    """
    Input:
        R: A binary relation on 𝒫₂(X) ⨉ 𝒫₂(X).
        R_plus: The reflexive-, transitive-, 2-consistent closure of a binary relation R on 𝒫₂(X) ⨉ 𝒫₂(X).
    Output:
        True if the relation R satisifies condition X2 given in definition 5.1

    The condition says that for all a,b,x,y in X: 
        if {a,b}tc(R){x,y} but it's not the case that {x,y}tc(R){a,b}, 
        then it's not the case that {x,y}R_plus{a,b}.

    tc(R) is the tranisitive closure of R.
    """
    r_tc = copy.deepcopy(R)     # Use a copy of R as base for the transitive closure
    while True:
        change_made = R2(r_tc)  # When no change is made, it's a transitive closure
        if not change_made:    
            break
    
    for ab, xys in r_tc.items():
        for xy in xys:
            if (xy not in r_tc) or (ab not in r_tc[xy]):
                if (ab in R_plus[xy]):
                    return False
            
    return True


def get_equiv_r_plus(r_plus: ptwo_bin_rel) -> ptwo_bin_rel:
    """
    Input:
        R_plus: The reflexive-, transitive-, 2-consistent closure of a binary relation R on 𝒫₂(X) ⨉ 𝒫₂(X).
    Output:
        An equivalence relation equiv_r_plus as described in Lemma 5.3
        
    The relation equiv_r_plus is a binary relation on the extended support of R, 
    and defined such that p equiv_r_plus q iff p R_plus q and q R_plus p.
    """

    equiv_rel: ptwo_bin_rel = {p: {p} for p in r_plus.keys()}

    for p, qs in r_plus.items():
        for q in qs:
            if p in r_plus[q]:
                equiv_rel[p].add(q)

    return equiv_rel


def get_q_set(equiv_r_plus: ptwo_bin_rel) -> set[ptwo]:
    """
    Input:
        equiv_r_plus: An equivalence relation on the extended support of some relation R.
    Output:
        A set Q_set containing one element from each equivalence class in equiv_r_plus.
    """

    to_keep = set()
    to_remove = set()
    for a in equiv_r_plus.keys():
        if a in to_remove:
            continue
        to_keep.add(a)
        bs = equiv_r_plus[a]
        to_remove.update(bs)
        
    return to_keep


def get_order_r_plus(Q_set: set[ptwo], R_plus: ptwo_bin_rel) -> ptwo_bin_rel: 
    """
    Input:
        Q_set: the set of equivalence classes in equiv_r_plus
        R_plus: The reflexive-, transitive-, 2-consistent closure of a binary relation R on 𝒫₂(X) ⨉ 𝒫₂(X)
    Output:
        An ordering of the equivalnce classes in Q_set as defined in Lemma 5.4

    For two classes [p] and [q] in Q_set, we have that [p] <= [q] iff p R_plus q. 
    """

    order_r_plus = {p: {p} for p in Q_set} # These will all be removed again later

    for p in Q_set:
        for q in Q_set:
            if q in R_plus[p]:
                order_r_plus[p].add(q)

    return order_r_plus

def get_canoncial_dag(order_r_plus: ptwo_bin_rel) -> nx.DiGraph:
    """
    Input:
        order_r_plus: An ordering of the equivalence classes of R_plus
    Output:
        The canonical dag of R as a networkx DiGraph

    See Definition 5.5. 
    """

    g_r = nx.DiGraph(order_r_plus).reverse()        # reverse it so we get edges q -> p instead of p -> q
                                                    # This performs step 1 and 2 in definition 5.5, but also adds a self-loop to each class
    g_r.remove_edges_from(nx.selfloop_edges(g_r))   # This removes the self-loops
    
    # Find all classes [aa]
    leaf_list = []
    for node in g_r.nodes:
        if node[0] == node[1]: #type: ignore
            leaf_list.append(node)
    
    leaf_dict = {node: node[0] for node in leaf_list} # Rename each (a,a) to a
    g_r = nx.relabel_nodes(g_r, leaf_dict)

    return g_r

def get_canoncial_network(g_r) -> nx.DiGraph:
    """
    Input:
        G_r: the canonical DAG for a relation R on 𝒫₂(X) ⨉ 𝒫₂(X)
    Output:
        The caonincal network of R as networkx DiGraph

    See definition 6.3
    """

    n_r = nx.transitive_reduction(g_r)  # Removing all shortcuts from G_r

    roots = find_roots(n_r)             # Find all roots of G_r

    if len(roots) != 1:                                                             # TODO: Can we ever have zero roots in G_r?
        for node in roots:
            n_r.add_edge("rho", node)   # If there are multiple roots, connect them all as children to a new root rho.

    return n_r

def find_roots(G: nx.DiGraph) -> set[ptwo]:
    """
    Input:
        G: A networkx DiGraph
    Output:
        The set of root nodes in G 
    """
    roots = set()
    for node, degree in G.in_degree:
        if degree == 0:                 # A node is a root if it has indegree 0
            roots.add(node)
    return roots
        

def Algorithm_1(X: set[leaf], R: ptwo_bin_rel) -> bool | tuple[nx.DiGraph, nx.DiGraph]:
    R = unify_representation(R)                         # Make sure the representations of elements {a,b} = {b,a} are consistent.
    supp_plus_R = get_extended_support(X, R)            # 1
    R_plus = get_R_plus(R, supp_plus_R)                 # 2
    if X1(R_plus) and X2(R, R_plus):                    # 3
        equiv_R_plus = get_equiv_r_plus(R_plus)         # 4
        Q_set = get_q_set(equiv_R_plus)                 # 5
        order_R_plus = get_order_r_plus(Q_set, R_plus)  # 6
        G_r = get_canoncial_dag(order_R_plus)           # 7
        N_r = get_canoncial_network(G_r)                # 8
        return G_r, N_r                                 # 9
    return False                                        # 10


# TODO: Will the labels/nodes/leafs/groundset always be strings or some other comparable type? 
# Then will probably be easier to just say that every set {a,b} = {b,a} is represented as (a,b) iff a < b instead of this.
def unify_representation(R: ptwo_bin_rel) -> ptwo_bin_rel:
    """
    Input:
        R: a dictionary representing a binary relation over 𝒫₂(X) ⨉ 𝒫₂(X), with some sets {a,b} possibly represented as both (a,b) and (b,a).
    Output:
        Modifies R and returns it as a new representation of the relation where any element {a,b} of 𝒫₂(X) now has a consistent representation.
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



def draw_DAG(G: nx.DiGraph) -> None:
    """
    Input:
        G: A networkx DiGraph where internal nodes are tuples of leaves
    Output:
        None
    
    Draws the network with layers aligned using matplotlib and networkx multipartite layout.
    """

    for layer, nodes in enumerate(nx.topological_generations(G)):
        for node in nodes:
            G.nodes[node]["layer"] = layer

    leaf_list = []
    for node in G.nodes:
        if type(node) != tuple and node != "rho":   # When constructing the canonical dag/network for R, we relabeled all leaves to not be tuples anymore.
            leaf_list.append(node)
    leaf_layer = max(G.nodes[node]["layer"] for node in leaf_list)

    for node in leaf_list:
        G.nodes[node]["layer"] = leaf_layer

    pos = nx.multipartite_layout(G, subset_key="layer", align="horizontal")

    nx.draw(G, pos)
    nx.draw_networkx_labels(G, pos)
    plt.show()





def main():
    X, r = read_constraints_csv("test_file.csv")

    # X: set[leaf] = {"x", "y", "a", "b", "q", "t", "u", "v"}
    # r: ptwo_bin_rel = {
    #     ("x","y"): {("a","b"), ("q","t")},
    #     ("u","v"): {("b","a")}
    # }

    res = Algorithm_1(X, r)
    if type(res) == tuple:
        g_r, n_r = res
    else:
        print("Not realizable")
        return
    
    draw_DAG(g_r)
    draw_DAG(n_r)
    
    

if __name__ == "__main__":
    main()