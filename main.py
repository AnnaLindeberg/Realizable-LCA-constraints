# -*- coding: utf-8 -*-

import networkx as nx
import collections.abc # just for typechecking that we have things that can be used as nodes
import copy
import sys

# For printing, testing and debugging
import matplotlib.pyplot as plt

from parse_input import read_constraints_csv
from graph_drawing import draw_DAG


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
Similarily since (a,a) != (a,), we always represent {a} = {a,a} as (a,a).
"""

type ptwo_bin_rel = dict[ptwo, set[ptwo]] 
"""
The type ptwo_bin_rel represent a binary relation R over 𝒫₂(X) ⨉ 𝒫₂(X).
A relation R is represented by a dictionary D such that {a,b}R{x,y} if and only if D[(a,b)] contains (x,y).  
"""


# TODO: add a step/function to verify that the relation R actually is defined over 𝒫_2(X) ⨉ 𝒫_2(X)? Currently done in parse_input.

def get_extended_support(X: set[leaf], R: ptwo_bin_rel) -> set[ptwo]:
    """
    Input: 
        X: A set of leaf nodes.
        R: A dict representing a binaryrelation on 𝒫₂(X) ⨉ 𝒫₂(X).
    Output: 
        The extended support of R

    Given a binary relation R, get the extended support supp_plus_R.

    supp_R = { p in 𝒫₂(X) | there is some q in 𝒫₂(X) with pRq or qRp }

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
        Let S = R
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

        change_r2 = R2(S)

        change_r3 = R3(S, supp_plus_R)
        
        if not change_r2 and not change_r3:
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
            p_size = len(qs)            # Store original size to track if a change was made

            to_add: set[ptwo] = set()   # All r's such that pRqRr for some q
            for q in qs:
                if q not in S:          # To make sure we don't get key errors
                    continue
                rs = S[q]
                to_add.update(rs)
            if len(to_add) > 0:         # Needed to not add empty sets as elements when no q was related to anything
                S[p].update(to_add) 

            if p_size != len(S[p]):     # The size has a changed if something was added
                change_made = True
    return change_made


# TODO: The notation with ab and cd here doesn't match up that well with the article.
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
    for ab in S:
        for cd in S:
            # Supported combinations are all of (a,c), (a,d), (b,c), (b,d) that are in the extended support.
            # Note that the nested loop structure makes it so we also check (c, a), (d, a) etc later.
            supported_combinations = get_supported_combinations(ab, cd, supp_plus)
            if len(supported_combinations) == 0:
                continue

            overlap = S[ab].intersection(S[cd])  # All xy such that abSxy and cdSxy 
            if len(overlap) == 0:
                continue

            for p in supported_combinations:
                pre_len = len(S[p])             # Save the length before updates to see if a change was made
                S[p].update(overlap)
                if pre_len != len(S[p]):
                    change_made = True
    return change_made


def get_supported_combinations(ab: ptwo, cd: ptwo, supp_plus: set[ptwo]) -> set[ptwo]:
    """
    Input:
        ab: A tuple (a,b)
        cd: A tuple (c,d)
        supp_plus: The extended support for some relation
    Output:
        All combinations (g,h) where 
            g in (a,b) 
            and h in (c,d) 
            and (g,h) in supp_plus.
    """

    a, b = ab
    c, d = cd
    
    supported_combinations = set()
    if (a, c) in supp_plus:
        supported_combinations.add((a, c))
    if (a, d) in supp_plus:
        supported_combinations.add((a, d))
    if (b, c) in supp_plus:
        supported_combinations.add((b, c))
    if (b, d) in supp_plus:
        supported_combinations.add((b, d))
    
    return supported_combinations


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
        for (p, q) in pqs:                      # Check for each abRpq:
            if p == q and (p != a or p != b):   # if pq = pp, and ab != pp
                return False                    # then the condtion is broken
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

    tc_R = get_transitive_closure(R)
    
    for ab, xys in tc_R.items():
        for xy in xys:                                      # Check for each (ab, xy) in tc(R).
            if (xy not in tc_R) or (ab not in tc_R[xy]):    # If (xy, ab) not in tc(R),
                if (ab in R_plus[xy]):                      # and (xy, ab) in R_plus
                    return False                            # then the condition is broken.
            
    return True

def get_transitive_closure(R: ptwo_bin_rel) -> ptwo_bin_rel:
    """
    Input:
        R: A binary relation on 𝒫₂(X) ⨉ 𝒫₂(X).
    Output:
        The transitive closure of R.
    
    Makes a deepcopy of R and repeatedly applies the rule R2 on the copy to get the transitive closure.
    """

    tc_R = copy.deepcopy(R)     # Use a copy of R as base for the transitive closure
    while True:
        change_made = R2(tc_R)  # When no change is made, it's the transitive closure
        if not change_made:    
            return tc_R



def get_equiv_r_plus(R_plus: ptwo_bin_rel) -> ptwo_bin_rel:
    """
    Input:
        R_plus: The reflexive-, transitive-, 2-consistent closure of a binary relation R on 𝒫₂(X) ⨉ 𝒫₂(X).
    Output:
        An equivalence relation equiv_r_plus as described in Lemma 5.3
        
    The relation equiv_r_plus is a binary relation on the extended support of R, 
    and defined such that p equiv_r_plus q iff p R_plus q and q R_plus p.
    """

    equiv_rel: ptwo_bin_rel = {p: {p} for p in R_plus.keys()}   # R_plus is reflexive, so (p,p) is in R_plus for all p in the extended support.
                                                                # This is mainly for convienice to avoid key errors in the next part.

    for p, qs in R_plus.items():
        for q in qs:                        # For each (p,q) in R_plus:
            if p in R_plus[q]:              # if (q,p) in R_plus,
                equiv_rel[p].add(q)         # then add (p,q) to equiv_R_plus

    return equiv_rel


def get_q_set(equiv_R_plus: ptwo_bin_rel) -> set[ptwo]:
    """
    Input:
        equiv_r_plus: An equivalence relation on the extended support of some relation R.
    Output:
        A set Q_set containing one element from each equivalence class in equiv_r_plus.
    """

    to_keep = set()
    to_remove = set()
    for a in equiv_R_plus.keys():
        if a in to_remove:
            continue
        to_keep.add(a)
        bs = equiv_R_plus[a]
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

    order_R_plus = {p: {p} for p in Q_set}  # The reflexive pairs are part of the ordering, 
                                            # but will not be represented in the canonical DAG 
                                            # since we only add edges between distinct classes

    for p in Q_set:
        for q in Q_set:
            if q in R_plus[p]:
                order_R_plus[p].add(q)

    return order_R_plus

def get_canoncial_dag(order_R_plus: ptwo_bin_rel) -> nx.DiGraph:
    """
    Input:
        order_r_plus: An ordering of the equivalence classes of R_plus
    Output:
        The canonical dag of R as a networkx DiGraph

    See Definition 5.5. 
    """

    G_r = nx.DiGraph(order_R_plus).reverse()        # reverse it so we get edges q -> p instead of p -> q
                                                    # This performs step 1 and 2 in definition 5.5, but also adds a self-loop to each class
    G_r.remove_edges_from(nx.selfloop_edges(G_r))   # This removes the self-loops
    
    # Find all classes [aa]
    leaf_list = []
    for node in G_r.nodes:
        if node[0] == node[1]: #type: ignore
            leaf_list.append(node)
    
    leaf_dict = {node: node[0] for node in leaf_list} # Rename each (a,a) to a
    G_r = nx.relabel_nodes(G_r, leaf_dict)

    return G_r

def get_canoncial_network(G_r) -> nx.DiGraph:
    """
    Input:
        G_r: the canonical DAG for a relation R on 𝒫₂(X) ⨉ 𝒫₂(X)
    Output:
        The canonincal network of R as networkx DiGraph

    See definition 6.3
    """

    N_r = nx.transitive_reduction(G_r)  # Removing all shortcuts from G_r

    roots = find_roots(N_r)             # Find all roots of G_r

    if len(roots) != 1:                                                             # TODO: Can we ever have zero roots in G_r?
        for node in roots:
            N_r.add_edge("rho", node)   # If there are multiple roots, connect them all as children to a new root rho.

    return N_r

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


def main():

    # Get csv file with constraints either as commandline argument or as user input
    if len(sys.argv) == 2:
        constraint_file = sys.argv[1]
    else:
        constraint_file = input("Please write the name of a csv file defining a relation: ")

    try:
        X, R = read_constraints_csv(constraint_file)
    except ValueError as e:
        print(e)
        return
    except FileNotFoundError:
        print("The file", constraint_file, "could not be found")
        return

    res = Algorithm_1(X, R)

    if type(res) == tuple: # The type is tuple[nx.DiGraph, nx.DiGraph] if the relation was realizable
        G_r, N_r = res
    else:
        print("The relation is not realizable")
        return
    
    draw_DAG(G_r)
    draw_DAG(N_r)


if __name__ == "__main__":
    main()