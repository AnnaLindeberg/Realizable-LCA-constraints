# -*- coding: utf-8 -*-

import networkx as nx
import sys

from lca_types import leaf, ptwo, ptwo_bin_rel
from parse_input import read_constraints_csv
from graph_drawing import draw_DAG
from helper_functions import (
    unify_representation,
    get_extended_support,
    get_R_plus,
    X1,
    X2,
    get_equiv_r_plus,
    get_q_set,
    get_order_r_plus,
    get_canoncial_dag,
    get_canoncial_network
)
   
        

def Algorithm_1(X: set[leaf], R: ptwo_bin_rel) -> bool | tuple[nx.DiGraph, nx.DiGraph]:
    """
    Input:
        X: Some ground set
        R: A binary relation on 𝒫₂(X) ⨉ 𝒫₂(X).
    Output:
        False if the relation is not realizable, otherwise returns the Canonical DAG and Canonical Network of R.
    """
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




def Algorithm_1_full_output(X: set[leaf], R: ptwo_bin_rel) -> tuple[bool, list]:
    """
    Input:
        X: Some ground set
        R: A binary relation on 𝒫₂(X) ⨉ 𝒫₂(X).
    Output:
        False and a list with [R, supp_plus_R, R_plus] if the relations is not realizable.
        True and a list with [R, supp_plus_R, R_plus, equiv_R_plus, Q_set, order_R_plus, G_r, N_r] if the relation is realizable.
        
    Does the same thing as Algorithm_1 but also returns the result of each step in the final output.
    """
    R = unify_representation(R)                         # Make sure the representations of elements {a,b} = {b,a} are consistent.
    supp_plus_R = get_extended_support(X, R)            # 1
    R_plus = get_R_plus(R, supp_plus_R)                 # 2
    if X1(R_plus) and X2(R, R_plus):                    # 3
        equiv_R_plus = get_equiv_r_plus(R_plus)         # 4
        Q_set = get_q_set(equiv_R_plus)                 # 5
        order_R_plus = get_order_r_plus(Q_set, R_plus)  # 6
        G_r = get_canoncial_dag(order_R_plus)           # 7
        N_r = get_canoncial_network(G_r)                # 8
        return True, [R, supp_plus_R, R_plus, equiv_R_plus, Q_set, order_R_plus, G_r, N_r]            
    return False, [R, supp_plus_R, R_plus]         

def main():
    """
    TODO: Add docstring
    """

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