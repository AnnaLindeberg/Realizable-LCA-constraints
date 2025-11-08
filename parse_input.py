import csv


"""
The expected format for a csv file is that the first line lists all of the leaves, and every line after that uses 4 values to describe a constraint

So the relation R on P_2(X) with X = {a,b,c,x,y,z} and constraints abRxy, acRxy and xyRbc would look like this:

a,b,c,x,y,z
a,b,x,y
a,c,x,y
x,y,b,c
"""


def read_constraints_csv(filename: str) -> tuple[set, dict]:
    """
    Input:
        filename: A string with the name of a csv-file describing a leafset X and a relation on 𝒫₂(X) ⨉ 𝒫₂(X).
    Output:
        The leafset as a set of strings and the relation as a dictionary R where R[p] contains q means th pRq.

    The file should have the leafs on the first line and on entry in the relation per following line.
    """
    with open(filename) as file:

        rdr = csv.reader(file)

        leafs = rdr.__next__()
        leafs = {leaf.strip() for leaf in leafs}

        R = {}

        for i, constraint in enumerate(rdr):
            if len(constraint) < 4:
                raise ValueError(f"Wrong formatting for constraint nr {i}, too few values. Constraint abRcd should be written a,b,c,d")
            if len(constraint) > 4:
                raise ValueError(f"Wrong formatting for constraint nr {i}, too many values. Constraint abRcd should be written a,b,c,d")
            
            a, b, c, d = constraint
            a, b, c, d = a.strip(), b.strip(), c.strip(), d.strip()

            if (a not in leafs 
                or b not in leafs 
                or c not in leafs 
                or d not in leafs): 
                raise ValueError(f"Constraint nr {i} mention a leaf not in the leaf set. List all leaves on the first line. Constraint abRcd should be written a,b,c,d")
            
            if (a,b) not in R:
                R[(a,b)] = set()
            R[(a,b)].add((c,d))
            

        return leafs, R


def main():
    print(read_constraints_csv("test_file.csv"))
        
    

if __name__ == "__main__":
    main()