import csv



def read_constraints_csv(filename: str) -> tuple[set, dict]:
    with open(filename) as file:

        rdr = csv.reader(file)

        leafs = rdr.__next__()
        leafs = {leaf.strip() for leaf in leafs}

        constraints = []
        for line in enumerate(rdr):
            constraints.append(line)
    
        R = {}

        for i, constraint in enumerate(constraints):
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