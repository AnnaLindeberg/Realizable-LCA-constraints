import csv



def read_constraints_csv(filename: str) -> tuple[set, dict]:

    with open(filename) as file:

        rdr = csv.reader(file)

        leafs = rdr.__next__()
        leafs = {leaf.strip() for leaf in leafs}
        
        constraints = []
        for line in rdr:
            constraints.append(line)
    
        R = {}
        for a, b, c, d in constraints:
            a, b, c, d = a.strip(), b.strip(), c.strip(), d.strip()
            if (a,b) not in R:
                R[(a,b)] = set()
            R[(a,b)].add((c,d))

        return leafs, R


def main():
    print(read_constraints_csv("test_file.csv"))
        
    

if __name__ == "__main__":
    main()