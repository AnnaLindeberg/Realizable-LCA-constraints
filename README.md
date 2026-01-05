# Realizable-LCA-constraints



## Installation

The program requires Python 3.12 or higher.

#### Dependencies
* [NetworkX](https://networkx.github.io/)
* [matplotlib](https://matplotlib.org/)

## Usage and description

In a 

The program can be run using the `main.py`. A filename can be given as an extra argument when running the program, and if not the program will prompt the user for a file. 

The expected input format is a csv-file with the leafset defined on the first row, followed by one lca constraint per row. A constraint should be given as 4 values, where the first 2 values are the left hand side and the last two values are the right hand side.

So the relation R on P_2(X) with X = {a,b,c,x,y,z} and constraints abRxy, acRxy and xyRbc should be defined like this:

a,b,c,x,y,z  
a,b,x,y  
a,c,x,y  
x,y,b,c  



## Citation and references
#Todo: will this be a library? Figure out what to call it. Also find how to do proper citation.

If you use this library in your project or code, please consider citing:
  
  * __Inferring DAGs and Phylogenetic Networks from Least Common Ancestors, A. Lindeberg, A. Alfonsson, V. Moulton, G. E. Scholz, M. Hellmuth (2025)__
