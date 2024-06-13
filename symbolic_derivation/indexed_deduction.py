#!/usr/bin/env python3

import pprint
from sympy import *

# Indexed arrays
C = IndexedBase('C') # PM_C
S = IndexedBase('S') # PM_S
P = IndexedBase('P') # PC

# Indices
n, k = symbols('n k', cls=Idx)

# Indexed Equations
i_equations = {
    1 : - C[n+1] + (C[n] * S[n] - P[n]) / (S[n] - 1),
    2 : - S[n+1] + S[n] - 1,
    3 : - C[n+1] + (C[n] * S[n] + P[n]) / (S[n] + 1),
    4 : - S[n+1] + S[n] + 1
} # function equations

def traverse_sum(node):
    # TODO: gather IndexedBases and their indices
    # TODO: substitute gathered bases and indices into sum-notation
    # TODO: retain rest of node (constants)
    # TODO: return adapted if adapted != node else None
    print(f"node: {node.args}")
    return None

def traverse_expression(node):
    if len(node.args) == 0:
        return
    # check if subexpression is a sum
    if node.func == Add:
        sum = traverse_sum(node)
        if sum is not None:
            # TODO: substitute subexpression with sum
            pass
    # Continue traversal
    for arg in node.args:
        traverse_expression(arg)
    

def simplify_to_sum(expr, x, n, start: int = 0, end: int = None, step: int = 1):
    # Start with original expression
    origin = (x, solve(expr, x)[0])
    generated = []
    chain = []

    # Create a base case k = start and generate n = k and n = k+1
    for k in range(start, start+3, step):
        next = (x.subs(n, k), solve(expr.subs(n, k), x.subs(n, k))[0])
        generated.append(next)
    
    # Create a substitution chain if possible
    next = generated[0]
    chain.append(next)
    for gen in generated[1:]:
        next = (gen[0], gen[1].subs(next[0], next[1]))
        chain.append(next)

    # Traverse the expressions in the chain
    traverse_expression(chain[1][1])

    # pprint(generated)
    # pprint(chain)
    


def indexed_deduction():
    # Step 1: solve (4) for S[n+1]
    f = i_equations[4]
    x = S[n+1]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

    # Step 2: substitute step 1 result in (3)
    f = i_equations[3].subs(res, x)
    x = P[n]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

    # Step 3: substitute C[]S[] with A[]
    A = IndexedBase('A')
    f = -x + res.subs([(C[n+1]*S[n+1], A[n+1]), (C[n]*S[n], A[n])])
    x = A[n+1]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

    # Step 4: perform iterative generation to simplify to sum
    res = simplify_to_sum(res-x, x, n)
    pprint(res)

if __name__ == "__main__":
    init_printing(use_unicode=True)
    indexed_deduction()