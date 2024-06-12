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

# TODO: implement
def simplify_to_sum(expr, x):
    pass

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

    # Step 4: perform inductive step
    res = simplify_to_sum(res-x, x)
    pprint(res)

if __name__ == "__main__":
    init_printing(use_unicode=True)
    indexed_deduction()