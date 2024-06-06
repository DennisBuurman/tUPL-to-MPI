#!/usr/bin/env python3

# 
# ...
# 
# Author: Dennis Buurman, Leiden University

import sys

from sympy import *
"""
1) PM_C[PM[i]] = (PM_C[PM[i]] * PM_S[PM[i]] - PC[i])/(PM_S[PM[i]] - 1)\n
2) PM_S[PM[i]] = PM_S[PM[i]] - 1\n
3) PM_C[m]     = (PM_C[m] * PM_S[m] + PC[i])/(PM_S[m] + 1)\n
4) PM_S[m]     = PM_S[m] + 1\n
5) PM[i]       = m
"""
pmc, pms, pm, pc = symbols('pmc pms pm pc', integer=True, positive=True) # non-function symbols
equations = {
    1 : - pmc + (pmc * pms - pc) / (pms - 1),
    2 : - pms + pms - 1,
    3 : - pmc + (pmc * pms + pc) / (pms + 1),
    4 : - pms + pms + 1
} # non-function equations

expressions = {
    1 : (pmc * pms - pc) / (pms - 1),
    2 : pms - 1,
    3 : (pmc * pms + pc) / (pms + 1),
    4 : pms + 1
} # non-function expressions

c = Function('c') # PM_C
s = Function('s') # PM_S
m = Function('m') # PM
p = Function('p') # PC
n = symbols('n', integer=True, positive=True)
i = symbols('i', integer=True, positive=True)

f_equations = {
    1 : - c(n) + (c(n-1) * s(n-1) - p(n-1)) / (s(n-1) - 1),
    2 : - s(n) + s(n-1) - 1,
    3 : - c(n) + (c(n-1) * s(n-1) + p(n-1)) / (s(n-1) + 1),
    4 : - s(n) + s(n-1) + 1
} # function equations

f_expressions = {
    1 : (c(n) * s(n) - p(n)) / (s(n) - 1),
    2 : s(n) - 1,
    3 : (c(n) * s(n) + p(n)) / (s(n) + 1),
    4 : s(n) + 1
} # function expressions

# TODO: rewrite print statement
def print_equation(f, x, a):
    """ Prints a function 'f' and the result 'a' of solving 'f' for 'x'."""
    print(f"f({x}): {f} = 0 -> {a}")

# TODO: rewrite print statement
def print_linear_equation(F, x, a):
    """ Prints a set of functions 'F' and the result 'a' of solving 'F' for 'x'."""
    for f in F:
        print("{", f)
    print(f"f({x}) -> {a}")

def regular_solve(debug=True):
    """ Tries to solve the equations using the sympy.solve function."""
    ans = []

    # What does it mean if pmc = pc? This is a clue for iterative versions
    f = equations[1]
    x = pmc
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_equation(f, pmc, a)

    f = equations[2]
    x = pms
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_equation(f, x, a)

    return ans

def linear_solve(debug=True):
    """ Tries to solve the expressions as sets using the sympy.solve function. """
    ans = []
    
    # Promising: expanding this with steps i might provide the desired answer
    f = [expressions[1], expressions[2]]
    x = pmc
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_linear_equation(f, x, a)
    # Promising: expanding this with steps i might result in the 'm'-level recalc.
    x = pc
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_linear_equation(f, x, a)

    # Weirdly enough, this equation set is empty
    f = [expressions[3], expressions[4]]
    x = pmc
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_linear_equation(f, x, a)

    return ans

def function_solve(debug=True):
    """ Tries to solve the f_equations using the sympy.solve function."""
    ans = []

    # This shows us that s(n) can be substituted with s(n-1)-1 
    # Notice that this only applies to equations 1 and 2, because they address using 'PM[i]'
    f = f_equations[2]
    x = s(n)
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_equation(f, x, a)

    # This shows us that s(n) can be substituted with s(n-1)+1 
    # Notice that this only applies to equations 1 and 2, because they address using 'm'
    f = f_equations[4]
    x = s(n)
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_equation(f, x, a)

    # Solving for c(n) changes nothing
    f = f_equations[3]
    x = c(n)
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_equation(f, x, a)

    # Substituting s(n-1)+1 with s(n) changes nothing when solving for c(n)
    f = f_equations[3].subs(s(n-1)+1, s(n))
    x = c(n)
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_equation(f, x, a)
    # Promising: solving for p(n-1) results in the base step I used in my deduction!
    x = p(n-1)
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_equation(f, x, a)
    
    # # Duplicate of [3] with reversed sign; may be useful to prove reverse 
    # f = f_equations[1]
    # x = c(n)
    # a = solve(f, x)
    # ans.append((f, x, a))
    # if debug:
    #     print_equation(f, x, a)

    return ans

def function_linear_solve(debug=True):
    """ Tries to solve the f_equations as sets using the sympy.solve function. """
    ans = []

    # Promising: c(n) = p(n)/s(n); Can we prove this means for [0, 1, ..., n]?
    f = [f_expressions[1], f_expressions[2]]
    x = c(n)
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_linear_equation(f, x, a)
    # Promising: p(n) = c(n)*s(n); Can we prove this means for [0, 1, ..., n]?
    x = p(n)
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_linear_equation(f, x, a)

    # Promising: c(n) = -p(n)/s(n); Can we also prove the reverse?
    f = [f_expressions[3], f_expressions[4]]
    x = c(n)
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_linear_equation(f, x, a)

    return ans

def recurrent_solve(debug=True):
    """ Tries to solve the f_equations using the sympy.rsolve function."""
    ans = []

    # Recurrent solve is useful to determine that s(n) changes from a constant i, which is the init amount
    f = f_equations[2]
    x = s(n)
    a = rsolve(f, x, {s(0): i})
    ans.append((f, x, a))
    if debug:
        print_equation(f, x, a)

    # Recurrent solve is useful to determine that s(n) changes from a constant i, which is the init amount
    f = f_equations[4]
    x = s(n)
    a = rsolve(f, x, {s(0): i})
    ans.append((f, x, a))
    if debug:
        print_equation(f, x, a)

    # FIXME: ValueError: The independent term should be a sum of hypergeometric functions, got '-p(n - 1)/(s(n - 1) - 1)'
    # f = f_equations[1]
    # x = c(n)
    # a = rsolve(f, x, {c(0): i})
    # ans.append(a)
    # if debug:
    #    print_equation(f, x, a)

    ## FIXME: ValueError: The independent term should be a sum of hypergeometric functions, got '-p(n - 1)/(s(n - 1) - 1)'
    # f = f_equations[3]
    # x = c(n)
    # a = rsolve(f, x, {c(0): i})
    # ans.append(a)
    # if debug:
    #     print_equation(f, x, a)

    # FIXME: SymPy cannot solve using multiple generators
    f = f_equations[3].subs(s(n-1), pms).subs(p(n-1), pc)
    x = c(n)
    a = rsolve(f, x) # result depends heavily on starting conditions
    ans.append(a)
    if debug:
        print_equation(f, x, a)
    
    # FIXME: SymPy cannot solve using multiple generators
    f = f_equations[3].subs(s(n-1), pms).subs(c(n), pmc).subs(c(n-1), pmc)
    x = p(n)
    a = rsolve(f, x) # result depends heavily on starting conditions
    ans.append(a)
    if debug:
        print_equation(f, x, a)

    return ans

def kmeans_deduction() -> None:
    """ Creates a deduction using write operations from the tUPL k-means 'core' loop body. """
    debug = False
    ans = []

    # Provides a starting clue (pmc = pc)
    ans += regular_solve(debug)
    # Shows promising results that may help substitution
    ans += linear_solve(debug)
    # Most useful: shows substitutions and leads to base case!
    ans += function_solve(debug)
    # Shows promising results that may help substitution
    ans += function_linear_solve(debug)
    # Only single generator equations can be handled.
    # This leaves only eq. 2 and 4
    ans += recurrent_solve(debug)

    ### Deduction steps
    
    # (1) s(n) = s(n-1) + 1 substitution
    f = f_equations[4]
    x = s(n)
    a = solve(f, x)
    print_equation(f, x, a)

    # (2) applying (1) to equation 3 and solve for p(n-1)
    f = f_equations[3].subs(s(n-1)+1, s(n))
    x = p(n-1)
    a = solve(f, x)
    ans.append((f, x, a))
    print_equation(f, x, a)

    # (3) ...

    ### END

    if debug:
        print(ans)

def main():
    init_printing(use_unicode=True)
    kmeans_deduction()

if __name__ == "__main__":
    sys.exit(main())