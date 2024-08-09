from data import *
from sympyutil import print_equation, print_linear_equation

from typing import List

from sympy import solve, rsolve, init_printing

def regular_solve(debug=True) -> List[any]:
    """ Tries to solve the equations using the sympy.solve function."""
    ans: List[any] = []

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

def linear_solve(debug=True) -> List[any]:
    """ Tries to solve the expressions as sets using the sympy.solve function. """
    ans: List[any] = []
    
    # Promising: expanding this with steps i might provide the desired answer
    f: List[any] = [expressions[1], expressions[2]]
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
    f: List[any] = [expressions[3], expressions[4]]
    x = pmc
    a = solve(f, x)
    ans.append((f, x, a))
    if debug:
        print_linear_equation(f, x, a)

    return ans

def function_solve(debug=True) -> List[any]:
    """ Tries to solve the f_equations using the sympy.solve function."""
    ans: List[any] = []

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

def function_linear_solve(debug=True) -> List[any]:
    """ Tries to solve the f_equations as sets using the sympy.solve function. """
    ans: List[any] = []

    # Promising: c(n) = p(n)/s(n); Can we prove this means for [0, 1, ..., n]?
    f: List[any] = [f_expressions[1], f_expressions[2]]
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

def recurrent_solve(debug=True) -> List[any]:
    """ Tries to solve the f_equations using the sympy.rsolve function."""
    ans: List[any] = []
    i = symbols('i', integer=True, positive=True)

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

if __name__ == "__main__":
    init_printing(use_unicode=True)
    regular_solve()
    linear_solve()
    function_solve()
    function_linear_solve()
    recurrent_solve()