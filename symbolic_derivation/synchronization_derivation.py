#!/usr/bin/env python3

# 
# This file contains a proof-of-concept k-means synchronization deduction.
# The idea is that deductions required to perform sync. delay transformations can be done programatically.
# If done programatically, there is a high chance that it can be automated.
# 
# Author: Dennis Buurman, Leiden University

# TODO: look at IndexedBase, Indexed, and Idx for better proof!
# https://stackoverflow.com/questions/75493900/using-sympy-to-sum-over-a-range-of-symbols

import sys

from sympy import *
from sympy.simplify.simplify import sum_combine

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
j = symbols("j", positive=True, integer=True)

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

def print_equation(f, x, a):
    """ Prints a function 'f' and the result 'a' of solving 'f' for 'x'."""
    s = f"{f} = 0 -> {x} = {a}"
    print("-"*len(s))
    print(s)
    print("-"*len(s))

def print_linear_equation(F, x, a):
    """ Prints a set of functions 'F' and the result 'a' of solving 'F' for 'x'."""
    s = f"-> {x} = {a}"
    print("-"*len(s))
    for f in F:
        print("{", f)
    print(s)
    print("-"*len(s))

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

def solving_trial() -> None:
    """ Tries to solve the write operations from the tUPL k-means 'core' loop body for different terms. """
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

    if debug:
        print(ans)

def get_subs_chain(base, generated):
    chain = [base]
    prev_x, prev_expr = base
    for x, expr in generated:
        prev_x, prev_expr = (x, expr.subs(prev_x, prev_expr))
        chain.append((prev_x, prev_expr))
    return chain

def get_chain_changes(chain):
    if len(chain) < 1:
        return []
    changes = [chain[0][1]]
    if len(chain) < 2:
        return changes
    prev = chain[0][1]
    for x in chain[1:]:
        changes.append(x[1] - prev)
        prev = x[1]
    return changes

def extract_args(expr):
    args = []
    if len(expr.args) > 1:
        for subexpr in expr.args:
            for arg in extract_args(subexpr):
                args.append(arg)
    elif len(expr.args) == 1:
        args = [expr.args[0]]
    return args

# TODO: compare expr.atoms() with expr.args extraction methods
def get_components(changes):
    components = {}
    for c in changes:
        args = extract_args(c)
        sets = [(c.args[x], args[x]) if len(c.args) > 1 else (c, args[x]) for x in range(len(args))]
        for s in sets:
            arg = s[0]
            arg_i = arg.subs([(value, i) for _, value in sets])
            if arg_i not in components:
                components[arg_i] = [arg]
            else:
                components[arg_i].append(arg)
        # for atom in c.atoms(Function):
        #     print(f"atom: {atom}")
    return components

def gather_sum(component_key, values):
    if len(values) == 1:
        return values[0]
    elif len(values) > 1:
        summable = True
        indexes = [extract_args(x)[0] for x in values]
        low = min(indexes)
        up = max(indexes)
        step = indexes[1] - indexes[0]
        if (False in [True if indexes[x] - indexes[x-1] == step else False for x in range(1, len(indexes))]):
            summable = False
        if summable:
            if step == 1:
                return Sum(component_key.subs(i, j), (j, low, up))
            else:
                print("WARNING: step size > 1 not implemented!")
                return 0
    return 0

def create_sum(origin, base, k, step_k, components):
    expr = step_k[1] - step_k[0] # base for substitution
    low = min(extract_args(expr))
    up = max(extract_args(expr))
    lower_bound = n-i
    upper_bound = n
    for key in components:
        values = components[key]
        expanded_sum = 0
        for v in values:
            expanded_sum += v
        expr = expr.subs(expanded_sum, gather_sum(key, values))
    # TODO: the following should only be done for components, not for constants
    expr = expr.subs([(low, lower_bound), (up, upper_bound)])
    expr = expr.subs([(x, n-(up-x)) for x in range(low+1, up)])
    print(f"result: {expr}, {low}, {up}")
    return expr

def simplify_to_sum(expr, x, n, start: int, end: int = None):
    result = expr
    k: int = start
    origin = (x, solve(expr, x)[0]) # original
    base = (x.subs(n, k), solve(expr.subs(n, k), x.subs(n, k))[0]) # base case with k = start

    generated = [] # list of generated expressions where k > start
    stop = False
    tries = 0 # tries to find summification
    limit = 3 # limit on tries
    while not stop:
        k += 1
        next = (x.subs(n, k), solve(expr.subs(n, k), x.subs(n, k))[0])
        generated.append(next)
        if k > start + 2:
            tries += 1
            chain = get_subs_chain(base, generated) # forward substitution chain
            changes = get_chain_changes(chain) # see what is changed/generated each iteration
            components = get_components(changes) # collect components from changes
            # Relate components back to base case in order to subs back the starting values and n/k
            result = create_sum(origin, base, k, chain[-1], components)
            stop = True
    return result

def simple_deduction():
    """ Creates a deduction using functions and substitutions. """
    # (1) s(n) = s(n-1) + 1 substitution
    f = f_equations[4]
    x = s(n)
    a = solve(f, x)
    print("(1):")
    print_equation(f, x, a)

    # (2) applying (1) to equation 3 and solve for p(n-1)
    f = f_equations[3].subs(a[0], x)
    x = p(n-1)
    a = solve(f, x)
    print("(2):")
    print_equation(f, x, a)

    # (3) Substitute c()*s() with A() to simplify
    A = Function("A") # substitute function for composite c()*s()
    f = a[0] - x
    f = f.subs(c(n)*s(n), A(n)).subs(c(n-1)*s(n-1), A(n-1))
    x = p(n-1) # solve for c(n)*s(n)
    a = solve(f, x)
    print("(3):")
    print_equation(f, x, a)

    # (4) substitute n with 1 for base case
    f1 = f.subs(n, 1)
    x1 = A(1)
    a1 = solve(f1, x1)
    print("(4):")
    print_equation(f1, x1, a1)

    # (5) substitute n with 2 for next step
    f2 = f.subs(n, 2)
    x2 = A(2)
    a2 = solve(f2, x2)
    print("(5):")
    print_equation(f2, x2, a2)

    # (6) substitute A(1) with A(0) + p(0)
    f3 = f2.subs(x1, a1[0])
    x3 = A(2)
    a3 = solve(f3, x3)
    print("(6):")
    print_equation(f3, x3, a3)

    # (7) do a third step for the fun of it
    f4 = f.subs(n, 3)
    x4 = A(3)
    a4 = solve(f4, x4)
    print("(7):")
    print_equation(f4, x4, a4)

    # (8) substitute A(2) ...
    f5 = f4.subs(x3, a3[0])
    x5 = A(3)
    a5 = solve(f5, x5)
    print("(8):")
    print_equation(f5, x5, a5)

    # (9) revert steps back to n
    f6 = f5.subs([(A(3), A(n)), (A(0), A(n-3)), (p(0), p(n-3)), (p(1), p(n-2)), (p(2), p(n-1))])
    x6 = A(n)
    a6 = solve(f6, x6)
    print("(9):")
    print_equation(f6, x6, a6)

    # (10) revert -3 back to -i simplify into summation
    psum = Sum(p(j),(j, i, n-1))
    f7 = f6.subs([(A(n-3), A(n-i)), (a6[0].subs(A(n-3), 0), psum)])
    x7 = A(n)
    a7 = solve(f7, x7)
    print("(10):")
    print_equation(f7, x7, a7)

    msg  = "Now, we only need to assume that A(0) = 0 in order to receive a summation from 0 to n - 1.\n"
    msg += "Then we would rewrite the summation to A(n-1) + p(n-1) = Sum(p(i), (i, 0, n-1)) and substitute it in the original equation.\n"
    msg += "This results in the 'i'-level variant. Substitute the summation back to the form c(n)*s(n) to receive the 'm' level variant. \n"
    msg += "As the denominator only contains one generator we can simply use rsolve, which results in 'n' === Sum(1, (i, 0, n-1))."
    print(msg)

def semi_automatic_deduction():
    """ Creates a deduction using functions and substitutions. """
    # (1) s(n) = s(n-1) + 1 substitution
    f = f_equations[4]
    x = s(n)
    a = solve(f, x)
    print_equation(f, x, a)

    # (2) applying (1) to equation 3 and solve for p(n-1)
    f = f_equations[3].subs(a[0], x)
    x = p(n-1)
    a = solve(f, x)
    print_equation(f, x, a)

    # (3) Substitute c()*s() with A() to simplify
    A = Function("A") # substitute function for composite c()*s()
    f = a[0] - x
    f = f.subs(c(n)*s(n), A(n)).subs(c(n-1)*s(n-1), A(n-1))
    x = p(n-1) # solve for c(n)*s(n)
    a = solve(f, x)
    print_equation(f, x, a)

    a0 = solve(f, A(n))
    ans = simplify_to_sum(a0[0]-A(n), A(n), n, 1)
    print_equation(ans, A(n), solve(ans, A(n)))

# TODO: implement
def indexed_deduction():
    pass

def main():
    init_printing(use_unicode=True)
    # solving_trial()
    # simple_deduction()
    semi_automatic_deduction()
    indexed_deduction()

if __name__ == "__main__":
    sys.exit(main())