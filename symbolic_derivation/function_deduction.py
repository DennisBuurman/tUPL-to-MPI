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

from sympyutil import *
from data import *

from sympy import *
from sympy.simplify.simplify import sum_combine

i, j = symbols('i j', integer=True, positive=True)

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

def func_deduction():
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

    # (4) Inductive step: simplify to sum
    a = solve(f, A(n))
    ans = simplify_to_sum(a[0]-A(n), A(n), n, 1)
    print_equation(ans, A(n), solve(ans, A(n)))

def main():
    init_printing(use_unicode=True)
    func_deduction()

if __name__ == "__main__":
    sys.exit(main())