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

# TODO: break into multiple functions
def expr_to_sum(expr, x, symbol, start: int = None, end: int = None):
    """ Tries to find and substitute recurrence with sum-notation using mathematical induction. \n 
        This is done by filling in 'symbol' for [start, end] in 'expr'. \n
        returns the new expression. """
    ans = 0
    E = []
    stop = False
    i = symbols('i', integer=True, positive=True)
    j = symbols('j', integer=True, positive=True)
    k: int = start
    while not stop:
        # Fill in k and solve e_k for x_k
        e_k = expr.subs(n, k)
        x_k = x.subs(n, k)
        dummy = solve(e_k, x_k)[0]
        E.append((x_k, dummy))
        
        # Check if substitution chain results in summation
        if k > start + 2:
            e_prev = E[0] # set the first expr in E as previous
            changes = [e_prev[1]] # record changes after each subs, starting with the rhs of the first e
            for e in E[1:]:
                dummy = e[1].subs(e_prev[0], e_prev[1])
                changes.append(dummy - e_prev[1])
                e_prev = (e[0], dummy)
            # Gather changes and try to 'summify' them
            print(e_prev)
            print(changes)
            sums = {}
            for step in range(len(changes)):
                dummy = changes[step].subs(step, i)
                # TODO: make recursively go through args
                if len(dummy.args) > 1:
                    for d in dummy.args:
                        if d not in sums:
                            sums[d] = [d.subs(i, step)]
                        else:
                            sums[d].append(d.subs(i, step))
                elif d not in sums:
                    sums[d] = [d.subs(i, step)]
                else:
                    sums[d].append(d.subs(i, step))
            print(sums)
            for key in sums:
                sum_list = sums[key]
                if len(sum_list) > 1:
                    i_start = sum_list[0].args[0] # extract i from p(i)
                    i_end = sum_list[-1].args[0] # extract i from p(i)
                    step_size = sum_list[1].args[0] - sum_list[0].args[0] # difference between steps
                    equal_step_size = True
                    if len(sum_list) > 2:
                        prev = sum_list[1].args[0]
                        for s in sum_list[2:]:
                            if s.args[0] - prev != step_size:
                                equal_step_size = False
                                break
                            prev = s.args[0]
                    if step_size == 1:
                        ans = ans + Sum(key.subs(i, j), (j, i_start, i_end))
                    else:
                        print(f"WARNING: expression has step size {step_size}. Summation may be possible, but is not automation is not implemented yet!")
                        return expr
                    print(f"start, end: {i_start}, {i_end}")
                else:
                    # add singular to result expression
                    ans = ans + sum_list[0]
            ans = ans - e_prev[0]
            print(f"ANS: {ans}")
            # Now substitute the integer start with var i and the current step k with n
            ans = ans.subs([(k, n), (k-1, n-1), (start-1, n-i)])
            print(f"ANS: {ans}")
            stop = True
        k += 1
    
    print(E)

    return ""

# TODO: refine summation deduction
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

    a0 = solve(f, A(n))
    print(expr_to_sum(a0[0]-A(n), A(n), n, 1))
    return

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
    j = symbols("j", positive=True, integer=True)
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

# TODO: implement
def indexed_deduction():
    pass

def main():
    init_printing(use_unicode=True)
    # solving_trial()
    simple_deduction()
    indexed_deduction()

if __name__ == "__main__":
    sys.exit(main())