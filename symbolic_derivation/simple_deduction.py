from data import *
from sympyutil import print_equation, print_linear_equation

from sympy import solve, Sum, init_printing

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
    i, j = symbols('i j', integer=True, positive=True)
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

if __name__ == "__main__":
    init_printing(use_unicode=True)
    simple_deduction()