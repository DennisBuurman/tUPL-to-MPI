"""
K-means core loop body operations:\n
1) PM_C[PM[i]] = (PM_C[PM[i]] * PM_S[PM[i]] - PC[i])/(PM_S[PM[i]] - 1)\n
2) PM_S[PM[i]] = PM_S[PM[i]] - 1\n
3) PM_C[m]     = (PM_C[m] * PM_S[m] + PC[i])/(PM_S[m] + 1)\n
4) PM_S[m]     = PM_S[m] + 1\n
5) PM[i]       = m
"""
from sympy import symbols, Function, init_printing

# Non-function symbols
pmc, pms, pm, pc = symbols('pmc pms pm pc', integer=True, positive=True)

# Function symbols
c = Function('c') # PM_C
s = Function('s') # PM_S
m = Function('m') # PM
p = Function('p') # PC

# Function value symbols
n = symbols('n', integer=True, positive=True)

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

if __name__ == "__main__":
    init_printing(use_unicode=True)
    print(equations)
    print(f_equations)