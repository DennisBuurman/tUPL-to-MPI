#!/usr/bin/env python3

import pprint
from sympy import *
from collections import defaultdict

from typing import List, Tuple, Dict

# Indexed arrays
C = IndexedBase('C') # PM_C
S = IndexedBase('S') # PM_S
P = IndexedBase('P') # PC

# Indices
n, k, i = symbols('n k i', cls=Idx)

# Indexed Equations
i_equations: Dict = {
    1 : - C[n+1] + (C[n] * S[n] - P[n]) / (S[n] - 1),
    2 : - S[n+1] + S[n] - 1,
    3 : - C[n+1] + (C[n] * S[n] + P[n]) / (S[n] + 1),
    4 : - S[n+1] + S[n] + 1
} # Indexed equations

c = Function('c') # PM_C
s = Function('s') # PM_S
p = Function('p') # PC

# Function value symbols
N, K, I = symbols('N K I', integer=True, positive=True)

# Function equations
f_equations = {
    1 : - c(N) + (c(N-1) * s(N-1) - p(N-1)) / (s(N-1) - 1),
    2 : - s(N) + s(N-1) - 1,
    3 : - c(N) + (c(N-1) * s(N-1) + p(N-1)) / (s(N-1) + 1),
    4 : - s(N) + s(N-1) + 1
} # function equations

# Map Indexed and Function notations
index_map = {
    C: c,
    S: s,
    P: p,
    n: N,
    k: K,
    i: I,
    c: C,
    s: S,
    p: P,
    N: n,
    K: k,
    I: i
}

def get_step_size(indices: List[int]) -> int:
    """ Returns step size integer if step size is equal between any two indices; else None. """
    res = None
    if len(indices) > 1:
        step = indices[1] - indices[0]
        if not (False in [True if indices[x] - indices[x-1] == step else False for x in range(1, len(indices))]):
            res = step
    return res

def get_index_dict(args: List) -> Dict:
    """ Returns a dict with IndexedBase as keys and a list of Indexes as values. """
    d: Dict = defaultdict(list)
    if type(args) == Indexed:
        args = [args]
    for a in args:
        if len(a.args) > 0 and a.func == Indexed:
            base = a.args[0]
            idx = a.args[1]
            d[base].append(idx)
    return d

def get_dict_head(d: Dict[any, List]):
    """ Takes a dict with a list as values and sets the values to the first element (head, idx 0) of the lists. """
    for key in d:
        if type(d[key]) == list:
            if len(d[key]) > 1:
                # TODO: account for this warning
                print(f"WARNING: key '{key}' contains multiple ({len(d[key])}) indices!")
            d[key] = d[key][0]
    return d

def transform_to_sum(node):
    """ Tries to transform a node to sum-notation. """
    adapted = node
    index_dict: Dict = get_index_dict(node.args)
    for base in index_dict:
        indices = index_dict[base]
        step_size = get_step_size(indices)
        # TODO: implement step sizes > 1
        if step_size == 1:
            start = min(indices)
            end = max(indices)
            # TODO: replace target with sum.doit()?
            target = Add(*[base[i] for i in indices])
            adapted = adapted.subs(target, Sum(base[k], (k, start, end)))
    return adapted

def traverse_expression(node):
    """ Traverses an expression, trying to find Add() nodes to transform into sum-notation. """
    if len(node.args) == 0:
        return node
    # check if subexpression is a sum
    if node.func == Add:
        node = transform_to_sum(node)
    # Continue traversal while keeping node structure
    return node.func(*[traverse_expression(arg) for arg in node.args])

def substitute_traversal(node, target_base: IndexedBase, target_idx: Idx, new_idx, skip: List = []):
    """ Traverses an expression, looking to substitute each target index of the provided IndexedBase with new index. """
    if len(node.args) == 0:
        return node
    # Check if node is Indexed with targeted base
    if node.func == Indexed:
        base = node.args[0]
        idx = node.args[1]
        if base == target_base and idx not in skip:
            node = node.subs(target_idx, new_idx)
    # Otherwise, check if targeted base is in a Sum object
    # TODO: check if skip list is required for Sum objects
    elif node.func == Sum:
        var = node.args[0]
        limits = node.args[1]
        if var.func == Indexed:
            if var.args[0] == target_base:
                node = node.func(var, limits.subs(target_idx, new_idx))
        else:
            # TODO: accomodate for other function types, i.e., summations where var includes mul or div operations
            print(f"WARNING: variable of sum node '{var}' contains function '{var.func}'.")
    # Continue traversal, keeping the rest in tact
    return node.func(*[substitute_traversal(arg, target_base, target_idx, new_idx, skip) for arg in node.args])

def chain_to_origin(chain: List[Tuple], origin: Tuple, start: int) -> List[Tuple]:
    """ Substitute the indices from the original expression back into the chain. """
    # Get indices of original expression
    x_idx: Dict = get_dict_head(get_index_dict(origin[0]))
    expr_idx: Dict = get_dict_head(get_index_dict(origin[1].args))

    # Substitute original indexes back into chain
    target_symbol = n # TODO: infer target symbol; either here or in parent
    for count in range(len(chain)):
        step: int = start + count
        expr = chain[count]
        # Substitute lower and upper bounds
        for target_base in x_idx:
            skip: List = []
            target_idx = x_idx[target_base].subs(target_symbol, step) # upper bound index
            new_idx = x_idx[target_base] # set the original index back as upper bound
            expr = (substitute_traversal(expr[0], target_base, target_idx, new_idx, skip), expr[1])
            skip.append(new_idx)
            
            target_idx = x_idx[target_base].subs(target_symbol, start) # lower bound index
            new_idx = target_symbol - i
            expr = (substitute_traversal(expr[0], target_base, target_idx, new_idx, skip), expr[1])

        for target_base in expr_idx:
            skip: List = []
            target_idx = expr_idx[target_base].subs(target_symbol, step) # upper bound index
            new_idx = expr_idx[target_base]
            expr = (expr[0], substitute_traversal(expr[1], target_base, target_idx, new_idx, skip))
            skip.append(new_idx)
            
            target_idx = expr_idx[target_base].subs(target_symbol, start) # lower bound index
            new_idx = target_symbol - i
            expr = (expr[0], substitute_traversal(expr[1], target_base, target_idx, new_idx, skip))
        
        chain[count] = expr
    
    return chain

def simplify_to_sum(expr, x, n, start: int = 0, end: int = None, step: int = 1):
    # Start with original expression
    # FIXME: only works for array expanding recurrences. Constant expanding recurrences, increments/decrements for example, do not work!
    origin: Tuple = (x, solve(expr, x)[0])
    generated: List[Tuple] = []
    chain: List[Tuple] = []

    # Create a base case k=start and generate n=k and n=k+1, n=k+2
    for k in range(start, start+3, step):
        next: Tuple = (x.subs(n, k), solve(expr.subs(n, k), x.subs(n, k))[0])
        generated.append(next)
    
    print(f"Generated: {generated}")
    
    # Create a substitution chain if possible
    next = generated[0]
    chain.append(next)
    for gen in generated[1:]:
        next = (gen[0], gen[1].subs(next[0], next[1]))
        chain.append(next)

    print(f"Subs chain: {chain}")

    # Traverse and transform the expressions in the chain
    chain = [(c[0], traverse_expression(c[1])) for c in chain]

    print(f"Traversed chain: {chain}")
    
    # Substitute origin indices back into chain
    chain = chain_to_origin(chain, origin, start)
    
    return chain[-1]

def indexed_deduction() -> None:
    print(">>> Step 1: solve (4) for S[n+1] <<<")
    f = i_equations[4]
    x = S[n+1]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

    print(">>> Step 2: substitute step 1 result in (3) <<<")
    f = i_equations[3].subs(res, x)
    x = P[n]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

    print(">>> Step 3: substitute C[]S[] with A[] <<<")
    # Note: composite substitution is not necessary! 
    # Solving for C[n+1]*S[n+1] is equivalent
    A = IndexedBase('A')
    f = (-x + res).subs([(C[n+1]*S[n+1], A[n+1]), (C[n]*S[n], A[n])])
    x = A[n+1]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

    print(">>> Step 4: perform iterative generation and simplify to sum <<<")
    sum_x, sum_res = simplify_to_sum(res-x, x, n)
    pprint(Eq(sum_x, sum_res))
    print("")

    print(">>> Step 5: move sum to one side and fill in n = i = 0 <<<")
    f = (-sum_x + sum_res).subs([(n, 0), (i, 0)]).doit()
    x = A[1]-A[0]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

    print(">>> Step 6: Assume A[0] = 0 for base case <<<")
    f = (-x + res).subs(A[0], 0)
    x = A[1]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

    print(">>> Step 7: substitute base case n = i = 0, A[1] = P[0] into step 4 result")
    f = (-sum_x + sum_res).subs([(-i+n, 0), (A[0], 0)])
    x = A[n+1]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

    print(">>> Step 8: substitute n with n-1 <<<")
    f = (-x + res).subs(n, n-1)
    x = A[n]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

    print(">>> Step 9: substitute A[] with C[]S[] <<<")
    f = (-x + res).subs(A[n], C[n]*S[n])
    x = C[n]*S[n]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

    print(">>> Step 10: substitute equation from step 9 into original (3) <<<")
    f = i_equations[3].subs(x, res)
    x = C[n+1]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

    print(">>> Step 11: gather terms in numerator <<<")
    # TODO: write (incremental) gather and scatter sum functions
    f = (-x + res).subs(res.args[1], Sum(P[k], (k, 0, n)))
    x = C[n+1]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

    print(">>> Step 12: transform denominator using recurrent solver on original (4) <<<")
    # FIXME: Indexed recurrence results in 'sympy.tensor.indexed.IndexException: Indexed needs at least one index.'
    # Use functional notation to use recurrent solver; note that generator must be a symbol ('n', not 'n+1')
    # TODO: map equations using function
    f_d = f_equations[4]
    x_d = s(N)
    res_d = rsolve(f_d, x_d, {s(0):0})
    if res_d in index_map:
        res_d = index_map[res_d] # TODO: map
    x_d = S[n] # TODO: map
    if res_d == x_d.args[1]:
        # TODO: account for 1: constants and 2: operations involving n
        # Goal is to deduce a sum notation if possible
        res_d = Sum(1, (k, 0, n-1))
    pprint(Eq(x_d, res_d))
    print("")

    print(">>> Step 13: substitute step 12 into step 11 result <<<")
    f = (-x + res).subs(x_d, res_d)
    x = C[n+1]
    res = solve(f, x)[0]
    pprint(Eq(x, res))
    print("")

if __name__ == "__main__":
    init_printing(use_unicode=True)
    indexed_deduction()