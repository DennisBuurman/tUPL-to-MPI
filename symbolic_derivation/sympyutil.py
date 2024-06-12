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

def extract_args(expr):
    args = []
    if len(expr.args) > 1:
        for subexpr in expr.args:
            for arg in extract_args(subexpr):
                args.append(arg)
    elif len(expr.args) == 1:
        args = [expr.args[0]]
    return args

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

def get_subs_chain(base, generated):
    chain = [base]
    prev_x, prev_expr = base
    for x, expr in generated:
        prev_x, prev_expr = (x, expr.subs(prev_x, prev_expr))
        chain.append((prev_x, prev_expr))
    return chain

if __name__ == "__main__":
    print("Utilities file for sympy deduction")