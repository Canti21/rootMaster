import sympy as sp
from sympy import Abs
from mathe import *

MAX_IT_REACHED_MSG = 'Warning: Method failed after reaching max iterations. Result is not precise'
ZERO_DIV = "Warning: Stopped after encounter a division by zero. Result may not be precise"
EPSILON = "Warning: Machine epsilon achieved. Execution stopped"

def newton_raphson(function, variable, x_value, tolerancy, max_iterations):
    try:
        iterations = []
        i = 0
        p = x_value
        iteration_data = {
            'iteration': 0,
            'xk': p,
        }
        iterations.append(iteration_data)
        df = sp.diff(function, variable)
        while i < max_iterations:
            fx = evaluate(function, variable, p)
            dfx = evaluate(df, variable, p)
            if dfx == 0:
                return iterations, p, ZERO_DIV
            np = sub(p, div(fx, dfx))
            iteration_data = {
                'iteration': i,
                'xk': format_non_sci_number(np),
            }
            iterations.append(iteration_data)
            if (np != 0 and ((1 + absolute(np)) == 1)):
                return iterations, format_non_sci_number(np), EPSILON
            if (Abs(sub(np, p)) < tolerancy):
                return iterations, format_non_sci_number(np), None
            i = i + 1
            p = np
        return iterations, format_non_sci_number(p), MAX_IT_REACHED_MSG
    except ValueError as e:
        return iterations, format_non_sci_number(p), e

def bisection(function, variable, a_value, b_value, tolerancy, max_iterations):
    try:
        iterations = []
        i = 0
        a = a_value
        b = b_value
        tol = tolerancy
        max_it = max_iterations
        while i < max_it:
            p = add(a, div(sub(b, a), 2))
            fp = evaluate(function, variable, p)

            # For table generation

            iteration_data = {
                'iteration': i,
                'p': format_non_sci_number(p),
            }

            iterations.append(iteration_data)

            # End of table generation

            if (p != 0 and ((1 + absolute(p)) == 1)):
                return iterations, format_non_sci_number(p), EPSILON
            if (fp == 0) or (absolute(div(sub(b, a), 2)) < tol):
                return iterations, format_non_sci_number(p), None
            i = i+1
            fa = evaluate(function, variable, a)
            if (fa.evalf() * fp.evalf()) > 0:
                a = p
                fa = fp
            else:
                b = p
        return iterations, format_non_sci_number(p), MAX_IT_REACHED_MSG
    except ValueError as e:
        return iterations, format_non_sci_number(p), e

def secant(function, variable, a_value, b_value, tolerancy, max_iterations):
    try:
        iterations = []
        i = 2
        p0 = a_value
        p1 = b_value
        iteration_data = {
            'iteration': 0,
            'xk': p0,
        }
        iterations.append(iteration_data)
        iteration_data = {
            'iteration': 1,
            'xk': p1,
        }
        iterations.append(iteration_data)
        q0 = evaluate(function, variable, p0)
        q1 = evaluate(function, variable, p1)
        while i <= max_iterations:
            if sub(q1, q0) == 0:
                return iterations, format_non_sci_number(p), ZERO_DIV
            p = sub(p1, mul(q1, div(sub(p1, p0), sub(q1, q0))))
            iteration_data = {
                'iteration': i,
                'xk': format_non_sci_number(p),
            }
            iterations.append(iteration_data)
            if (p != 0 and ((1 + absolute(p)) == 1)):
                return iterations, format_non_sci_number(p), EPSILON
            if absolute(sub(p, p1)) < tolerancy:
                return iterations, format_non_sci_number(p), None
            i = i + 1
            p0 = p1
            q0 = q1
            p1 = p
            q1 = evaluate(function, variable, p)
        return iterations, format_non_sci_number(p), MAX_IT_REACHED_MSG
    except ValueError as e:
        return iterations, format_non_sci_number(p), e

def fixed_point(function, variable, p0_value, tolerancy, max_iterations):
    try:
        iterations = []
        i = 1
        iteration_data = {
            'iteration': 0,
            'xk': p0_value,
        }
        iterations.append(iteration_data)
        while i < max_iterations:
            p = evaluate(function, variable, p0_value)
            iteration_data = {
                'iteration': i,
                'xk': format_non_sci_number(p),
            }
            iterations.append(iteration_data)
            if (p != 0 and ((1 + absolute(p)) == 1)):
                return iterations, format_non_sci_number(p), EPSILON
            if absolute(p - p0_value) < tolerancy:
                return iterations, format_non_sci_number(p), None
            i = i + 1
            p0_value = p
        return iterations, format_non_sci_number(p), MAX_IT_REACHED_MSG
    except ValueError as e:
        return iterations, format_non_sci_number(p), e

def muller(function, variable, p0, p1, p2, tolerancy, max_iterations):
    try:
        iterations = []
        h1 = sub(p1, p0)
        h2 = sub(p2, p1)
        if h1 == 0 or h2 == 0:
            return None, 'Error', ZERO_DIV
        fp0 = sp.simplify(evaluate(function, variable, p0))
        fp1 = sp.simplify(evaluate(function, variable, p1))
        fp2 = sp.simplify(evaluate(function, variable, p2))
        d1 = div(sub(fp1, fp0), h1)
        d2 = div(sub(fp2, fp1), h2)
        if add(h2, h1) == 0:
            return None, 'Error', ZERO_DIV
        d = div(sub(d2, d1), add(h2, h1))
        i = 3
        iteration_data = {
            'iteration': 0,
            'xk': p0,
        }
        iterations.append(iteration_data)
        iteration_data = {
            'iteration': 1,
            'xk': p1,
        }
        iterations.append(iteration_data)
        iteration_data = {
            'iteration': 2,
            'xk': p2,
        }
        iterations.append(iteration_data)
        while i <= max_iterations:
            b = sp.simplify(add(d2, mul(h2, d)))
            D = square_root(sub(mul(b, b), mul(4, mul(fp2, d))))
            if absolute(sub(b, D)) < absolute(add(b, D)):
                E = b + D
            else:
                E = b - D
            hmul = mul((-2), fp2)
            h = sp.simplify(div(hmul, E))
            p = add(p2, h)
            iteration_data = {
                'iteration': i,
                'xk': format_non_sci_number(p),
            }

            iterations.append(iteration_data)
            if (h != 0 and ((1 + absolute(h)) == 1)):
                return iterations, format_non_sci_number(p), EPSILON
            if absolute(h) < tolerancy:
                return iterations, format_non_sci_number(p), None
            p0 = p1
            p1 = p2
            p2 = p
            h1 = sub(p1, p0)
            h2 = sub(p2, p1)
            if h1 == 0 or h2 == 0:
                return iterations, format_non_sci_number(p), ZERO_DIV
            fp0 = sp.simplify(evaluate(function, variable, p0))
            fp1 = sp.simplify(evaluate(function, variable, p1))
            fp2 = sp.simplify(evaluate(function, variable, p2))
            d1 = div(sub(fp1, fp0), h1)
            d2 = div(sub(fp2, fp1), h2)
            if add(h2, h1) == 0:
                return iterations, format_non_sci_number(p), ZERO_DIV
            d = sp.simplify(div(sub(d2, d1), add(h2, h1)))
            i = i + 1
        return iterations, format_non_sci_number(p), MAX_IT_REACHED_MSG
    except ValueError as e:
        return iterations, format_non_sci_number(p), e

def calculateTolerancy(presition):
    expr = "(1/2)*(10**(-k))"
    k = "k"
    tolerancy = evaluate(expr, k, presition)
    return tolerancy

def evaluate(function, variable, value):
    # Converts function string into a syumbolic expresion
    expr = sp.sympify(function)
    
    # Define the symbolic variable
    x = sp.symbols(variable)

    # Check if the denominator is zero before substitution
    if sp.denom(expr).subs(x, value) == 0:
        raise ValueError("Division by zero detected.")
    
    # Substitute the variable
    expr_eval = expr.subs(x, value)
    
    # Evaluate the expresion
    result = expr_eval.evalf()
    
    return result