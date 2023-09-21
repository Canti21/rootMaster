import sympy as sp
from sympy import Abs
import numpy as np

MAX_IT_REACHED_MSG = 'Warning: Method failed after reaching max iterations. Result is not precise'
DIVERG_MSG = "Diverge"
ZERO_DIV = "Warning: Stoped after encounter a division by zero. Result may not be precise"
EPSILON = "Warning: Machine epsilon achieved. Execution stopped"

def newton_raphson(function, variable, x_value, tolerancy, max_iterations):
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
        np = p - (fx / dfx)
        iteration_data = {
            'iteration': i,
            'xk': extend(np),
        }
        iterations.append(iteration_data)
        if (np != 0 and ((1 + abs(np)) == 1)):
            return iterations, extend(np), EPSILON
        if (Abs(np - p) < tolerancy):
            return iterations, extend(np), None
        i = i + 1
        p = np
    return iterations, extend(p), MAX_IT_REACHED_MSG

def bisection(function, variable, a_value, b_value, tolerancy, max_iterations):
    iterations = []
    i = 0
    a = a_value
    b = b_value
    tol = tolerancy
    max_it = max_iterations
    while i < max_it:
        p = a + (b-a)/2
        fp = evaluate(function, variable, p)

        # For table generation

        iteration_data = {
            'iteration': i,
            'p': extend(p),
        }

        iterations.append(iteration_data)

        # End of table generation

        if (p != 0 and ((1 + abs(p)) == 1)):
            return iterations, extend(p), EPSILON
        if (fp == 0) or ((abs(b - a)/2) < tol):
            return iterations, extend(p), None
        i = i+1
        fa = evaluate(function, variable, a)
        if (fa.evalf() * fp.evalf()) > 0:
            a = p
            fa = fp
        else:
            b = p
    return iterations, extend(p), MAX_IT_REACHED_MSG

def secant(function, variable, a_value, b_value, tolerancy, max_iterations):
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
        if (q1 - q0) == 0:
            return iterations, extend(p), ZERO_DIV
        p = p1 - q1 * (p1 - p0)/(q1 - q0)
        iteration_data = {
            'iteration': i,
            'xk': extend(p),
        }
        iterations.append(iteration_data)
        if (p != 0 and ((1 + abs(p)) == 1)):
            return iterations, extend(p), EPSILON
        if abs(p - p1) < tolerancy:
            return iterations, extend(p), None
        i = i + 1
        p0 = p1
        q0 = q1
        p1 = p
        q1 = evaluate(function, variable, p)
    return iterations, extend(p), MAX_IT_REACHED_MSG

def fixed_point(function, variable, p0_value, tolerancy, max_iterations):
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
            'xk': extend(p),
        }
        iterations.append(iteration_data)
        if (p != 0 and ((1 + abs(p)) == 1)):
            return iterations, extend(p), EPSILON
        if abs(p - p0_value) < tolerancy:
            return iterations, extend(p), None
        i = i + 1
        p0_value = p
    return iterations, extend(p), MAX_IT_REACHED_MSG

def muller(function, variable, p0, p1, p2, tolerancy, max_iterations):
    iterations = []
    h1 = p1 - p0
    h2 = p2 - p1
    if h1 == 0 or h2 == 0:
        return None, 'Error', ZERO_DIV
    fp0 = evaluate(function, variable, p0)
    fp1 = evaluate(function, variable, p1)
    fp2 = evaluate(function, variable, p2)
    d1 = (fp1 - fp0)/h1
    d2 = (fp2 - fp1)/h2
    d = (d2 - d1)/(h2 + h1)
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
        b = d2 + (h2 * d)
        D = ((b**2) - 4 * fp2 * d)**(0.5)
        if abs(b - D) < abs (b + D):
            E = b + D
        else:
            E = b - D
        h = ((-2) * fp2) / E
        p = p2 + h
        iteration_data = {
            'iteration': i,
            'xk': extend(p),
        }

        iterations.append(iteration_data)
        if (h != 0 and ((1 + abs(h)) == 1)):
            return iterations, extend(p), EPSILON
        if abs(h) < tolerancy:
            return iterations, extend(p), None
        p0 = p1
        p1 = p2
        p2 = p
        h1 = p1 - p0
        h2 = p2 - p1
        fp0 = evaluate(function, variable, p0)
        fp1 = evaluate(function, variable, p1)
        fp2 = evaluate(function, variable, p2)
        d1 = (fp1 - fp0)/h1
        d2 = (fp2 - fp1)/h2
        d = (d2 - d1)/(h2 + h1)
        i = i + 1
    return iterations, extend(p), MAX_IT_REACHED_MSG

def calculateTolerancy(presition):
    expr = "(1/2)*(10**(-k))"
    k = "k"
    tolerancy = evaluate(expr, k, presition)
    return tolerancy

def evaluate(function, variable, value):
    # Convierte la cadena de la función en una expresión simbólica
    expr = sp.sympify(function)
    
    # Define la variable simbólica
    x = sp.symbols(variable)
    
    # Sustituye la variable por el valor proporcionado
    expr_evaluada = expr.subs(x, value)
    
    # Evalúa la expresión
    resultado = expr_evaluada.evalf()
    
    return resultado

def extend(number):
    return np.format_float_positional(number, trim='-')