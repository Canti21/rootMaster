import sympy as sp
import numpy as np
import cmath as cm

def format_non_sci_number(number):
    """
    Format a number as a string without scientific notation, handling complex numbers.

    This function takes a numeric input 'number' and formats it as a string without
    scientific notation. If 'number' is a complex number, it separates and formats
    its real and imaginary parts as strings. Otherwise, it formats the input as a
    regular number.

    Parameters:
    number (numeric or complex): The input number to be formatted.

    Returns:
    str: A string representation of the formatted number.

    Example:
    >>> format_non_sci_number(3.14159265359)
    '3.14159265359'

    >>> format_non_sci_number(2 + 3j)
    '2.00000000000000 + 3.00000000000000i'
    """
    # Check if the number is complex
    if isinstance(number, complex) or sp.I in sp.sympify(number).atoms():
        # If its complex, then remove scientific notation from real and imaginary part
        real_part = np.format_float_positional(sp.re(number), trim='-')
        imag_coeff = sp.im(number)
        imag_part = np.format_float_positional(abs(imag_coeff), trim='-')

        # Determine the sign of the imaginary part
        sign = '-' if imag_coeff < 0 else '+'
        return f"{real_part} {sign} {imag_part}i"
    else:
        # Remove scientific notation
        return np.format_float_positional(number, trim='-')

def absolute(number):
    if isinstance(number, complex) or sp.I in sp.sympify(number).atoms():
        real_part = sp.re(number)
        imag_coeff = sp.im(number)
        imag_part = abs(imag_coeff)
        sign = '-' if imag_coeff < 0 else '+'
        absol = np.sqrt(float((real_part * real_part) + (imag_part * imag_part)))
    else:
        absol = abs(number)
    return absol

def add(a, b):
    """
    Add two numbers, handling complex and real values.

    Parameters:
    a (numeric or complex): The first number.
    b (numeric or complex): The second number.

    Returns:
    numeric or complex: The result of adding 'a' and 'b'.
    """
    if isinstance(a, complex) or sp.I in sp.sympify(a).atoms():
        # a is complex
        real_part_a = sp.re(a)
        imag_part_a = sp.im(a)
        if isinstance(b, complex) or sp.I in sp.sympify(b).atoms():
            # b is also complex
            real_part_b = sp.re(b)
            imag_part_b = sp.im(b)
            real = (real_part_a + real_part_b)
            imag = (imag_part_a + imag_part_b) 
            return complex(real, imag)
        else:
            # b is real
            real = b + real_part_a
            imag = imag_part_a
            return complex(real, imag)
    else:
        # a is real
        if isinstance(b, complex) or sp.I in sp.sympify(b).atoms():
            # b is complex
            real_part_b = sp.re(b)
            imag_part_b = sp.im(b)
            real = a + real_part_b
            imag = imag_part_b
            return complex(real, imag)
        else:
            # b is also real
            return a + b
        
def sub(a, b):
    """
    Subtract two numbers, handling complex and real values.

    Parameters:
    a (numeric or complex): The first number.
    b (numeric or complex): The second number.

    Returns:
    numeric or complex: The result of subtracting 'b' from 'a'.
    """
    if isinstance(a, complex) or sp.I in sp.sympify(a).atoms():
        # a is complex
        real_part_a = sp.re(a)
        imag_part_a = sp.im(a)
        if isinstance(b, complex) or sp.I in sp.sympify(b).atoms():
            # b is also complex
            real_part_b = sp.re(b)
            imag_part_b = sp.im(b)
            real = (real_part_a - real_part_b)
            imag = (imag_part_a - imag_part_b) 
            return complex(real, imag)
        else:
            # b is real
            real = real_part_a - b
            imag = imag_part_a
            return complex(real, imag)
    else:
        # a is real
        if isinstance(b, complex) or sp.I in sp.sympify(b).atoms():
            # b is complex
            real_part_b = sp.re(b)
            imag_part_b = sp.im(b)
            real = a - real_part_b
            imag = imag_part_b
            return complex(real, imag)
        else:
            # b is also real
            return a - b

def mul(a, b):
    """
    Multiply two numbers, handling complex and real values.

    Parameters:
    a (numeric or complex): The first number.
    b (numeric or complex): The second number.

    Returns:
    numeric or complex: The result of multiplying 'a' and 'b'.
    """
    if isinstance(a, complex) or sp.I in sp.sympify(a).atoms():
        real_part_a = sp.re(a)
        imag_part_a = sp.im(a)
        if isinstance(b, complex) or sp.I in sp.sympify(b).atoms():
            real_part_b = sp.re(b)
            imag_part_b = sp.im(b)
            real = (real_part_a * real_part_b) - (imag_part_a * imag_part_b)
            imag = (real_part_a * imag_part_b) + (imag_part_a * real_part_b)
            return complex(real, imag)
        else:
            real = b * real_part_a
            imag = b * imag_part_a
            return complex(real, imag)
    else:
        if isinstance(b, complex) or sp.I in sp.sympify(b).atoms():
            real_part_b = sp.re(b)
            imag_part_b = sp.im(b)
            real = a * real_part_b
            imag = a * imag_part_b
            return complex(real, imag)
        else:
            return a * b
        
def div(a, b):
    """
    Divide two numbers, handling complex and real values.

    Parameters:
    a (numeric or complex): The numerator.
    b (numeric or complex): The denominator.

    Returns:
    numeric or complex: The result of dividing 'a' by 'b'.
    """
    if isinstance(a, complex) or sp.I in sp.sympify(a).atoms():
        # a is complex
        real_part_a = sp.re(a)
        imag_part_a = sp.im(a)
        if isinstance(b, complex) or sp.I in sp.sympify(b).atoms():
            # b is also complex
            real_part_b = sp.re(b)
            imag_part_b = sp.im(b)
            denom = (real_part_b * real_part_b) + (imag_part_b * imag_part_b)
            real = ((real_part_a * real_part_b) + (imag_part_a * imag_part_b)) / denom
            imag = ((imag_part_a * real_part_b) - (real_part_a * imag_part_b)) / denom
            return complex(real, imag)
        else:
            # b is real
            real = real_part_a / b
            imag = imag_part_a / b
            return complex(real, imag)
    else:
        # a is real
        if isinstance(b, complex) or sp.I in sp.sympify(b).atoms():
            # b is complex
            real_part_b = sp.re(b)
            imag_part_b = sp.im(b)
            numerator = mul(a, b)
            denominator = mul(b, b)
            return div(numerator, denominator)
        else:
            # b is also real
            return a / b
        
def square_root(num):
    """
    Calculate the square root of a number, handling complex and real values.

    Parameters:
    num (numeric or complex): The number for which the square root is calculated.

    Returns:
    numeric or complex: The square root of 'num'.
    """
    # Verify if the number is complex
    if isinstance(num, complex) or sp.I in sp.sympify(num).atoms():
        # Calculate the square root of the complex number
        square_root = cm.sqrt(num)
        return square_root
    else:
        # If its not complex then
        square_root = sp.sqrt(float(num))
        return square_root