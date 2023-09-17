from flask import Flask, render_template, request
from functions import *

# Folio Estatuto 200054087

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/method/bisection', methods=['GET', 'POST'])
def bisection_page():
    iterations = []
    resultado = None
    if request.method == 'POST':
        function = request.form['funcion']
        variable = request.form['variable']
        valor_a = float(request.form['valor_a'])
        valor_b = float(request.form['valor_b'])
        precision = float(request.form['precision'])
        max_iterations = float(request.form['max_iterations'])
        tolerancy = calculateTolerancy(precision)

        iterations, resultado = bisection(function, variable, valor_a, valor_b, tolerancy, max_iterations)

    return render_template('bisection.html', resultado=resultado, iterations=iterations)

@app.route('/method/newton', methods=['GET', 'POST'])
def newton_page():
    resultado = None
    if request.method == 'POST':
        function = request.form['funcion']
        variable = request.form['variable']
        x_value = float(request.form['x_value'])
        precision = float(request.form['precision'])
        tolerancy = calculateTolerancy(precision)
        max_iterations = float(request.form['max_iterations'])

        resultado = newton_raphson(function, variable, x_value, tolerancy, max_iterations)

    return render_template('newraph.html', resultado=resultado)

@app.route('/method/secant', methods=['GET', 'POST'])
def secant_page():
    resultado = None
    if request.method == 'POST':
        function = request.form['funcion']
        variable = request.form['variable']
        a_value = float(request.form['valor_a'])
        b_value = float(request.form['valor_b'])
        precision = float(request.form['precision'])
        max_iterations = float(request.form['max_iterations'])
        tolerancy = calculateTolerancy(precision)

        resultado = secant(function, variable, a_value, b_value, tolerancy, max_iterations)

    return render_template('secant.html', resultado=resultado)

@app.route('/method/muller', methods=['GET', 'POST'])
def muller_page():
    iterations = []
    resultado = None
    if request.method == 'POST':
        function = request.form['funcion']
        variable = request.form['variable']
        x0_value = float(request.form['x0_value'])
        x1_value = float(request.form['x1_value'])
        x2_value = float(request.form['x2_value'])
        precision = float(request.form['precision'])
        max_iterations = float(request.form['max_iterations'])
        tolerancy = calculateTolerancy(precision)

        iterations, resultado = muller(function, variable, x0_value, x1_value, x2_value, tolerancy, max_iterations)

    return render_template('muller.html', resultado=resultado, iterations=iterations)


if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0')