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
    alert = None
    resultado = None
    if request.method == 'POST':
        function = request.form['funcion']
        variable = request.form['variable']
        valor_a = float(request.form['valor_a'])
        valor_b = float(request.form['valor_b'])
        precision = float(request.form['precision'])
        max_iterations = float(request.form['max_iterations'])
        tolerancy = calculateTolerancy(precision)

        iterations, resultado, alert = bisection(function, variable, valor_a, valor_b, tolerancy, max_iterations)

    return render_template('bisection.html', resultado=resultado, iterations=iterations, alert=alert)

@app.route('/method/newton', methods=['GET', 'POST'])
def newton_page():
    iterations = []
    alert = None
    resultado = None
    if request.method == 'POST':
        function = request.form['funcion']
        variable = request.form['variable']
        x_value = float(request.form['x_value'])
        precision = float(request.form['precision'])
        tolerancy = calculateTolerancy(precision)
        max_iterations = float(request.form['max_iterations'])

        iterations, resultado, alert = newton_raphson(function, variable, x_value, tolerancy, max_iterations)

    return render_template('newraph.html', resultado=resultado, iterations=iterations, alert=alert)

@app.route('/method/secant', methods=['GET', 'POST'])
def secant_page():
    iterations = []
    alert = None
    resultado = None
    if request.method == 'POST':
        function = request.form['funcion']
        variable = request.form['variable']
        a_value = float(request.form['valor_a'])
        b_value = float(request.form['valor_b'])
        precision = float(request.form['precision'])
        max_iterations = float(request.form['max_iterations'])
        tolerancy = calculateTolerancy(precision)

        iterations, resultado, alert = secant(function, variable, a_value, b_value, tolerancy, max_iterations)

    return render_template('secant.html', resultado=resultado, iterations=iterations, alert=alert)

@app.route('/method/fixed-point', methods=['GET', 'POST'])
def fixed_point_page():
    iterations = []
    alert = None
    resultado = None
    if request.method == 'POST':
        function = request.form['funcion']
        variable = request.form['variable']
        p0 = float(request.form['p0_value'])
        precision = float(request.form['precision'])
        max_iterations = float(request.form['max_iterations'])
        tolerancy = calculateTolerancy(precision)

        iterations, resultado, alert = fixed_point(function, variable, p0, tolerancy, max_iterations)

    return render_template('fixed-point.html', resultado=resultado, iterations=iterations, alert=alert)

@app.route('/method/muller', methods=['GET', 'POST'])
def muller_page():
    iterations = []
    alert = None
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

        iterations, resultado, alert = muller(function, variable, x0_value, x1_value, x2_value, tolerancy, max_iterations)

    return render_template('muller.html', resultado=resultado, iterations=iterations, alert=alert)

@app.errorhandler(404)
def page_not_found(error):
    return render_template('404.html'), 404

if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0')