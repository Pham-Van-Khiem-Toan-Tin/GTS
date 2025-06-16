from decimal import Decimal, getcontext
from sympy import *
from tabulate import tabulate
import numpy as np
def calculate_jacobian_xk(jacobi, xk):
    x1, x2, x3 = symbols("x1 x2 x3")
    jacobi_xk = np.zeros((3,3), dtype=float)
    for i in range(3):
        for j in range(3):
            jacobi_xk[i][j] = jacobi[i][j].subs({x1: xk[0], x2: xk[1], x3: xk[2]})
    return jacobi_xk
def newton_method(fs, jacobi, x0, maxLoop):
    x1, x2, x3 = symbols("x1 x2 x3")
    iterationLoop = 0
    print("Kết quả phương pháp Newton")
    headers = np.array(["k", "x1", "x2", "x3"])
    data = np.array([[Integer(iterationLoop), x0[0], x0[1], x0[2]]])
    print(tabulate(data, headers=headers, tablefmt="grid", floatfmt=".12f"))
    xk = np.array(x0)
    while iterationLoop < maxLoop:
        fk = fs[0].subs({x1: xk[0], x2: xk[1], x3: xk[2]})
        gk = fs[1].subs({x1: xk[0], x2: xk[1], x3: xk[2]})
        hk = fs[2].subs({x1: xk[0], x2: xk[1], x3: xk[2]})
        Fxk = np.array([fk, gk, hk])
        jacobi_xk = calculate_jacobian_xk(jacobi, xk)
        jacobi_xk_inv = np.linalg.inv(jacobi_xk)
        h = np.dot(jacobi_xk_inv, Fxk)
        xk = xk - h
        iterationLoop += 1
        data = np.array([[ iterationLoop, xk[0], xk[1], xk[2]]])
        print(tabulate(data, headers=headers, tablefmt="grid", floatfmt=".8f"))
    return xk
def main():
    pi_numrical = pi.evalf() # số pi
    x1, x2, x3 = symbols("x1 x2 x3")
    f = 3*x1 - cos(x2 * x3) - 0.5
    g = 4*x1**2 - 625*x2**2 + 2*x2 -1
    h = E**(-x1*x2) + 20*x3 + (10*pi_numrical - 3)/3
    X_0 = [0.0, 0.0, 0.0]
    matrix_j = np.array([
        [diff(f, x1), diff(f, x2), diff(f, x3)],
        [diff(g, x1), diff(g, x2), diff(g, x3)],
        [diff(h, x1), diff(h, x2), diff(h, x3)]
    ])
    print("Ma trận Jacobian:")
    print(tabulate(matrix_j, tablefmt="grid", floatfmt=".8f"))
    newton_method([f, g, h], matrix_j, X_0, 5)
if __name__ == "__main__":
    main()