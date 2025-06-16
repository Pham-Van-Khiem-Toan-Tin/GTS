from decimal import Decimal, getcontext
from sympy import *
import sys
from tabulate import tabulate
import numpy as np

def check_sign_change(f, a, b):
    # Kiểm tra sự thay đổi dấu giữa hai điểm a và b
    x = symbols("x")
    fa = f.subs(x, a).evalf()
    fb = f.subs(x, b).evalf()
    return fa * fb < 0

def check_monotonicity(f, a_0, b_0):
    # Kiểm tra tính đơn điệu của hàm số trên đoạn [a, b]
    x = symbols("x")
    a = a_0
    b = b_0
    sol_set = solveset(diff(f, x ), x, Interval(a, b))
    sol_set = Union(sol_set, solveset(diff(f, x, 2), x, Interval(a, b)));
    if(sol_set.is_empty): 
        return True
    else: 
        return False



def bisection_method(f, a, b, epsilon, maxLoop):
    x = symbols("x")
    x_0 = d = 0
    fa_value = f.subs(x, a)
    fa_prime = diff(f, x)
    fa_double_prime = diff(fa_prime, x)
    fa_double_prime_value = fa_double_prime.subs(x, a)
    if fa_value * fa_double_prime_value < 0:
        x_0 = a
        d = b
    else:
        x_0 = b
        d = a
    min1 = MAX1 = fa_prime.subs(x, a)
    sol_set = FiniteSet(a, b)
    sol_set = Union(sol_set, solveset(diff(f, x), x, Interval(a,b)))
    if(sol_set.is_iterable): 
            for args in sol_set:
                f_prime_value = fa_prime.subs(x, args)
                abs_f_prime_value = Abs(f_prime_value)
                MAX1 = max(MAX1, abs_f_prime_value);
                min1 = min(min1, abs_f_prime_value);
    exp = min1*epsilon / (MAX1 - min1)
    iterationLoop = 1
    dx = -f.subs(x, x_0)*(x_0 - d)/(f.subs(x, x_0) - f.subs(x, d))
    data0 = np.array([[0, x_0, dx]])
    headers = np.array(["k", "xk", "dxk"])
    print("Lần lặp thứ 0")
    print(tabulate(data0, headers=headers, tablefmt="grid", floatfmt=".15f"))
    xk = x_0 + dx
    x_array = np.array([x_0, xk])
    while abs(dx) >= exp and iterationLoop < maxLoop:
        print("Lần lặp thứ ", iterationLoop)
        dx = -f.subs(x, xk)*(xk - d)/(f.subs(x, xk) - f.subs(x, d))
        data = np.array([[iterationLoop, xk, dx]])
        headers = np.array(["k", "xk", "dxk"])
        print(tabulate(data, headers=headers, tablefmt="grid", floatfmt=".15f"))
        xk = xk + dx
        x_array = np.append(x_array, xk)
        iterationLoop += 1
    print(x_array[iterationLoop - 1], x_array[iterationLoop - 2], MAX1, min1)
    expAnalyst = (MAX1 - min1)*abs(x_array[iterationLoop - 1] - x_array[iterationLoop - 2])/min1
    print("Đánh giá sai số (M1 - m1)|x", iterationLoop-1, " - x", iterationLoop - 2, "|/m1 =",  expAnalyst)
    return x_array[iterationLoop - 1]
def main():
    x = simplify("x")
    expr = 1.2*x**5 - 2.57*x + 2
    func = sympify(expr)
    a = float(input("Nhập giá trị a: "))
    b = float(input("Nhập giá trị b: "))
    if not check_sign_change(func, a, b):
        print("Khoảng cách ly không đúng.")
        return
    if not check_monotonicity(func, a, b):
        print("Hàm số không đơn điệu trên khoảng cách ly.")
        return
    epsilon = Decimal(input("Nhập sai số epsilon: "))
    maxLoop = Integer(input("Nhập số lần lặp tối đa: "))
    result = bisection_method(func, a, b, epsilon, maxLoop)
    if result is not None:
        print("Nghiệm của phương trình: ", result)
if __name__ == "__main__":
    main()
