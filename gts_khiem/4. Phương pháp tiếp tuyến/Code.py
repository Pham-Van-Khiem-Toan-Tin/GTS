from decimal import Decimal, getcontext
from sympy import *
import sys
from tabulate import tabulate
import numpy as np
# Kiểm tra liên tục và đơn điệu
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
#Kiểm tra fa * fb < 0
def check_sign_change(f, a, b):
    # Kiểm tra sự thay đổi dấu giữa hai điểm a và b
    x = symbols("x")
    fa = f.subs(x, a).evalf()
    fb = f.subs(x, b).evalf()
    return fa * fb < 0

def secant_method(f, a, b, epsilon, max_iterations):
    # Áp dụng phương pháp tiếp tuyến để tìm nghiệm
    x = symbols("x")
    x0 = a
    x1 = b
    iteration = 0
    
    while abs(x1 - x0) > epsilon and iteration < max_iterations:
        print("Lần lặp thứ ", iteration)
        data = np.array([[iteration, x1, abs(x1 - x0)]])
        headers = np.array(["k", "xk", "|xk - xk-1|"])
        print(tabulate(data, headers=headers, tablefmt="grid", floatfmt=".8f"))
        if f.subs(x,x0) == f.subs(x,x1):
            return None
        f_prime = diff(f, x)
        x_next = x1 - f.subs(x, x1) / f_prime.subs(x, x1)
        x0 = x1
        x1 = x_next
        iteration += 1

    if abs(x1 - x0) <= epsilon:
        return x1
    else:
        return None

# Hàm chính
def main():
    # Đặt số chữ số đáng tin
    # n = int(input("Nhập số chữ số đáng tin (n): "))
    # getcontext().prec = n + 2  # Số chữ số đáng tin + 2 (tăng thêm 2 chữ số để đảm bảo độ chính xác)
    x = sympify("x")
    # Đọc dữ liệu từ file text
    expr = x**5 - log(x)/log(E) - 12 # Nhập hàm tại đây
    func = sympify(expr)

    # Trích xuất hệ số từ dữ liệu
    # coefficients = [Decimal(coeff) for coeff in data]

    # Trích xuất giá trị khoảng cách ly từ dữ liệu
    a = Decimal(input("Nhập giá trị a: "))
    b = Decimal(input("Nhập giá trị b: "))

    # Kiểm tra tính đơn điệu của hàm số
    if not check_monotonicity(func, a, b):
        print("Hàm số không đơn điệu trên khoảng cách ly.")
        return

    # Kiểm tra khoảng cách ly
    if not check_sign_change(func, a, b):
        print("Khoảng cách ly không đúng.")
        return

    # Nhập sai số từ bàn phím
    epsilon = Decimal(input("Nhập sai số epsilon: "))

    # Tìm nghiệm bằng phương pháp tiếp tuyến
    result = secant_method(func, a, b, epsilon, max_iterations=100)

    if result is not None:
        print("Nghiệm của phương trình: ", result)

# Chạy chương trình
if __name__ == "__main__":
    main()
