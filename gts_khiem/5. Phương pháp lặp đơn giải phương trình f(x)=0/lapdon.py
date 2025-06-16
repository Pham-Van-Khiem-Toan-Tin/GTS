from decimal import Decimal, getcontext
from sympy import *
from sympy.calculus.util import continuous_domain
import sys
from tabulate import tabulate
import numpy as np
#kiểm tra liên tục
def is_continuous_on_interval(func, a, b):
    x = Symbol("x")
    # Tìm miền liên tục trên toàn R
    cont_domain = continuous_domain(func, x, S.Reals)
    # Định nghĩa đoạn [a, b]
    interval = Interval(a, b)
    # Kiểm tra xem [a, b] có nằm hoàn toàn trong miền liên tục không
    if interval.is_subset(cont_domain):
        return True, f"Hàm liên tục trên [{a}, {b}]"
    else:
        # Tìm các điểm không thuộc miền liên tục
        discontinuities = interval - cont_domain
        return False, f"Hàm không liên tục tại các điểm: {discontinuities}"
def check_sign_change(f, a, b):
    # Kiểm tra sự thay đổi dấu giữa hai điểm a và b
    x = symbols("x")
    fa = f.subs(x, a).evalf()
    fb = f.subs(x, b).evalf()
    return fa * fb < 0
def single_loop_method(func, funcx, a, b,x_0, exp, maxLoop):
    x = symbols("x")
    # chứng minh ánh xạ compact
    corners = [(a) for a in [a, b]]
    values = [
        funcx.subs({x: a}).evalf() for (a) in corners
    ]
    max_value = max(values)
    max_index = values.index(max_value)  # Lấy chỉ số của giá trị max trong values
    max_point = corners[max_index]
    print(f"Max = {max_value} khi x = {max_point}")

    min_value = min(values)
    min_index = values.index(min_value)  # Lấy chỉ số của giá trị min trong values
    min_point = corners[min_index]
    print(f"Min = {min_value} khi x = {min_point}")
    function_string = funcx.__str__()
    if min_value >= a and max_value <= b:
        print(f"Hàm số {function_string} có thể ánh xạ từ [{a}, {b}] sang [ {a}, {b}]")
    else:
        print(f"Hàm số {function_string} không thể ánh xạ từ [ {a}, {b}] sang [ {a}, {b}]")
        return
    func_prime = Abs(diff(funcx, x))
    corners_prime = [(a) for a in [a, b]]
    values_x_prime = [
        func_prime.subs({x: a}).evalf() for (a) in corners_prime
    ]
    max_value_x_prime = max(values_x_prime)
    max_index_x_prime = values_x_prime.index(
        max_value_x_prime
    )  # Lấy chỉ số của giá trị max trong values
    max_point_x_prime = corners_prime[max_index_x_prime]
    print(f"Max trị tuyệt đối đạo hàm theo x = {max_value_x_prime} khi (x) = {max_point_x_prime}")
    if(max_value_x_prime > 1):
        print(f"Không thể sử dụng lặp đơn do kmax = {max_value_x_prime} > 1")
        return 
    exp_0 = exp * (1- max_value_x_prime) / max_value_x_prime # hậu nghiệm
    # exp_0 = exp * (1 - max_value_x_prime) # tiên nghiệm
    print(f"Sai số ban đầu: {exp_0}")
    x_pre = x_0
    xk = x_0
    iteration = 1
    check_exp = True
    headers = np.array(["k", "xk", "||xk - xk-1||extremely"])
    data = np.array([[0, str(x_pre), None]])
    # x_1 = funcx.subs(x, x_0) # tiên nghiệm
    print(tabulate(data, headers=headers, tablefmt="grid", floatfmt=".8f"))
    while check_exp and iteration < maxLoop:
        x_pre = xk
        xk = funcx.subs(x, x_pre)
        # exp_0 = exp_0 / max_value_x_prime**(iteration) # tiên nghiệm
        # if(abs(x_1 - x_0) < exp_0): # tiên nghiệm
        #     check_exp = False # tiên nghiệm
        if(abs(xk - x_pre) < exp_0): # hậu nghiệm
            check_exp = False # hậu nghiệm
        data = np.array([[iteration, str(xk), str(abs(xk - x_pre))]])
        print(tabulate(data, headers=headers, tablefmt="grid", floatfmt=".8f"))
        iteration += 1
    print("Giá trị xk cuối cùng là: ", xk)
def main():
    # Đặt số chữ số đáng tin
    # n = int(input("Nhập số chữ số đáng tin (n): "))
    # getcontext().prec = n + 2  # Số chữ số đáng tin + 2 (tăng thêm 2 chữ số để đảm bảo độ chính xác)
    x = sympify("x")
    pi_eval = pi.evalf() # số pi
    # Đọc dữ liệu từ file text
    expr = x + 3/(5*x**2 -20) # Nhập hàm số ở đây
    func = sympify(expr)
    exprx = -3/(5*x**2 -20) # hàm phi x ví dụ x = ax + b; x = g(x)
    funcx = sympify(exprx) 
    
    # Trích xuất hệ số từ dữ liệu
    # coefficients = [Decimal(coeff) for coeff in data]

    # Trích xuất giá trị khoảng cách ly từ dữ liệu
    a = Decimal(input("Nhập giá trị a với dạng số thực: "))
    b = Decimal(input("Nhập giá trị b với dạng số thực: "))

    # Kiểm tra tính đơn điệu của hàm số
    check_continuous, continuous_message = is_continuous_on_interval(func, a, b)
    print(continuous_message)
    if(not check_continuous):
        return
    # Kiểm tra khoảng cách ly
    # if not check_sign_change(func, a, b):
    #     print("Khoảng cách ly không đúng.")
    #     return

    # Nhập sai số từ bàn phím
    epsilon = Decimal(input("Nhập sai số epsilon: "))
    x_o = float(input("Nhập giá trị x_0: ")) 
    # Tìm nghiệm bằng phương pháp lặp đơn
    single_loop_method(func, funcx, a, b, x_o, epsilon, maxLoop=100)


# Chạy chương trình
if __name__ == "__main__":
    main()