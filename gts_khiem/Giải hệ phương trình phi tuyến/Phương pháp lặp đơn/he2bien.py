from decimal import Decimal, getcontext
from sympy import *
import sys
from tabulate import tabulate
import numpy as np


def check_target_mapping(fs, D):
    x1, x2 = symbols("x1 x2")
    result = np.array([False, False])
    for i in range(len(fs)):
        # Tính giá trị hàm tại các góc
        corners = [(a, b) for a in D[0] for b in D[1]]
        values = [fs[i].subs({x1: a, x2: b}).evalf() for (a, b) in corners]
        max_value = max(values)
        max_index = values.index(max_value)  # Lấy chỉ số của giá trị max trong values
        max_point = corners[max_index]
        print(f"Max = {max_value} khi (x1, x2) = {max_point}")

        min_value = min(values)
        min_index = values.index(min_value)  # Lấy chỉ số của giá trị min trong values
        min_point = corners[min_index]
        print(f"Min = {min_value} khi (x1, x2) = {min_point}")
        function_string = fs[i].__str__()

        if min_value >= D[i][0] and max_value <= D[i][1]:
            print("Hàm số", function_string, " có thể ánh xạ từ D sang D")
            result[i] = True
        else:
            print("Hàm số", function_string, " không thể ánh xạ từ D sang D")
    if all(result):
        return True
    else:
        return False


def find_maximum_prime_row(f, D):
    x1, x2 = symbols("x1 x2")
    result_sum = 0.0 

    f_prime_x1 = diff(f, x1)
    corners_x1 = [(a, b) for a in D[0] for b in D[1]]
    values_x1 = [f_prime_x1.subs({x1: a, x2: b}).evalf() for (a, b) in corners_x1]
    max_value_x1 = max(values_x1)
    max_index_x1 = values_x1.index(
        max_value_x1
    )  # Lấy chỉ số của giá trị max trong values
    max_point_x1 = corners_x1[max_index_x1]
    result_sum += max_value_x1
    print(f"Max trị tuyệt đối đạo hàm theo x1 = {max_value_x1} khi (x1, x2, x3) = {max_point_x1}")
    f_prime_x2 = diff(f, x2)
    corners_x2 = [(a, b) for a in D[0] for b in D[1]]
    values_x2 = [f_prime_x2.subs({x1: a, x2: b}).evalf() for (a, b) in corners_x2]
    max_value_x2 = max(values_x2)
    max_index_x2 = values_x2.index(
        max_value_x2
    )  # Lấy chỉ số của giá trị max trong values
    max_point_x2 = corners_x2[max_index_x2]
    result_sum += max_value_x2
    print(f"Max trị tuyệt đối đạo hàm theo x2 = {max_value_x2} khi (x1, x2, x3) = {max_point_x2}")
    return result_sum


def check_standard(fs, D):
    standards = np.array([0.0, 0.0])
    for i in range(len(fs)):
        max_value_i = find_maximum_prime_row(fs[i], D)
        standards[i] = max_value_i
    max_standards = max(standards)
    if max_standards < 1:
        return max_standards, True
    else:
        return max_standards, False


def main():
    x1, x2 = symbols("x1 x2")
    f = sin(x1*x2)/3
    g = cos(x1**2 + x2**2)/4
    D = [(0.0, 0.5), (0.0, 0.5)]
    X_0 = [0.0, 0.25]  # Nhập X_0 ở đây
    result_check_mapping = check_target_mapping([f, g], D)
    if result_check_mapping:
        print("Hệ phương trình phi tuyến có thể ánh xạ từ D sang D")
    else:
        print("Hệ phương trình phi tuyến không thể ánh xạ từ D sang D")
        return
    max_standards, result_check_standard = check_standard([f, g], D)
    if result_check_standard:
        print("Hệ phương trình phi tuyến hội tụ")
    else:
        print("Hệ phương trình phi tuyến không hội tụ")
        return
    print(f"k= {max_standards}")
    exp = Float(input("Nhập giá trị sai số:"))
    exp_0 = (exp * (1 - max_standards) / max_standards)
    # exp_0 = float(exp * (1 - max_standards) / max_standards**0) # Sai số tiền nghiệm
    print("Sai số ban đầu là:", exp_0)
    iteration = 1
    x_pre = X_0
    xk=X_0
    headers = np.array(["k", "xk", "||xk - xk-1||extremely"])
    data = np.array([[0, str(x_pre), None]])
    print(tabulate(data, headers=headers, tablefmt="grid", floatfmt=".8f"))
    check_exp = True
    check_list = np.array([0.0, 0.0])
    
    while check_exp and iteration < 100:
        x_pre = xk
        xk = np.array(
            [
                f.subs({x1: xk[0], x2: xk[1]}),
                g.subs({x1: xk[0], x2: xk[1]}),
            ]
        )
        for i in range(len(x_pre)):
            check_data = xk[i] - x_pre[i]
            check_list[i] = abs(check_data)
        if check_list[0] < exp_0 and check_list[1] < exp_0:
            check_exp = False
        # exp_0 = float(exp * (1 - max_standards) / max_standards**iteration) # Sai số tiền nghiệm
        data = np.array([[iteration, str(xk), str(check_list)]])
        print(tabulate(data, headers=headers, tablefmt="grid", floatfmt=".8f"))
        iteration += 1
    print("Giá trị xk cuối cùng là: ", xk)


if __name__ == "__main__":
    main()
