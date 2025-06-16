from decimal import Decimal, getcontext
from sympy import *
import sys
from tabulate import tabulate
import numpy as np
import math

def read_augmented_matrix_from_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        rows = len(lines)
        cols = len(lines[0].split()) - 1  # Trừ 1 cột cuối cùng là cột tự do
        matrix = np.zeros((rows, cols))
        b = np.zeros(rows)
        for i in range(rows):
            line_values = [float(x) for x in lines[i].split()]
            matrix[i] = line_values[:-1]
            b[i] = line_values[-1]
        return matrix, b
# chuyển đổi ma trận mở rộng về ma trận B của X = BX + d
def convert_to_BXd_matrix(A, b):
    m = len(A)
    for i in range(m):
        b[i] = b[i] / A[i,i]
        for j in range(len(A[i])):
            if i != j:
                A[i,j] = -A[i,j]/A[i,i]
        A[i,i] = 0
    augmented_matrix = np.concatenate((A, b.reshape(m,1)), axis = 1)
    print(tabulate(augmented_matrix, tablefmt="grid", floatfmt=".8f"))
    return augmented_matrix

def solve_single_loop_method(BXd, X_0, maxLoop, exp, TH):
    B = BXd[:,:-1]
    d = BXd[:, -1].reshape(3, 1)
    standard_array = np.zeros(len(B), dtype=float)
    for i in range(len(B)):
        sum_row = 0.0
        for j in range(len(B[i])):
            sum_row += abs(B[i,j])
        print(f"Tổng trị tuyệt đối hàng {i + 1}: {sum_row}")
        standard_array[i] = sum_row
    q = round(max(standard_array), 8)
    print(f"q = {q}")
    if q > 1:
        print("Bài toán không hội tụ")
        return
    k = 0
    check_exp = True
    X_pre = X_0
    X_k = X_0
    n = len(B[0])
    headers = np.array(("k"))
    for i in range(n):
        headers = np.append(headers, "x_" + str(i))
    headers = np.append(headers, "exp")
    data = np.array([np.concatenate(([0], X_pre.reshape(-1), [None]))])
    print(tabulate(data, headers=headers, tablefmt="grid", floatfmt=".8f"))
    if TH == 1:
        X_1 = np.dot(B, X_k) + d
        maxLoop = math.log(exp * (1-q)/np.max(np.abs(X_1 - X_k)))/math.log(q)
        print(f"Số lần lặp tối đa theo công thức sai số tiên nghiệm: {maxLoop}")
    while k < maxLoop and check_exp:
        X_pre = X_k
        X_k = np.dot(B, X_k) + d
        delta_X = q * np.max(np.abs(X_k - X_pre)) / (1-q)
        omega_X_k = 0.0
        if TH == 2:
            omega_X_k = delta_X * 100 / np.max(np.abs(X_k))
            if exp is not None and omega_X_k < exp:
                check_exp = False
        if TH == 4:
            omega_X_k = np.max(np.abs(X_k - X_pre))
            if(omega_X_k < exp / 2):
                check_exp = False
        if TH == 3:
            omega_X_k = q*np.max(np.abs(X_k - X_pre))/(1-q)
            omega_X_k = omega_X_k/np.max(np.abs(X_k)) 
        data = np.array([np.concatenate(([int(k + 1)], X_k.reshape(-1), [omega_X_k] ))])
        print(tabulate(data, tablefmt="grid", floatfmt=".8f"))
        k += 1
    
def main():
    filename = "input.txt"
    A,b = read_augmented_matrix_from_file(filename)
    augmented_matrix = convert_to_BXd_matrix(A, b)
    X_0 = np.array([[78.0, 87.0, 78.0]]).T
    exp = 1
    maxLoop = 5
    TH = 3 # 1 là tiên nghiệm; 2 là hậu nghiệm, 3 là sai tương đối, 4 là sai số tuyệt đối
    solve_single_loop_method(augmented_matrix, X_0, maxLoop, exp, TH)
if __name__ == "__main__":
    main()