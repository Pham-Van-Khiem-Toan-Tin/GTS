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
def check_diagonal_dominance_bool(matrix):
    A = np.array(matrix)
    if A.shape[0] != A.shape[1]:
        raise ValueError("Ma trận không vuông, không thể kiểm tra chéo trội")

    n = A.shape[0]
    diag = np.abs(np.diag(A))

    # Chéo trội theo hàng
    row_sums = np.sum(np.abs(A), axis=1) - diag
    is_row_dominant = np.all(diag > row_sums)

    # Chéo trội theo cột
    col_sums = np.sum(np.abs(A), axis=0) - diag
    is_col_dominant = np.all(diag > col_sums)

    return is_row_dominant, is_col_dominant


def solve_by_jacobi_row(A, b, X_0, maxLoop, exp, TH):
    m = len(A)
    # Biến đổi A thành Bx + d
    for i in range(m):
        for j in range(m):
            if i != j:
                A[i, j] = -A[i,j] / A[i,i]    
        b[i] = b[i] / A[i,i]
        A[i, i] = 0
    augmented_matrix = np.concatenate((A, b.reshape(m,1)), axis = 1)
    print("Ma trận sau biến đổi thành Bx + d")
    print(tabulate(augmented_matrix, tablefmt="grid", floatfmt=".8f"))
    #Giải hệ
    B = augmented_matrix[:,:-1]
    d = augmented_matrix[:, -1].reshape(m, 1)
    row_sums = np.sum(np.abs(B), axis=1)
    for i in range(len(row_sums)):
        print(f"Tổng trị tuyệt đối hàng {i} là: {row_sums[i]}")
    q = norm_inf = np.linalg.norm(B, ord=np.inf)
    print(f"Chuẩn vô cùng của B là: {norm_inf}")
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
        headers = np.append(headers, "x_" + str(i + 1))
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
def solve_by_jacobi_column(A, b, X_0, maxLoop, exp, TH):
    m = len(A)
    #Biến đổi hệ thành Y = (I-AT)Y + b
    diag_element = np.diag(A)
    inv_diag = 1.0 / diag_element
    T = np.diag(inv_diag)
    print("Ma trận T:")
    print(tabulate(T, tablefmt="grid", floatfmt=".8f"))
    print(T)
    I = np.eye(m)
    B_I_AT = I - np.dot(A, T)
    print("Ma trận I - AT:")
    print(tabulate(B_I_AT, tablefmt="grid", floatfmt=".8f"))
    q = np.linalg.norm(B_I_AT, ord=1)
    print(f"Chuẩn tuyệt đối của ma trận I - AT là q: {q}")
    if q > 1:
        print("Bài toán không hội tụ")
        return
    B = I - np.dot(T, A)
    print("Ma trận B = I - TA: ")
    print(tabulate(B, tablefmt="grid", floatfmt=".8f"))
    d = np.dot(T, b)
    print("Ma trận d = Tb:")
    print(tabulate([d], tablefmt="grid", floatfmt=".8f"))
    lam_da = np.max(np.abs(np.diag(A))) / np.min(np.abs(np.diag(A)))
    print(f"Lamda: {lam_da}")
    k = 0
    check_exp = True
    X_pre = X_0
    X_k = X_0
    n = len(B[0])
    headers = np.array(("k"))
    for i in range(n):
        headers = np.append(headers, "x_" + str(i + 1))
    headers = np.append(headers, "exp")
    data = np.array([np.concatenate(([0], X_pre.reshape(-1), [None]))])
    print(tabulate(data, headers=headers, tablefmt="grid", floatfmt=".8f"))
    X_1 = np.dot(B, X_k) + d
    exp_X1_X0 = np.linalg.norm(X_1 - X_0, ord=1)
    if TH == 1:
        maxLoop = math.log(exp * (1-q)/exp_X1_X0/lam_da)/math.log(q)
        print(f"Số lần lặp tối đa theo công thức sai số tiên nghiệm: {maxLoop}")
    while k < maxLoop and check_exp:
        X_pre = X_k
        X_k = np.dot(B, X_k) + d
        exp_hn = exp*(1-q)/(lam_da * q)
        norm_Xk = np.linalg.norm(X_k - X_pre, ord=1)
        if norm_Xk <= exp_hn and TH == 2:
            check_exp = False
        
        data = np.array([np.concatenate(([int(k + 1)], X_k.reshape(-1),[norm_Xk] ))])
        print(tabulate(data, tablefmt="grid", floatfmt=".8f"))
        k += 1
    norm_Xk = np.linalg.norm(X_k - X_pre, ord=1)
    print(f"Sai số tuyệt đối delta_Xk: {norm_Xk}")
    exp_tgd = norm_Xk / np.linalg.norm(X_k, ord=1)
    print(f"Sai số tương đối: {round(exp_tgd * 100, 3)}%")
def main():
    filename = "jacobi.txt"
    X_0 = np.array([[1.0, 2.0, 1.0, -1.0]]).T
    exp = 0.01
    maxLoop = 10
    TH = 2 # 1 là tiên nghiệm; 2 là hậu nghiệm
    A,b = read_augmented_matrix_from_file(filename)
    row_dom, col_dom = check_diagonal_dominance_bool(A)
    if row_dom and col_dom:
        print("Ma trận chéo trội theo cả hàng và cột.")
        print("Bạn có thể chọn phương pháp giải:")
        print("1. Phương pháp Jacobi (theo hàng)")
        print("2. Phương pháp Jacobi (theo cột)")
        choice = input("Chọn phương pháp (1 hoặc 2): ")
        if choice == '1':
            solve_by_jacobi_row(A, b, X_0, maxLoop, exp, TH)
        elif choice == '2':
            solve_by_jacobi_column(A, b.reshape(len(A), -1), X_0, maxLoop, exp, TH)
        else:
            print("Lựa chọn không hợp lệ.")
    elif row_dom:
        print("Ma trận chéo trội theo hàng → Áp dụng phương pháp Jacobi theo hàng.")
        solve_by_jacobi_row(A, b, X_0, maxLoop, exp, TH)
    elif col_dom:
        print("Ma trận chéo trội theo cột → Áp dụng phương pháp Jacobi theo cột.")
        solve_by_jacobi_column(A, b.reshape(len(A), -1), X_0, maxLoop, exp, TH)

    else:
        print("⚠️ Ma trận không chéo trội. Không đảm bảo hội tụ khi dùng phương pháp lặp.")
if __name__ == "__main__":
    main()