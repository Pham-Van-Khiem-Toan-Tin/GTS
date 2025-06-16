import numpy as np
from tabulate import tabulate

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
def check_diagonal_dominance(A):
    m, n = A.shape
    for i in range(m):
        sum_row = 0
        for j in range(n):
            if i != j:
                sum_row += abs(A[i, j])
        if abs(A[i, i]) < sum_row:
            return False
    return True
def find_matrix_B_and_d(A, b):
    m, n = A.shape
    b_matrix = np.copy(b.reshape(m, 1))
    B = np.copy(A)
    for i in range(m):
        b_matrix[i,0] = b_matrix[i,0] / A[i,i]
        for j in range(n):
            if i != j:
                B[i, j] = -B[i, j] / A[i, i]
        B[i, i] = 0
    return B, b_matrix
def find_matrix_L_and_U(A, m, n):
    L = np.zeros_like(A, dtype=float)
    U = np.zeros_like(A, dtype=float)
    for i in range(m):
        for j in range(n):
            if(j < i):
                L[i, j] = A[i, j]
            elif(j > i):
                U[i,j] = A[i,j]
    return L, U
def find_standard_array(B):
    standard_array = np.zeros(len(B), dtype=float)
    for i in range(len(B)):
        sum_L = 0.0
        sum_U = 0.0
        for j in range(len(B[i])):
            if j < i:
                sum_L += abs(B[i, j])
            elif j > i:
                sum_U += abs(B[i, j])
        standard_array[i] = sum_L / (1 - sum_U)
    print(standard_array)
    return max(standard_array)
def gauss_seidel(L, U, d, X_0, m, standard, maxLoop, exp):
    X_k = X_0
    X_pre = X_0
    k = 0
    check_exp = True
    I = np.eye(L.shape[0])
    I_L_inv = np.linalg.inv(I - L)
    while k < maxLoop and check_exp:
        X_pre = X_k.copy()
        X_k = np.dot(I_L_inv, np.dot(U, X_pre)) + np.dot(I_L_inv, d)
        data = np.array([np.concatenate(([int(k + 1)], X_k.reshape(-1)))])
        print(tabulate(data, tablefmt="grid", floatfmt=".8f"))
        delta_X_matrix =abs(abs(X_k) - abs(X_pre))
        if standard*max(delta_X_matrix)/(1 - standard) < exp:
            check_exp = False
        k += 1
    print("Giá trị nghiệm x:")
    print(X_k)
def main():
    X_0 = np.array([1.0, 1.0, 1.0, 1.0]) # Giá trị khởi tạo cho x
    filename = "input.txt"  # Tên tệp chứa ma trận mở rộng
    A, b = read_augmented_matrix_from_file(filename)
    m, n = A.shape
    print("Ma trận mở rộng:")
    print(tabulate(np.concatenate((A, b.reshape(m, 1)), axis=1), tablefmt="grid", floatfmt=".8f"))
    B, d = find_matrix_B_and_d(A, b)
    print("Ma trận B:")
    print(tabulate(B, tablefmt="grid", floatfmt=".8f"))
    print("Ma trận d:")
    print(tabulate(d, tablefmt="grid", floatfmt=".8f"))
    L, U = find_matrix_L_and_U(B, m, n)
    print("Ma trận L:")
    print(tabulate(L, tablefmt="grid", floatfmt=".8f"))
    print("Ma trận U:")
    print(tabulate(U, tablefmt="grid", floatfmt=".8f"))
    standard = find_standard_array(B)
    print("Hệ số co:")
    print(standard)
    # Số lần lặp tối đa
    maxLoop = int(input("Nhập số lần lặp tối đa: "))
    #Sai số
    exp = float(input("Nhập sai số: "))
    # Giải hệ phương trình Ax = b bằng phương pháp Gauss
    gauss_seidel(L, U, d, X_0.reshape(n, 1), m, standard, maxLoop, exp)

if __name__ == "__main__":
    main()