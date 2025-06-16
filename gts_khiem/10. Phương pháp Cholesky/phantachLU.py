import numpy as np
from tabulate import tabulate


def read_augmented_matrix_from_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        rows = len(lines)
        cols = len(lines[0].split())  # Trừ 1 cột cuối cùng là cột tự do
        matrix = np.zeros((rows, cols))
        for i in range(rows):
            line_values = [float(x) for x in lines[i].split()]
            matrix[i] = line_values
        return matrix
def extract_matrix_LU(A, m, n):
    L = np.zeros((m, n), dtype=float)
    U = np.zeros((m, n), dtype=float)
    for t in range(n):
        for i in range(t, m):
            sum_L = 0.0
            for j in range(t):
                sum_L += L[i, j] * U[j, t]
            L[i, t] = A[i, t] - sum_L
        for j in range(t, n):
            sum_U = 0.0
            for k in range(t):
                sum_U += L[t, k] * U[k, j]
            U[t, j] = (A[t, j] - sum_U) / L[t, t]
    return L, U
def main():
    filename = "input.txt"
    A = read_augmented_matrix_from_file(filename)
    m, n = A.shape
    L, U = extract_matrix_LU(A, m, n)
    print("Matrix A:")
    print(A)
    print("Matrix L:")
    print(L)
    print("Matrix U:")
    print(U)
    print(np.dot(L, U))
    # print(check_cross_dominant(A))
if __name__ == "__main__":
    main()