import numpy as np
from tabulate import tabulate
import math

# kiểm tra ma trận đối xứng không suy biến
def is_symmetric_nonsingular(A, tol=1e-10):
    check_symmetric = True
    # Kiểm tra ma trận đối xứng
    if not np.allclose(A, A.T, atol=tol):
        check_symmetric = False
    # Kiểm tra ma trận không suy biến
    det_A = np.linalg.det(A)
    check_nonsingular = True
    if abs(det_A) < tol:
        check_nonsingular = False
    return check_symmetric, check_nonsingular

def read_matrix_A_from_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        rows = len(lines)
        cols = len(lines[0].split()) # Trừ 1 cột cuối cùng là cột tự do
        matrix = np.zeros((rows, cols))
        for i in range(rows):
            line_values = [float(x) for x in lines[i].split()]
            matrix[i] = line_values
        return matrix

def format_complex(c, precision=8):
    real_part = round(c.real, precision)
    imag_part = round(c.imag, precision)
    # Loại bỏ phần thập phân .0 nếu là số nguyên
    real_str = f"{real_part:.{precision}f}".rstrip("0").rstrip(".")
    imag_str = f"{imag_part:.{precision}f}".rstrip("0").rstrip(".")
    # Xử lý dấu của phần ảo
    if imag_part >= 0:
        return f"{real_str}+{imag_str}j"
    else:
        return f"{real_str}{imag_str}j"

def extrac_matrix_U(A, m, n):
    U = np.zeros((m, n), dtype=complex)
    for i in range(m):
        sum_row = 0.0
        for k in range(i):
            sum_row += U[k, i]**2
        U[i, i] = np.sqrt(A[i, i] - sum_row)
        # Đảm bảo phần ảo của U[i, i] là 0 nếu nó quá nhỏ
        if np.abs(U[i, i].imag) < 1e-10:
            U[i, i] = U[i, i].real
        for j in range(n):
            if(i != j):
                sum = 0.0
                for k in range(i):
                    sum += U[k, i] * U[k, j]
                U[i, j] = (A[i, j] - sum) / U[i, i]
                if np.abs(U[i, j].imag) < 1e-10:
                    U[i, j] = U[i, j].real
    return U

def gauss_of_St_method(matrix, m, n):
    augment_matrix = np.copy(matrix)
    result_y = np.zeros(n, dtype=complex)
    for i in range(m):
        sum_ = augment_matrix[i, -1]
        if i > 0:
            for j in range(i):
                sum_ -= augment_matrix[i, j] * result_y[j]
        result_y[i] = sum_ / augment_matrix[i,i]
    return result_y
def gauss_of_S_method(matrix, m,n):
    augment_matrix = np.copy(matrix)
    result = np.zeros(m, dtype=complex)
    for i in range(m - 1, -1, -1):
        sum_ = augment_matrix[i, -1]
        for j in range(i+1, n):
            sum_ -= augment_matrix[i, j] * result[j]
        result[i] = sum_ / augment_matrix[i, i]
        print(result[i])
    return result

def main():
    filename = "input.txt"
    A= read_matrix_A_from_file(filename)
    print("Ma trận A:")
    print(tabulate(A, tablefmt="grid"))
    check_symmetric, check_nonsingular = is_symmetric_nonsingular(A)
    print(f"Ma trận{'' if check_symmetric else ' không'} đối xứng")
    if not check_symmetric:
        print("Phải dùng At*A*x = At*B")
    print(f"Ma trận{' không' if check_nonsingular else ''} suy biến")
    m, n = A.shape
    U = extrac_matrix_U(A, m, n) if check_symmetric else A
    print("Ma trận U:")
    U_formatted = [[format_complex(U[i, j]) for j in range(U.shape[1])] for i in range(U.shape[0])]
    print(tabulate(U_formatted, tablefmt="grid"))
    print("Ma trận U^t:")
    U_T = U.T if check_symmetric else A.T
    U_T_formatted = [[format_complex(U_T[i, j]) for j in range(U_T.shape[1])] for i in range(U_T.shape[0])]
    print(tabulate(U_T_formatted, tablefmt="grid"))
    E_matrix = np.eye(m)
    result = None
    for j in range(n):
        b = E_matrix[:, j].reshape(-1, 1) if check_symmetric else np.dot(A.T, E_matrix[:, j].reshape(-1, 1))
        augment_ST_matrix = np.concatenate((U_T, b.astype(complex)), axis=1)
        augment_ST_matrix_formatted = [[format_complex(augment_ST_matrix[i, j]) for j in range(augment_ST_matrix.shape[1])] for i in range(augment_ST_matrix.shape[0])]
        print("Ma trận mở rộng của St*y = b:")
        print(tabulate(augment_ST_matrix_formatted, tablefmt="grid"))
        y = gauss_of_St_method(augment_ST_matrix, m, n)
        y_reshape = y.reshape(-1,1)
        y_formatted = [[format_complex(y_reshape[i, j]) for j in range(y_reshape.shape[1])] for i in range(y_reshape.shape[0])]
        print("Kết quả y:")
        print(tabulate(y_formatted, tablefmt="grid"))
        augment_S_matrix = np.concatenate((U, y_reshape), axis=1)
        print("Ma trận mở rộng của S*x = y:")
        augment_S_matrix_formatted = [[format_complex(augment_S_matrix[i, j]) for j in range(augment_S_matrix.shape[1])] for i in range(augment_S_matrix.shape[0])]
        print(tabulate(augment_S_matrix_formatted, tablefmt="grid"))
        result_complex_item = gauss_of_S_method(augment_S_matrix, m,n)
        if result is None:
            result = result_complex_item.reshape(-1,1)
        else:
            result = np.hstack((result, result_complex_item.reshape(-1,1)))
    print("Ma trận nghịch đảo là:")
    print(tabulate(result.astype(float), tablefmt="grid", floatfmt=".8f"))
if __name__ == "__main__":
    main()