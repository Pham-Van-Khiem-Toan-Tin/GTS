from fractions import Fraction

import numpy as np
from tabulate import tabulate

DECIMALS = 7

def read_input(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    n = int(lines[0])
    A = np.array([list(map(float, lines[i + 1].split())) for i in range(n)])
    b = np.array(list(map(float, lines[n + 1].split())))
    x0 = np.array(list(map(float, lines[n + 2].split())))
    TOL = float(lines[n + 3])
    N = int(lines[n + 4])
    return A, b, x0, TOL, N

# Bước 1.1 + 1.2: Kiểm tra tính chéo trội và tính q, lambda theo trường hợp
def kiem_tra_cheo_troi_va_tinh_q_lambda(A):
    import sys

    n = A.shape[0]

    row_dom = True
    col_dom = True
    strictly_row = False
    strictly_col = False

    row_q_list = []
    for i in range(n):
        sum_off = sum(abs(A[i][j]) for j in range(n) if j != i)
        if abs(A[i][i]) < sum_off:
            row_dom = False
        if abs(A[i][i]) > sum_off:
            strictly_row = True
        row_q_list.append(sum_off / abs(A[i][i]))

    col_q_list = []
    for j in range(n):
        sum_off = sum(abs(A[i][j]) for i in range(n) if i != j)
        if abs(A[j][j]) < sum_off:
            col_dom = False
        if abs(A[j][j]) > sum_off:
            strictly_col = True
        col_q_list.append(sum_off / abs(A[j][j]))

    # Trường hợp cả hai đều thỏa mãn
    if row_dom and strictly_row and col_dom and strictly_col:
        while True:
            choice = input("Ma trận chéo trội theo cả hàng và cột. Bạn muốn chọn kiểu nào? (row/col): ").strip().lower()
            if choice == 'row':
                print("row_q_list", row_q_list)
                q = max(row_q_list)
                lam = 1
                return q, lam, 'row'
            elif choice == 'col':
                print("col_q_list", col_q_list)
                q = max(col_q_list)
                diagonals = [abs(A[i][i]) for i in range(n)]
                lam = max(diagonals) / min(diagonals)
                return q, lam, 'col'
            else:
                print("Lựa chọn không hợp lệ. Vui lòng nhập 'row' hoặc 'col'.")
    elif row_dom and strictly_row:
        print("row_q_list", row_q_list)
        q = max(row_q_list)
        lam = 1
        return q, lam, 'row'
    elif col_dom and strictly_col:
        print("row_q_list", col_q_list)
        q = max(col_q_list)
        diagonals = [abs(A[i][i]) for i in range(n)]
        lam = max(diagonals) / min(diagonals)
        return q, lam, 'col'
    else:
        raise ValueError("Ma trận A không chéo trội theo hàng hoặc cột. Không đảm bảo hội tụ theo Jacobi.")

def float_to_fraction_str(x, max_denominator=10000):
    try:
        frac = Fraction(x).limit_denominator(max_denominator)
        # Nếu là số nguyên, chỉ in số nguyên
        if frac.denominator == 1:
            return str(frac.numerator)
        return f"{frac.numerator}/{frac.denominator}"
    except:
        return str(x)

def print_matrix_fraction(mat, name):
    print(f"\n{name} = ")
    rows = []
    for row in mat:
        rows.append([float_to_fraction_str(x) for x in row])
    print(tabulate(rows, tablefmt="fancy_grid"))

def print_vector_fraction(vec, name):
    print(f"\n{name} = ")
    for x in vec:
        print(float_to_fraction_str(x))

def ma_tran_sau_khi_bien_doi(A, b):
    n = len(b)
    B = np.zeros_like(A)
    d = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i == j:
                B[i][j] = 0
            else:
                B[i][j] = -A[i][j] / A[i][i]
        d[i] = b[i] / A[i][i]
    return B, d


def jacobi_theo_sach(A, b, x0, TOL, N):
    global err, x_new
    n = A.shape[0]

    # Bước 1.1 + 1.2
    q, lam, cheo = kiem_tra_cheo_troi_va_tinh_q_lambda(A)
    ord_type = np.inf if cheo == 'row' else 1

    print("\n📏 A là ma trận chéo trội theo:", "HÀNG" if cheo == 'row' else "CỘT")
    print(f"q = {q:.{DECIMALS}f}, lambda = {lam:.{DECIMALS}f}")

    B, d = ma_tran_sau_khi_bien_doi(A, b)
    print_matrix_fraction(B, "B")
    print_vector_fraction(d, "d")

    # Bước 2: TOL'
    tol = TOL * (1 - q) / (lam * q)
    print(f"epsilon = {tol:.{DECIMALS}f}")

    # Bước 3
    k = 1
    x_prev = x0.copy()
    x_first = None
    iteration_log = []

    while k <= N:
        x_new = np.zeros_like(x_prev)

        # Bước 3.1
        for i in range(n):
            s = sum(A[i][j] * x_prev[j] for j in range(n) if j != i)
            x_new[i] = (b[i] - s) / A[i][i]


        if k == 1:
            x_first = x_new.copy()
            # Dự đoán số vòng lặp
            du_doan_vong_lap = np.log((TOL * (1 - q)) / (np.linalg.norm(x_first - x0, ord=ord_type))) / np.log(q)
            print(f"Cần lớn hơn: {du_doan_vong_lap}")
            du_doan_vong_lap = int(np.ceil(du_doan_vong_lap))
            print(f"📈 Dự đoán số vòng lặp cần thiết: {du_doan_vong_lap}")

        # Bước 3.2
        err = np.linalg.norm(x_new - x_prev, ord=ord_type)

        iteration_log.append((k, *x_new, err))

        if err <= tol:
            break

        # Bước 3.3 + 3.4
        k += 1
        x_prev = x_new

    abs_err = err
    rel_err = abs_err / np.linalg.norm(x_new, ord=ord_type)
    post_err = lam * (q / (1 - q)) * abs_err
    pre_err = (q**k / (1 - q)) * np.linalg.norm(x_first - x0, ord=np.inf)

    return x_new, k, iteration_log, abs_err, rel_err, post_err, pre_err
def main():
    A, b, x0, TOL, N = read_input("input.txt")
    try:
        x_final, k, logs, abs_err, rel_err, post_err, pre_err = jacobi_theo_sach(A, b, x0, TOL, N)


        print("\n📘 Quá trình lặp Jacobi:")
        print(tabulate(logs, headers=["Lần lặp"] + [f"x{i + 1}" for i in range(len(x0))] + ["Sai số"], floatfmt=f".{DECIMALS}f",
                       tablefmt="fancy_grid"))

        print("\n🔎 Nghiệm gần đúng cuối cùng:", x_final)
        print("🔁 Số lần lặp:", k)

        print("\n📊 Sai số cuối cùng:")
        print(tabulate([
            ["Sai số tuyệt đối", abs_err],
            ["Sai số tương đối", rel_err],
            ["Sai số hậu nghiệm (3.3)", post_err],
            ["Sai số tiên nghiệm (3.2)", pre_err],
        ], headers=["Loại sai số", "Giá trị"], floatfmt=f".{DECIMALS}e", tablefmt="fancy_grid"))

    except ValueError as e:
        print("❌ Lỗi:", e)

if __name__ == "__main__":
    main()
