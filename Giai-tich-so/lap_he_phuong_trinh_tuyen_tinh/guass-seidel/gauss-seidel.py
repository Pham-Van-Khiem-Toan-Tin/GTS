from fractions import Fraction

import numpy as np
from tabulate import tabulate

DECIMALS = 6

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

def kiem_tra_cheo_troi_va_tinh_s_mu(A):
    n = A.shape[0]
    row_dom = True
    col_dom = True
    strictly_row = False
    strictly_col = False

    # Kiểm tra chéo trội hàng
    for i in range(n):
        sum_off = sum(abs(A[i][j]) for j in range(n) if j != i)
        if abs(A[i][i]) < sum_off:
            row_dom = False
        if abs(A[i][i]) > sum_off:
            strictly_row = True

    # Kiểm tra chéo trội cột
    for j in range(n):
        sum_off = sum(abs(A[i][j]) for i in range(n) if i != j)
        if abs(A[j][j]) < sum_off:
            col_dom = False
        if abs(A[j][j]) > sum_off:
            strictly_col = True

    # Nếu cả hai cùng chéo trội
    if row_dom and strictly_row and col_dom and strictly_col:
        while True:
            choice = input("Ma trận chéo trội theo cả hàng và cột. Bạn muốn chọn kiểu nào? (row/col): ").strip().lower()
            if choice == 'row':
                # Trội hàng
                s = 0
                mu_list = []
                for i in range(n):
                    sum_less = sum(abs(A[i][j]) for j in range(i))
                    sum_greater = sum(abs(A[i][j]) for j in range(i+1, n))
                    mu = sum_less / (abs(A[i][i]) - sum_greater) if abs(A[i][i]) > sum_greater else float('inf')
                    mu_list.append(mu)
                mu = max(mu_list)
                return s, mu, 'row'
            elif choice == 'col':
                # Trội cột
                s_list = []
                mu_list = []
                for j in range(n):
                    sum_greater = sum(abs(A[i][j]) for i in range(j+1, n))
                    sum_less = sum(abs(A[i][j]) for i in range(j))
                    s_j = (1 / abs(A[j][j])) * sum_greater if abs(A[j][j]) != 0 else float('inf')
                    mu_j = sum_less / (abs(A[j][j]) - sum_greater) if abs(A[j][j]) > sum_greater else float('inf')
                    s_list.append(s_j)
                    mu_list.append(mu_j)
                s = max(s_list)
                mu = max(mu_list)
                return s, mu, 'col'
            else:
                print("Lựa chọn không hợp lệ. Vui lòng nhập 'row' hoặc 'col'.")
    # Chỉ chéo trội hàng
    elif row_dom and strictly_row:
        s = 0
        mu_list = []
        for i in range(n):
            sum_less = sum(abs(A[i][j]) for j in range(i))
            sum_greater = sum(abs(A[i][j]) for j in range(i+1, n))
            mu = sum_less / (abs(A[i][i]) - sum_greater) if abs(A[i][i]) > sum_greater else float('inf')
            mu_list.append(mu)
        mu = max(mu_list)
        return s, mu, 'row'
    # Chỉ chéo trội cột
    elif col_dom and strictly_col:
        s_list = []
        mu_list = []
        for j in range(n):
            sum_greater = sum(abs(A[i][j]) for i in range(j+1, n))
            sum_less = sum(abs(A[i][j]) for i in range(j))
            s_j = (1 / abs(A[j][j])) * sum_greater if abs(A[j][j]) != 0 else float('inf')
            mu_j = sum_less / (abs(A[j][j]) - sum_greater) if abs(A[j][j]) > sum_greater else float('inf')
            s_list.append(s_j)
            mu_list.append(mu_j)
        s = max(s_list)
        mu = max(mu_list)
        return s, mu, 'col'
    else:
        raise ValueError("Ma trận A không chéo trội theo hàng hoặc cột. Không đảm bảo hội tụ theo Gauss-Seidel.")

def tinh_saiso_gauss_seidel(s, mu, cheo, x_new, x_prev, x_first, x0, k, ord_type):
    if cheo == 'row':
        # Theo sách: ord_type = np.inf
        tien_nghiem = (mu**k) / (1 - mu) * np.linalg.norm(x_first - x0, ord=ord_type)
        hau_nghiem = mu / (1 - mu) * np.linalg.norm(x_new - x_prev, ord=ord_type)
    else:  # cheo == 'col'
        # Theo sách: ord_type = 1
        tien_nghiem = (mu**k) / ((1 - s)*(1 - mu)) * np.linalg.norm(x_first - x0, ord=ord_type)
        hau_nghiem = mu / ((1 - s)*(1 - mu)) * np.linalg.norm(x_new - x_prev, ord=ord_type)
    return tien_nghiem, hau_nghiem


def tach_B1_B2_d(A, b):
    n = len(b)
    B1 = np.zeros_like(A)
    B2 = np.zeros_like(A)
    d = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j < i:
                B1[i, j] = -A[i, j] / A[i, i]
            elif j > i:
                B2[i, j] = -A[i, j] / A[i, i]
        d[i] = b[i] / A[i, i]
    return B1, B2, d


def gauss_seidel_theo_sach(A, b, x0, TOL, N):
    n = A.shape[0]
    s, mu, cheo = kiem_tra_cheo_troi_va_tinh_s_mu(A)
    ord_type = np.inf if cheo == 'row' else 1

    print("\n📏 A là ma trận chéo trội theo:", "HÀNG" if cheo == 'row' else "CỘT")
    print(f"s = {s:.{DECIMALS}f}, μ = {mu:.{DECIMALS}f}")




    # Tách và in ra ma trận lặp
    B1, B2, d = tach_B1_B2_d(A, b)
    print_matrix_fraction(B1, "B₁")
    print_matrix_fraction(B2, "B₂")
    print_vector_fraction(d, "d")

    tol = TOL * (1 - s) * (1 - mu)/mu
    print(f"epsilon = {tol:.{DECIMALS}f}")
    k = 1
    x_prev = x0.copy()
    x_first = None
    iteration_log = []

    while k <= N:
        x_new = x_prev.copy()

        for i in range(n):
            sum1 = sum(A[i][j] * x_new[j] for j in range(i))
            sum2 = sum(A[i][j] * x_prev[j] for j in range(i + 1, n))
            x_new[i] = (b[i] - sum1 - sum2) / A[i][i]

        if k == 1:
            x_first = x_new.copy()
            # Sai số tiên nghiệm dùng ||x_first - x0|| đúng chuẩn
            if mu < 1:
                if cheo == 'row':
                    du_doan_vong_lap = np.log((TOL * (1 - mu)) / (np.linalg.norm(x_first - x0, ord=ord_type))) / np.log(mu)
                else:
                    du_doan_vong_lap = np.log((TOL * (1 - s) * (1 - mu)) / (np.linalg.norm(x_first - x0, ord=ord_type))) / np.log(mu)
                du_doan_vong_lap = int(np.ceil(du_doan_vong_lap))
                print(f"📈 Dự đoán số vòng lặp cần thiết: {du_doan_vong_lap}")

        err = np.linalg.norm(x_new - x_prev, ord=ord_type)
        iteration_log.append((k, *x_new, err))

        if err <= tol:
            break

        k += 1
        x_prev = x_new.copy()

    abs_err = err
    rel_err = abs_err / np.linalg.norm(x_new, ord=ord_type)
    tien_nghiem, hau_nghiem = tinh_saiso_gauss_seidel(s, mu, cheo, x_new, x_prev, x_first, x0, k, ord_type)

    return x_new, k, iteration_log, abs_err, rel_err, hau_nghiem, tien_nghiem

def main():
    A, b, x0, TOL, N = read_input("input.txt")
    try:
        x_final, k, logs, abs_err, rel_err, post_err, pre_err = gauss_seidel_theo_sach(A, b, x0, TOL, N)

        print("\n📘 Quá trình lặp Gauss-Seidel:")
        print(tabulate(logs, headers=["Lần lặp"] + [f"x{i + 1}" for i in range(len(x0))] + ["Sai số"], floatfmt=f".{DECIMALS}f",
                       tablefmt="fancy_grid"))

        print("\n🔎 Nghiệm gần đúng cuối cùng:", x_final)
        print("🔁 Số lần lặp:", k)

        print("\n📊 Sai số cuối cùng:")
        print(tabulate([
            ["Sai số tuyệt đối", abs_err],
            ["Sai số tương đối", rel_err],
            ["Sai số hậu nghiệm", post_err],  # dùng biến post_err mới = hau_nghiem
            ["Sai số tiên nghiệm", pre_err],  # dùng biến pre_err mới = tien_nghiem
        ], headers=["Loại sai số", "Giá trị"], floatfmt=f".{DECIMALS}e", tablefmt="fancy_grid"))
    except ValueError as e:
        print("❌ Lỗi:", e)

if __name__ == "__main__":
    main()
