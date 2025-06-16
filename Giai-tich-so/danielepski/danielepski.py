import numpy as np
from fractions import Fraction
from tabulate import tabulate

def read_matrix(filename):
    with open(filename, "r") as f:
        data = [list(map(Fraction, line.strip().split())) for line in f if line.strip()]
    return np.array(data, dtype=object)

def print_matrix(mat, name="", precision=0):
    rows = []
    for row in mat:
        if precision > 0:
            rows.append([str(round(float(x), precision)) for x in row])
        else:
            rows.append([str(x) for x in row])
    print(f"\n{name}")
    print(tabulate(rows, tablefmt="fancy_grid"))

def danilevski_process(A):
    n = A.shape[0]
    B = A.copy()
    transforms = []
    print_matrix(B, f"A^{(1)} (Ban đầu)")
    for k in range(n-1, 0, -1):
        print(f"\n{'='*30}\nBước {n-k}: Đưa a[{k+1},{k}] về khác 0")
        if B[k, k-1] == 0:
            # Nếu a[k, k-1] = 0, tìm cột j < k-1 sao cho B[k, j] != 0 để hoán vị
            swap_col = None
            for j in range(k-1):
                if B[k, j] != 0:
                    swap_col = j
                    break
            if swap_col is not None:
                print(f"Hoán vị dòng/cột {k} <-> {swap_col+1} vì B[{k+1},{k}] = 0 và B[{k+1},{swap_col+1}] ≠ 0")
                B[[k-1, swap_col], :] = B[[swap_col, k-1], :]
                B[:, [k-1, swap_col]] = B[:, [swap_col, k-1]]
                print_matrix(B, f"A^{(n-k+1)} sau hoán vị")
            else:
                print("Không thể đưa về dạng Frobenius (dưới chéo đều 0)")
                return

        # Xây dựng ma trận M_k
        a = B[k, k-1]
        Mk = np.eye(n, dtype=object)
        Mk[k - 1, :] = B[k, :]
        Mk_inv = np.eye(n, dtype=object)
        Mk_inv[k - 1, :] = -B[k, :] / a
        Mk_inv[k - 1, k - 1] = Fraction(1, a)
        transforms.insert(0, Mk)
        print_matrix(Mk, f"M_{n-k+1}")
        print_matrix(Mk_inv, f"M_{n-k+1}^(-1)")

        # Tính A^{(k+1)} = Mk * B * Mk_inv
        B = Mk @ B @ Mk_inv
        print_matrix(B, f"A^{(n-k+2)} = M_{n-k+1} * A^{(n-k+1)} * M_{n-k+1}^(-1)")

    # In ra ma trận Frobenius cuối cùng
    print(f"\n{'='*30}\nMa trận Frobenius thu được:")
    print_matrix(B, f"Frobenius (A^{n+1})")
    return B, transforms

def characteristic_polynomial(F):
    n = F.shape[0]
    # Đa thức đặc trưng: lambda^n + a1*lambda^{n-1} + ... + an
    coeffs = [Fraction(1)]
    coeffs += [-F[0, j] for j in range(n)]
    print("\nĐa thức đặc trưng của A là:")
    s = "P(λ) = λ^" + str(n)
    for i in range(1, n+1):
        sign = '+' if coeffs[i] >= 0 else '-'
        val = abs(coeffs[i])
        if val == 1:
            val = ""
        if n-i == 0:
            s += f" {sign} {val}"
        elif n-i == 1:
            s += f" {sign} {val}λ"
        else:
            s += f" {sign} {val}λ^{n-i}"
    print(" ", s, "= 0")
    print("Hệ số (từ cao xuống thấp):", [str(x) for x in coeffs])
    return coeffs

def get_eigenvectors(transforms, eigenvalue):
    import sympy
    n = len(transforms) + 1
    Y = np.array([eigenvalue**i for i in range(n-1, -1, -1)], dtype=object).reshape((n, 1))
    P = np.eye(n, dtype=object)
    for Mk in transforms:
        P = P @ Mk
    P_sym = sympy.Matrix(P)
    P_inv = P_sym.inv()
    Y_sym = sympy.Matrix(Y)
    X_sym = P_inv * Y_sym
    X = np.array(X_sym.tolist(), dtype=object)
    X_flat = X.flatten()
    idx = np.argmax([abs(float(v)) > 1e-10 for v in X_flat])
    if X_flat[idx] != 0:
        X = X / X_flat[idx]

    return X



def main():
    print("PHƯƠNG PHÁP DANILEVSKI - TRÌNH BÀY CHUẨN SÁCH")
    print("="*40)
    A = read_matrix("input.txt")
    print_matrix(A, "Ma trận đầu vào (A)")
    F, transforms = danilevski_process(A)
    coeffs = characteristic_polynomial(F)

    # Tìm nghiệm thực của đa thức đặc trưng
    from sympy import Poly, Symbol, solveset, S
    x = Symbol('λ')
    poly = Poly([float(c) for c in coeffs], x)
    roots = poly.all_roots()
    roots_real = [r.evalf() for r in roots]
    print("\nCác nghiệm thực của đa thức đặc trưng (trị riêng):", roots_real)

    print("\n===> VECTOR RIÊNG THEO ĐÚNG SÁCH (cho mỗi trị riêng thực):")
    for lam in roots_real:
        v = get_eigenvectors(transforms, lam)
        print(f"\nTrị riêng λ = {lam}:")
        print_matrix(v, f"Vector riêng (ứng với λ = {lam})")

if __name__ == "__main__":
    main()
