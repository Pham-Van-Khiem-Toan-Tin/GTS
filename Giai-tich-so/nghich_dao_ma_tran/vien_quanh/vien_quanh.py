import numpy as np
from tabulate import tabulate
from fractions import Fraction

def read_input(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    n = int(lines[0])
    A = [list(map(Fraction, line.strip().split())) for line in lines[1:n+1]]
    return np.array(A, dtype=object)

def chia_khoi(A):
    n = A.shape[0]
    k = n - 1
    α11 = A[:k, :k]
    α12 = A[:k, k:]
    α21 = A[k:, :k]
    α22 = A[k:, k:]
    return α11, α12, α21, α22

def format_matrix(mat):
    # Tìm độ dài lớn nhất của từng cột để canh chỉnh
    str_mat = [[str(cell) for cell in row] for row in mat]
    col_widths = [max(len(row[i]) for row in str_mat) for i in range(len(str_mat[0]))]

    # Căn giữa từng phần tử
    for i in range(len(str_mat)):
        for j in range(len(str_mat[i])):
            str_mat[i][j] = str_mat[i][j].center(col_widths[j])

    return tabulate(str_mat, tablefmt="fancy_grid", stralign="center")


def nghich_dao_vien_quanh(A):
    n = A.shape[0]
    if n == 1:
        return np.array([[1 / A[0, 0]]], dtype=object)

    α11, α12, α21, α22 = chia_khoi(A)

    print("\n📋 Chia ma trận thành các khối:")
    print("α₁₁:\n", format_matrix(α11))
    print("α₁₂:\n", format_matrix(α12))
    print("α₂₁:\n", format_matrix(α21))
    print("α₂₂:\n", format_matrix(α22))

    α11_inv = nghich_dao_vien_quanh(α11)
    print("\nα₁₁⁻¹:\n", format_matrix(α11_inv))

    X = α11_inv @ α12
    Y = α21 @ α11_inv
    θ = α22 - α21 @ α11_inv @ α12
    θ_inv = np.linalg.inv(θ.astype(np.float64))
    θ_inv = np.array([[Fraction.from_float(x).limit_denominator() for x in row] for row in θ_inv], dtype=object)

    print("\nX = α₁₁⁻¹ · α₁₂:\n", format_matrix(X))
    print("Y = α₂₁ · α₁₁⁻¹:\n", format_matrix(Y))
    print("θ = α₂₂ - α₂₁ · α₁₁⁻¹ · α₁₂:\n", format_matrix(θ))
    print("θ⁻¹:\n", format_matrix(θ_inv))

    β11 = α11_inv + X @ θ_inv @ Y
    β12 = -X @ θ_inv
    β21 = -θ_inv @ Y
    β22 = θ_inv

    top = np.hstack((β11, β12))
    bottom = np.hstack((β21, β22))
    A_inv = np.vstack((top, bottom))

    print(f"\n🧮 Ma trận nghịch đảo ở cấp {n} :\n", format_matrix(A_inv))
    return A_inv

def main():
    A = read_input("input.txt")

    print("\n✅ Bắt đầu tính A⁻¹ bằng phương pháp viền quanh:")
    A_inv = nghich_dao_vien_quanh(A)

    print("\n📄 Ma trận nghịch đảo A⁻¹ cuối cùng:")
    print(format_matrix(A_inv))

    print("\n🧪 Kiểm tra A · A⁻¹ ≈ I:")
    approx_I = np.round(np.array(A, dtype=float) @ np.array(A_inv, dtype=float), 6)
    print(tabulate(approx_I, tablefmt="fancy_grid"))

if __name__ == "__main__":
    main()
