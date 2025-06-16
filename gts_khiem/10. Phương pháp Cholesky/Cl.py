import numpy as np

def is_symmetric(A):
    return np.allclose(A, A.T)

def cholesky_decomposition(A):
    n = A.shape[0]
    L = np.zeros_like(A)

    for i in range(n):
        for j in range(i+1):
            if i == j:
                L[i, j] = np.sqrt(A[i, i] - np.sum(L[i, :j]**2))
            else:
                L[i, j] = (A[i, j] - np.sum(L[i, :j] * L[j, :j])) / L[j, j]

    return L

def solve_cholesky(A, b):
    L = cholesky_decomposition(A)
    y = np.linalg.solve(L, b)
    x = np.linalg.solve(L.T, y)
    return x

# Đường dẫn tới file chứa ma trận mở rộng
filename = "augmented_matrix.txt"

# Đọc ma trận mở rộng từ file
augmented_matrix = np.loadtxt(filename)

# Tách ma trận hệ số và vector cột b từ ma trận mở rộng
A = augmented_matrix[:, :-1]
b = augmented_matrix[:, -1]

# Kiểm tra tính đối xứng của ma trận
if not is_symmetric(A):
    print("Ma trận không đối xứng. Áp dụng phương pháp khác để giải hệ phương trình.")
else:
    # Giải hệ phương trình bằng phương pháp Cholesky
    x = solve_cholesky(A, b)

    # In kết quả
    print("Kết quả:")
    for i, xi in enumerate(x):
        print("x", i+1, "=", xi)
