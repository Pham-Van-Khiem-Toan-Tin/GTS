import numpy as np
from tabulate import tabulate

from normalization import normalize_L2

# Đọc dữ liệu từ file input.txt.txt
with open("input.txt", "r") as f:
    n = int(f.readline())
    x0 = np.array([float(i) for i in f.readline().split()]).reshape(n, 1)
    A = np.array([float(i) for i in f.read().split()]).reshape(n, n)


# Phương pháp lũy thừa có in chi tiết
def power_method_verbose(A, x0, tol=1e-6, max_iter=1000):
    x = normalize_L2(x0)
    print("=== BẮT ĐẦU PHƯƠNG PHÁP LŨY THỪA ===")
    for i in range(max_iter):
        x_new = A @ x
        print(f"\nBước {i + 1}:")
        print("→ A * x =")
        print(x_new)
        x_new = normalize_L2(x_new)
        print("→ Chuẩn hóa (L2):")
        print(x_new)
        diff = np.linalg.norm(x_new - x)
        print(f"→ Sai số ||x_k+1 - x_k|| = {diff:.6e}")
        lam = float(((x.T @ A @ x) / (x.T @ x))[0, 0])
        print(f"lamda = {lam}")
        if diff < tol:
            print(f"✅ Hội tụ sau {i + 1} bước (sai số < {tol})")
            break
        x = x_new

    print("\n🎯 KẾT QUẢ LŨY THỪA")
    print(f"→ Giá trị riêng gần đúng: λ ≈ {round(lam, 6)}")
    print("→ Véc-tơ riêng tương ứng:")
    print(x)
    return lam, x


# Phương pháp xuống thang
def deflation(A, lam, v):
    return A - lam * (v @ v.T)


# Tổng quát: tìm tất cả trị riêng bằng lũy thừa + xuống thang
def compute_all_eigenpairs(A, x0, tol=1e-6, max_iter=5):
    n = A.shape[0]
    eigenvalues = []
    eigenvectors = []

    print("===== BẮT ĐẦU TÌM TOÀN BỘ TRỊ RIÊNG =====")
    for k in range(n):
        print(f"\n==============================")
        print(f"        LẦN THỨ {k+1}")
        print("==============================")
        lam, v = power_method_verbose(A, x0, tol, max_iter)
        eigenvalues.append(lam)
        eigenvectors.append(v)
        print(f"\n→ Ghi nhận trị riêng λ{k+1} ≈ {round(lam, 6)}")

        if k < n - 1:
            print("\n→ Ma trận trước khi xuống thang:")
            print(A)
            A = deflation(A, lam, v)
            print("→ Ma trận sau khi xuống thang (A - λ vvᵀ):")
            print(A)
    return eigenvalues, eigenvectors

# Gọi hàm chính
eigenvalues, eigenvectors = compute_all_eigenpairs(A, x0)

# In tổng kết cuối
print("\n===== TỔNG KẾT CUỐI CÙNG =====")
for i, (lam, v) in enumerate(zip(eigenvalues, eigenvectors), start=1):
    print(f"\nλ{i} ≈ {round(lam, 6)}")
    print(f"v{i}.T = {v.T}")

# In giá trị kỳ dị (singular values)
print("Giá trị kỳ dị là :", [np.sqrt(x) for x in eigenvalues])

# In eigenvectors đẹp
print(tabulate(eigenvectors, tablefmt="fancy_grid"))
