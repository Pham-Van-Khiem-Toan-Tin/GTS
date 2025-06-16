import numpy as np

with open("input.txt", "r") as f:
    n = int(f.readline())
    X = np.array([float(i) for i in f.readline().split()]).reshape(n, 1)
    A = np.array([float(i) for i in f.read().split()]).reshape(n, n)

def normalize_L2(x):
    return x / np.linalg.norm(x)

# Phương pháp lũy thừa (in chi tiết từng bước)
def power_method_verbose(A, x0, tol=1e-6, max_iter=1000):
    x = normalize_L2(x0)
    print(x)
    print("=== QUÁ TRÌNH LŨY THỪA ===")
    for i in range(max_iter):
        x_new = np.dot(A, x)
        print(f"Bước {i+1}: Ax = \n{x_new}")
        x_new = normalize_L2(x_new)
        print(f"         Sau chuẩn hóa L2: \n{x_new}\n")
        if np.linalg.norm(x_new - x) < tol:
            print(f"→ Dừng sau {i+1} bước do hội tụ (||x_k+1 - x_k|| < {tol})\n")
            break
        x = x_new
    lam = float(((np.dot(np.dot(x.T, A), x)) / (np.dot(X.T, x)))[0, 0])
    print(f"Giá trị riêng xấp xỉ: {round(lam, 6)}")
    print(f"Véc-tơ riêng chuẩn hóa: \n{x}")
    return lam, x

# Áp dụng
lambda1, v1 = power_method_verbose(A, X)

print("\nGiá trị riêng trội nhất (xấp xỉ):", round(lambda1, 6))
print("Véc-tơ riêng tương ứng (chuẩn hóa chuẩn 2):")
print(v1)