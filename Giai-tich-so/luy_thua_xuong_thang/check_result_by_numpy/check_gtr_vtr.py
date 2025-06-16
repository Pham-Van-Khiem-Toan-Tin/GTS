import numpy as np


with open("input.txt", "r") as f:
    n = int(f.readline())
    A = np.array([float(i) for i in f.read().split()]).reshape(n, n)



eigenvalues, eigenvectors = np.linalg.eig(A)

# In kết quả
print("Trị riêng (λ):")
print(eigenvalues)

print("\nVéc-tơ riêng (cột):")
print(eigenvectors)