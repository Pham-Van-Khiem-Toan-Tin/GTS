import numpy as np


with open('input.txt', 'r') as f:
    lines = [line.strip() for line in f if line.strip()]
    A = np.array([list(map(float, line.split())) for line in lines])


U, S, Vh = np.linalg.svd(A)

print("Giá trị kỳ dị (singular values):", S)
print("Vector kỳ dị trái (U):")
print(U)
print("Vector kỳ dị phải (Vh):")
print(Vh)