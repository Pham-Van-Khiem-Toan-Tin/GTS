import numpy as np
from tabulate import tabulate
# Hàm giải hệ phương trình tuyến tính bằng phương pháp Gauss-Jordan
def solve_gauss_jordan(A, b, c):
    n = len(A)
    augmented_matrix = np.concatenate((A, b.reshape(n, 1)), axis=1)

    for i in range(n):
        # Tìm pivot trong cột hiện tại
        pivot_row = i
        for j in range(i+1, n):
            if abs(augmented_matrix[j, i]) > abs(augmented_matrix[pivot_row, i]):
                pivot_row = j

        # Hoán đổi hàng pivot
        augmented_matrix[[i, pivot_row]] = augmented_matrix[[pivot_row, i]]

        # Chia hàng pivot cho pivot để đưa pivot về giá trị 1
        pivot = augmented_matrix[i, i]
        augmented_matrix[i] /= pivot

        # Loại bỏ các ẩn số khác hàng pivot
        for j in range(n):
            if j != i:
                factor = augmented_matrix[j, i]
                augmented_matrix[j] -= factor * augmented_matrix[i]

        if c:
            print("Phép biến đổi lần thứ ", i + 1 )
            A_loop = augmented_matrix[:, :-1]
            B_loop = augmented_matrix[:, -1]
            print(tabulate(np.round(augmented_matrix, 6), tablefmt="grid"))

    # Trích xuất ma trận hệ số và ma trận tự do từ ma trận mở rộng
    A = augmented_matrix[:, :-1]
    b = augmented_matrix[:, -1]

    return A, b

# Hàm đọc ma trận mở rộng từ file
def read_augmented_matrix_from_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        rows = len(lines)
        cols = len(lines[0].split()) - 1  # Trừ 1 cột cuối cùng là cột tự do
        matrix = np.zeros((rows, cols))
        b = np.zeros(rows)
        for i in range(rows):
            line_values = [float(x) for x in lines[i].split()]
            matrix[i] = line_values[:-1]
            b[i] = line_values[-1]
        return matrix, b

# Đường dẫn tới file chứa ma trận mở rộng
filename = "augmented_matrix.txt"

# Đọc ma trận mở rộng từ file
A, b = read_augmented_matrix_from_file(filename)

# Giải hệ phương trình bằng phương pháp Gauss-Jordan
A, b = solve_gauss_jordan(A, b, True)

# In ma trận hệ số sau khi áp dụng phương pháp Gauss-Jordan
print("\n🔹 Ma trận hệ số sau Gauss-Jordan:")
print(tabulate(np.round(A, 6), tablefmt="grid"))
# In ma trận tự do sau khi áp dụng phương pháp Gauss-Jordan
print("\n🔹 Ma trận tự do sau Gauss-Jordan:")
print(tabulate(np.round(b.reshape(-1, 1), 6), tablefmt="grid"))

# Tính rank của ma trận hệ số
rank_A = np.linalg.matrix_rank(A)

# Tính rank của ma trận mở rộng bậc thang
rank_augmented_matrix = np.linalg.matrix_rank(np.concatenate((A, b.reshape(len(b), 1)), axis=1))

# In ra rank của ma trận hệ số và ma trận mở rộng
print(f"\n🔹 Rank của ma trận hệ số: {rank_A}")
print(f"🔹 Rank của ma trận mở rộng: {rank_augmented_matrix}")
# Giải hệ phương trình bằng phương pháp Gauss-Jordan
A, b = solve_gauss_jordan(A, b, False)

# In nghiệm của hệ phương trình
print("\n🔹 Nghiệm của hệ phương trình:")
for i in range(len(b)):
    print(f"  x{i+1} = {b[i]:.6f}")

