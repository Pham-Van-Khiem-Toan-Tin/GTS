import numpy as np

def read_matrix_from_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        matrix = [list(map(float, line.strip().split())) for line in lines]
    return np.array(matrix)

def write_matrix_to_file(filename, matrix):
    with open(filename, 'w') as f:
        for row in matrix:
            f.write(' '.join(map(str, row)) + '\n')

# Đọc ma trận A từ file
A = read_matrix_from_file('matrix1.txt')

# Nhập hệ số a
a = float(input("Nhập hệ số a: "))

# Tạo ma trận đơn vị I cùng kích thước
if A.shape[0] != A.shape[1]:
    print("Không thể tạo ma trận đơn vị: A không phải là ma trận vuông.")
else:
    I = np.identity(A.shape[0])

    # Tính toán: C = A + a * I
    C = A + a * I

    # Ghi kết quả ra file
    write_matrix_to_file('result.txt', C)
    print("Đã thực hiện phép tính A + a*I và lưu vào result.txt")
