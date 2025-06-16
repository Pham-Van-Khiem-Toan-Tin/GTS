import numpy as np

# Đọc ma trận từ file
def read_matrix_from_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        matrix = [list(map(int, line.strip().split())) for line in lines]
    return np.array(matrix)

# Ghi ma trận ra file
def write_matrix_to_file(filename, matrix):
    with open(filename, 'w') as f:
        for row in matrix:
            f.write(' '.join(map(str, row)) + '\n')

# Đọc hai ma trận
matrix1 = read_matrix_from_file('matrix1.txt')
matrix2 = read_matrix_from_file('matrix2.txt')

# Cộng hai ma trận
result = matrix1 + matrix2
# trừ hai ma trận
# result = matrix1 - matrix2

# Ghi kết quả ra file
write_matrix_to_file('result.txt', result)

print("Đã cộng xong hai ma trận và lưu vào file result.txt")
