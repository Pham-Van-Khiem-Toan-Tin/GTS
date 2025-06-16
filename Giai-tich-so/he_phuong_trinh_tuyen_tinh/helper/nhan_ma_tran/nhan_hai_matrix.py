import numpy as np

def read_matrix_from_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        matrix = [list(map(int, line.strip().split())) for line in lines]
    return np.array(matrix)

def write_matrix_to_file(filename, matrix):
    with open(filename, 'w') as f:
        for row in matrix:
            f.write(' '.join(map(str, row)) + '\n')

# Đọc hai ma trận từ file
matrix1 = read_matrix_from_file('matrix1.txt')
matrix2 = read_matrix_from_file('matrix2.txt')

# Kiểm tra điều kiện nhân ma trận
if matrix1.shape[1] != matrix2.shape[0]:
    print("Không thể nhân hai ma trận: số cột của ma trận 1 phải bằng số dòng của ma trận 2.")
else:
    # Nhân hai ma trận
    result = np.dot(matrix1, matrix2)
    
    # Ghi kết quả ra file
    write_matrix_to_file('result.txt', result)
    print("Đã nhân xong hai ma trận và lưu vào file result.txt")
