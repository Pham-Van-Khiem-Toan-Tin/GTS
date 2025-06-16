import numpy as np
from tabulate import tabulate

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
# Đệ quy chuyển hàng có aii = 0 bằng hàng có aii khác 0 theo quy tắc từ trên xuống dưới từ trái qua phải
def change_position_row(extended_matrix, lenA, row_handle, col, max_col):
    augmented_matrix = extended_matrix
    row_point = row_handle
    for i in range(row_point, lenA - 1):
        pivot_row = i
        pivot_col = col
        if(augmented_matrix[i][pivot_col] == 0):
            for j in range(i + 1, lenA):
                if(augmented_matrix[j][pivot_col] != 0):
                    pivot_row = j
                    break
            if(pivot_row != i):
                print(f"Đổi chỗ hàng {i + 1} cho hàng {pivot_row + 1}")
                augmented_matrix[[i, pivot_row]] = augmented_matrix[[pivot_row, i]]
                pivot_col += 1
                return augmented_matrix
            else:
                if(pivot_col < max_col - 2):
                    augmented_matrix = change_position_row(extended_matrix, lenA, row_point, pivot_col + 1, max_col)
                else:
                    break
        else:
            break
    return augmented_matrix
        
# gauss thuận đưa hệ về dạng tam giác trên
def gauss_elimination(extended_matrix, m, n):
    augmented_matrix = extended_matrix
    
    for i in range( m-1):
        augmented_matrix = change_position_row(augmented_matrix, m, i, i, n + 1)
        col_factor = 0
        for p in range(n):
            if(abs(augmented_matrix[i,p]) > 1e-10 and augmented_matrix[i,p] != 0):
                col_factor = p
                break
        for j in range(i + 1, m):
            if(augmented_matrix[j, col_factor] != 0):
                factor = augmented_matrix[j, col_factor] / augmented_matrix[i,col_factor]
                print(f"Biến đổi: h{j + 1} - {augmented_matrix[j, col_factor]} * h{i + 1} / {augmented_matrix[i,col_factor]}")
                print(f"Hay Biến đổi: h{j + 1} - {factor} * h{i + 1}")
                augmented_matrix[j] -= factor * augmented_matrix[i]
        print(f"Ma trận sau bước biến đổi {i + 1}:")
        print(tabulate(augmented_matrix, tablefmt="grid", floatfmt=".8f"))
        
        # gauss_elimination(augmented_matrix, lenA, pivot_row + 1, pivot_col + 1)
    print("Ma trận sau quá trình gauss thuận")
    print(tabulate(augmented_matrix, tablefmt="grid", floatfmt=".8f"))
    return augmented_matrix
def infinitely_many_solutions(ladder_matrix, m, n, row_check):
    augmented_matrix = ladder_matrix
    col_split = 0
    row_an = m - row_check
    AB = augmented_matrix[:row_an, :]
    for i in range(len(AB[-1])):
        if(AB[-1, i] != 0):
            col_split = i
            break
    X = AB[:, :col_split + 1]
    X_cp = X
    B = AB[:, n:n + 1]
    B_an = AB[:, col_split + 1:n]
    B_an = -B_an
    B_e = np.concatenate((B, B_an), axis=1)
    AB = np.concatenate((X, B_e), axis=1)
    result = np.full((len(X), len(B_e[0])), 0.0)
    
    for i in range(len(X) - 1, -1, -1):
        matrix_sum = B_e[i]
        for j in range(i + 1, len(X[i])):
            if(j == len(X[i]) - i - 1):
                print(f"Biến đổi: h{i} - {X[i, j]}*x{j + 1}")
                matrix_sum -= X[i, j] * result[len(result) - i - 1]
                X_cp[i, j] = 0
        for j in range(len(X[i])):
            if(X[i,j] != 0):
                print(f"Biến đổi: h{i} / {X[i, j]}")
                result[i] = matrix_sum / X[i,j]
                X_cp[i] = X[i]/X[i,j]
                break
    AXB = np.concatenate((X_cp,result), axis=1)
    equation = np.full(len(AXB), "", dtype=object)
    for i in range(len(AXB)):
        equation_row = ""
        for j in range(len(AXB[i])):
            if(j < col_split):
                equation_row += str(AXB[i,j]) + "*x" + str(j+1) + " + "
            elif(j == col_split):
                equation_row += str(AXB[i,j]) + "*x" + str(j+1) + " = "
            elif(j == col_split + 1):
                equation_row += str(AXB[i,j]) + " + "
            elif(j == n):
                equation_row += str(AXB[i,j]) + "*x" + str(j)
            else:
                equation_row += str(AXB[i,j]) + "*x" + str(j) + " + "
        equation[i] = equation_row
    print("Hệ phương trình theo ẩn:")
    for i in range(len(equation)):
        print(equation[i])
def gauss_back(ladder_matrix, m, n):
    augmented_matrix = ladder_matrix
    result = np.full(n, 0.0)
    for i in range(m):
        check_row = np.full(n + 1, False)
        for j in range(n + 1):
            if(augmented_matrix[i, j] == 0 or abs(augmented_matrix[i, j]) < 1e-10):
                augmented_matrix[i, j] = 0
                check_row[j] = True
        row_x = check_row[:n]
        if(all(row_x) and check_row[n]):
            infinitely_many_solutions(augmented_matrix, m, n, i)
            return "Hệ phương trình có vô số nghiệm"
        elif(all(row_x) and not check_row[n]):
            return "Hệ phương trình vô nghiệm"
    for i in range(m - 1, -1, -1):
        # Lấy giá trị của x[i]
        sum_ = augmented_matrix[i, -1]
        for j in range(i+1, n):
            sum_ -= augmented_matrix[i, j] * result[j]
        
        result[i] = sum_ / augmented_matrix[i, i]
    return result
#hàm tính bằng phương pháp gauss
def solve_gauss(A, b):
    m = len(A)
    n = len(A[0])
    # Tạo ma trận mở rộng
    augmented_matrix = np.concatenate((A, b.reshape(m,1)), axis = 1)
    gauss_elimination_matrix = gauss_elimination(augmented_matrix, m, n)
    print("Kết quả sau quá trình gauss nghịch:")
    result = gauss_back(gauss_elimination_matrix, m, n)
    print(result)
def main():
    # Tên file input
    filename = "augmented_matrix.txt"
    A, b = read_augmented_matrix_from_file(filename)
    solve_gauss(A, b)
if __name__ == "__main__":
    main()