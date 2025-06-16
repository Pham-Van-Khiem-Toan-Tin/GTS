import numpy as np
from tabulate import tabulate


def read_augmented_matrix_from_file(filename):
    matrix_A = []
    matrix_B = []
    current_matrix = matrix_A
    with open(filename, "r") as file:
        for line in file:
            if line.strip() == "":
                current_matrix = matrix_B
                continue
            row = [float(num) for num in line.split()]
            current_matrix.append(row)
    matrix = np.array(matrix_A)
    b = np.array(matrix_B)
    return matrix, b


def find_pivot(Ab, m, nx, pivot):
    augmented_matrix = Ab
    pivot_row = 0
    pivot_col = 0
    max_element = 0
    break_all_loop = False
    for i in range(m):
        if break_all_loop:
            break
        if pivot[i]:
            continue
        for j in range(nx):
            if abs(augmented_matrix[i, j]) == 1:
                pivot[i] = True
                pivot_row = i
                pivot_col = j
                break_all_loop = True
                break
            if abs(augmented_matrix[i, j]) > max_element:
                pivot_row = i
                pivot_col = j
                max_element = abs(augmented_matrix[i, j])
    return pivot_row, pivot_col


def show_solution(Ab, nb):
    augmented_matrix = Ab
    augmented_matrix_nonzero = augmented_matrix[[i for i, x in enumerate(augmented_matrix) if x.any()]]
    A = augmented_matrix_nonzero[:, :-nb]
    b = augmented_matrix_nonzero[:, -nb:]
    for i in range(len(A)):
        equation_row = ""
        VT = ""
        VP = ""
        if len(b[i]) == 1:
            VP = b[i]
        else:
            for k in range(len(b[i])):
                if k == 0:
                    VP = "(" + str(b[i, k]) + ","
                elif k == len(b[i]) - 1:
                    VP += str(b[i, k]) + ")"
        for j in range(len(A[i])):
            if abs(A[i, j] - 1.0) < 1e-10:
                VT += "x" + str(j + 1)
            elif A[i, j] != 1 and A[i, j] != 0:
                VP += " + " + str(-1 * A[i, j]) + "*x" + str(j + 1)
        equation_row = VT + " = " + VP
        print(equation_row)


def check_have_row_zero(Ab, m, nx, nb):
    augmented_matrix = Ab
    A = augmented_matrix[:, :-nb]
    n = len(A[0])
    b = augmented_matrix[:, -nb:]
    check_b = np.full(len(b), False)
    for i in range(len(b)):
        check_row_b = np.full(len(b[i]), False)
        for j in range(len(b[i])):
            if b[i, j] == 0:
                check_row_b[j] = True
        if all(check_row_b):
            check_b[i] = True
    for i in range(m):
        check_row = np.full(nx, False)
        for j in range(n):
            if A[i, j] == 0:
                check_row[j] = True
        if all(check_row):
            if check_b[i]:
                print("Hệ phương trình có vô số nghiệm")
                return False, False
            else:
                print("Hệ phương trình vô nghiệm")
                return False, True
    return True, True


def solve_gauss_jd(A, b):
    m = len(A)  # Số hàng
    # Số phần tử x
    nx = len(A[0])
    nb = len(b[0])
    show_solution_ex = False
    augmented_matrix = np.concatenate((A, b), axis=1)
    pivot_row_list = np.full(m, False)
    pivot = np.full((m, 2), 0)
    print("Ma trận ban đầu:")
    print(tabulate(augmented_matrix, tablefmt="grid", floatfmt=".8f"))
    have_zero, no_solution = check_have_row_zero(augmented_matrix, m, nx, nb)
    if not have_zero and not no_solution:
        augmented_matrix = augmented_matrix[~np.all(augmented_matrix == 0, axis=1)]
        m = len(augmented_matrix)
        show_solution_ex = True
    elif not have_zero and no_solution:
        return
    for i in range(m):
        check_zero, check_solution = check_have_row_zero(augmented_matrix, m, nx, nb)
        if not check_zero:
            if not check_solution:
                show_solution(augmented_matrix, nb)
            return
        pivot_row, pivot_col = find_pivot(augmented_matrix, m, nx, pivot_row_list)
        pivot_row_list[pivot_row] = True
        pivot[i] = [pivot_row, pivot_col]
        for j in range(m):
            if j == pivot_row:
                continue
            print(
                f"Biến đổi h{j + 1} = h{j + 1} - ({augmented_matrix[j, pivot_col]}/{augmented_matrix[pivot_row, pivot_col]}) * h{pivot_row + 1}"
            )
            augmented_matrix[j] -= (
                augmented_matrix[j, pivot_col]
                * augmented_matrix[pivot_row]
                / augmented_matrix[pivot_row, pivot_col]
            )
            for k in range(len(augmented_matrix[j])):
                if(abs(augmented_matrix[j, k]) < 1e-10):
                    augmented_matrix[j, k] = 0.0
        if augmented_matrix[pivot_row, pivot_col] != 1:
            print(
                f"Chuẩn hoá: h{pivot_row + 1} = h{pivot_row + 1}/{augmented_matrix[pivot_row, pivot_col]}"
            )
            augmented_matrix[pivot_row] = (
                augmented_matrix[pivot_row] / augmented_matrix[pivot_row, pivot_col]
            )
        print(tabulate(augmented_matrix, tablefmt="grid", floatfmt=".8f"))
    if show_solution_ex:
        show_solution(augmented_matrix)
        return
    print("Nghiệm của hệ phương trình là: (Tự nhìn ma trận cuối cùng để kết luận nghiệm là x mấy)")
    result = augmented_matrix[:, -1].reshape(1, -1)
    print(result)


def main():
    filename = "augmented_matrix.txt"
    A, b = read_augmented_matrix_from_file(filename)
    solve_gauss_jd(A, b)


if __name__ == "__main__":
    main()
