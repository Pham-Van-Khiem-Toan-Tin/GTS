from decimal import Decimal, getcontext
from sympy import *
import sys
from tabulate import tabulate
import numpy as np


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
def check_cross_dominant(A):
    cross_dominant = True
    row_dominant = True
    column_dominant = True
    n = len(A)
    for i in range(n):
        sum_row = 0.0
        sum_column = 0.0
        for j in range(n):
            sum_row += abs(A[i,j])
            sum_column += abs(A[j, i])
        if abs(A[i,i]) < sum_row:
            row_dominant = False
            cross_dominant = False
        if abs(A[i,i]) < sum_column:
            column_dominant = False
            cross_dominant = False
    return cross_dominant, row_dominant, column_dominant
def main():
    filename = "input.txt"
    A,b = read_augmented_matrix_from_file(filename)
    print(check_cross_dominant(A))
if __name__ == "__main__":
    main()