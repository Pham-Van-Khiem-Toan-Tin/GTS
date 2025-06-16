import numpy as np

# Đọc ma trận bất kỳ từ file
with open('input.txt', 'r') as f:
    lines = [line.strip() for line in f if line.strip()]
    A = np.array([list(map(float, line.split())) for line in lines])

# Tính A_nhan_AT
result = np.dot(A.T, A)

# Ghi kết quả ra output.txt
with open('output.txt', 'w') as f:
    for row in result:
        f.write(' '.join(str(round(val, 6)) for val in row) + '\n')
