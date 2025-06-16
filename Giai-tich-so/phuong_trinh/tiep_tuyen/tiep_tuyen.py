import math
from tabulate import tabulate

# Hàm f(x) và đạo hàm cấp 1, cấp 2
def f(x):
    return 1.2 * x**5 - 2.57 * x + 2

def f1(x):
    return 6 * x**4 - 2.57

def f2(x):
    return 24 * x**3

# Đầu vào
a = -1.5
b = -1
eps = 1e-7
N0 = 15

# Kiểm tra điều kiện định lý 2.4
if f(a) * f2(a) > 0:
    x0 = a
    print("Chọn x0 = a = ", x0)
else:
    x0 = b
    print("Chọn x0 = b = ", x0)

# Tính m1 và M2
num_samples = 1000
xs = [a + i * (b - a) / num_samples for i in range(num_samples + 1)]
m1 = min(abs(f1(x)) for x in xs)
M2 = max(abs(f2(x)) for x in xs)

print(f"m1 = {m1}, M2 = {M2}")

# Tính epsilon như trong định lý
epsilon = math.sqrt((2 * m1 * eps) / M2)
print(f"d_(x_n-1) = {epsilon}")

# Lặp Newton-Raphson
i = 0
x_old = x0
table = []
while i < N0:
    x_new = x_old - f(x_old) / f1(x_old)
    dx = abs(x_new - x_old)
    
    table.append([i, x_old, dx])
    
    if dx < epsilon:
        x_old = x_new
        i += 1
        break
    i += 1
    x_old = x_new
x_new = x_old - f(x_old) / f1(x_old)
dx = abs(x_new - x_old)
table.append([i, x_new, dx])

# In bảng kết quả dạng grid
headers = ["Lần lặp (k)", "x_k", "d_xk"]
print(tabulate(table, headers=headers, tablefmt="grid", floatfmt=".15f"))

# In kết luận
print("\nNghiệm gần đúng là:", f"{x_old:.15f}")
print("Sai số là:", f"{abs(dx):.15f}")
