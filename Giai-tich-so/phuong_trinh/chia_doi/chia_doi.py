from tabulate import tabulate

def bisection_method(f, a, b, epsilon, N0):
    """
    Thuật toán chia đôi để tìm nghiệm gần đúng của phương trình f(x) = 0.

    Tham số:
    f       - Hàm số đầu vào f(x).
    a, b    - Biên của khoảng cách ly nghiệm.
    epsilon - Sai số dừng.
    N0      - Số lần lặp tối đa.

    Trả về:
    x       - Nghiệm gần đúng của phương trình.
    """
    if f(a) * f(b) > 0:
        print("Điều kiện khoảng cách ly nghiệm không thỏa mãn. Dừng chương trình.")
        return None
    if f(a) == 0:
        print("Nghiệm chính xác tại a:", a)
        return a
    if f(b) == 0:
        print("Nghiệm chính xác tại b:", b)
        return b

    data = []  # Lưu kết quả từng bước
    i = 0
    while i < N0:
        c = (a + b) / 2
        fc = f(c)
        error = abs(b - a) / 2
        
        data.append([round(a, 7), round(b, 7), round(c, 7), "-" if fc < 0 else "+", round(error, 7)])
        
        if abs(fc) == 0 or error < epsilon:
            break
        
        if f(a) * fc < 0:
            b = c
        else:
            a = c
        
        i += 1
    
    print(tabulate(data, headers=["a", "b", "c", "dấu f(c)", "ε"], tablefmt="grid", floatfmt=".7f"))
    return c

# Hàm ví dụ
def f(x):
    return x**3 -1.5*x**2 + 0.58*x - 0.057

# Thông số đầu vào
a = 0
b = 0.22
epsilon = 5e-3
N0 = 20

# Gọi hàm tìm nghiệm
result = bisection_method(f, a, b, epsilon, N0)
if result is not None:
    print("Nghiệm gần đúng:", round(result, 7))
