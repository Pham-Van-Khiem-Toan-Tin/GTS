from tabulate import tabulate

def secant_method(f, df, ddf, a, b, epsilon, N0):
    """
    Thuật toán dây cung để tìm nghiệm gần đúng của phương trình f(x) = 0.

    Tham số:
    f       - Hàm số đầu vào f(x).
    df      - Đạo hàm cấp 1 của f(x).
    ddf     - Đạo hàm cấp 2 của f(x).
    a, b    - Biên của khoảng cách ly nghiệm.
    epsilon - Sai số dừng.
    N0      - Số lần lặp tối đa.

    Trả về:
    x       - Nghiệm gần đúng của phương trình hoặc thông báo dừng nếu điều kiện không thỏa mãn.
    """
    # Bước 1: Kiểm tra điều kiện đầu trong Định lý 2.3
    if f(a) * f(b) > 0:
        print("Điều kiện khoảng cách ly nghiệm không thỏa mãn. Dừng chương trình.")
        return None
    if (df(a) * df(b) < 0) or (ddf(a) * ddf(b) < 0):
        print("Đạo hàm f' hoặc f'' đổi dấu trên [a, b]. Dừng chương trình.")
        return None
    
    # Xác định x0 và d
    if f(a) * ddf(a) < 0:
        x0, d = a, b
    else:
        x0, d = b, a
    
    # Bước 3: Tính m1 và M1
    m1 = min(abs(df(x)) for x in [a, b])
    M1 = max(abs(df(x)) for x in [a, b])
    print(f"m1 = {m1}, M1 = {M1}")
    # Gán epsilon'
    epsilon_prime = (m1 * epsilon) / (M1 - m1)
    
    i = 0  # Khởi tạo số lần lặp
    data = []  # Lưu kết quả từng bước
    d_xk = 0  # Khởi tạo d_xk
    while i <= N0:
        # Bước 4.1: Tính dx
        dx = - (f(x0) * (x0 - d)) / (f(x0) - f(d))
        x1 = x0 + dx
        
        # Lưu lại giá trị
        if i==N0:
            data.append([i, x0, ""])
        else:
            data.append([i, x0, dx])
        
        # Bước 4.2: Kiểm tra điều kiện dừng
        if abs(dx) < epsilon_prime:
            print(tabulate(data, headers=["k", "xk", "d_xk"], tablefmt="grid",floatfmt=".9f"))
            sai_so = abs(d_xk) * (M1 - m1) / m1
            return x0,sai_so  # Nghiệm gần đúng
        
        # Bước 4.3: Cập nhật giá trị
        i += 1
        if i <= N0:
            x0 = x1
            d_xk = dx

    # Bước 5: Dừng chương trình do số lần lặp đạt tối đa
    print(tabulate(data, headers=["k", "xk", "d_xk"], tablefmt="grid",floatfmt=".9f"))
    sai_so = abs(d_xk) * (M1 - m1) / m1
    return x0,sai_so  # Trả về nghiệm gần đúng sau N0 lần lặp

# Hàm ví dụ

def f(x):
    return x**3 - 2*x -5

def df(x):
    return 3*x**2 - 2

def ddf(x):
    return 6*x

# Thông số đầu vào
a = 2.0
b = 3.0
epsilon = 1e-8
N0 = 100

# Gọi hàm tìm nghiệm
result = secant_method(f, df, ddf, a, b, epsilon, N0)
if result is not None:
    x_approx, sai_so = result
    print("Nghiệm gần đúng:", round(x_approx, 10))
    print("Sai số:", sai_so)