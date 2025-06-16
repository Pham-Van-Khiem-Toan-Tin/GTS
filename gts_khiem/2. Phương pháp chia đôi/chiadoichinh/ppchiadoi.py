from tabulate import tabulate
import numpy as np
def bisection_method(f, a, b, tol, max_iter):
    """
    Tìm nghiệm của phương trình f(x) = 0 trên đoạn [a, b] bằng phương pháp chia đôi
    Args:
        f (function): Hàm số cần tìm nghiệm
        a (float): Điểm đầu của đoạn
        b (float): Điểm cuối của đoạn
        tol (float): Tolerance - Sai số cho phép
        max_iter (int): Số lượng vòng lặp tối đa

    Returns:
        root (float): Nghiệm của phương trình
    """
    if a >= b:
        print("Điểm đầu phải nhỏ hơn điểm cuối!")
        return None
    if f(a) * f(b) >= 0:
        print("Không thể tìm nghiệm với phương pháp chia đôi")
        return None
    if f(a) == 0:
        return a
    if f(b) == 0:
        return b
    i = 0
    for _ in range(max_iter):
        c = (a + b) / 2  # Tính điểm giữa
        print("Lần lặp thứ ", i)
        data = np.array([[a, b, c]])
        headers = np.array(["an", "bn", "cn+1"])
        print(tabulate(data, headers=headers, tablefmt="grid", floatfmt=".15f"))
        # print(a,b,c)
        
        if abs(b - a) / 2 < tol or f(c) == 0:  # Kiểm tra điều kiện dừng
            return c
        
        # Cập nhật lại a hoặc b
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
        i += 1
    print("Không tìm được nghiệm trong số lần lặp tối đa!")
    return None
def main():
    """
    Hàm chính
    """
    # Nhập hàm số
    def f(x):
        return x**5 - 0.2*x + 15 #Nhập hàm số ở đây

    a = float(input("Nhập điểm đầu a: "))  # Điểm đầu
    b = float(input("Nhập điểm cuối b: "))  # Điểm cuối
    tol = float(input("Nhập sai số: "))  # Sai số cho phép
    max_iter = int(input("Nhập số lượng vòng lặp tối đa: "))  # Số lượng vòng lặp tối đa

    result = bisection_method(f, a, b, tol, max_iter)
    print("Nghiệm của phương trình là:", result)
main()