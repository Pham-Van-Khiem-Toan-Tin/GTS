from sympy import symbols, solveset, diff, Interval, Union, S

def check_sign_changes(f, a_0, b_0):
    x = symbols("x")
    a = a_0
    b = b_0
    interval = Interval(a, b)
    
    # Tính đạo hàm bậc 1 và bậc 2
    f_prime = diff(f, x)    # f'(x)
    f_double_prime = diff(f, x, 2)  # f''(x)
    
    # Tìm các nghiệm của f'(x) = 0 và f''(x) = 0 trong [a, b]
    critical_points = solveset(f_prime, x, interval)  # Nghiệm của f'(x) = 0
    inflection_points = solveset(f_double_prime, x, interval)  # Nghiệm của f''(x) = 0
    
    # Gộp các điểm và thêm đầu mút
    all_points = Union(critical_points, inflection_points)
    all_points = [a] + sorted([p for p in all_points if p.is_real]) + [b]
    
    # Kiểm tra dấu của f'(x) trên từng khoảng
    f_prime_signs = []
    for i in range(len(all_points) - 1):
        test_point = (all_points[i] + all_points[i+1]) / 2  # Điểm giữa khoảng
        try:
            value = f_prime.subs(x, test_point)
            if value > 0:
                f_prime_signs.append(1)
            elif value < 0:
                f_prime_signs.append(-1)
            else:
                f_prime_signs.append(0)
        except:
            return "Không thể xác định dấu của f'(x)"
    
    # Kiểm tra dấu của f''(x) trên từng khoảng
    f_double_prime_signs = []
    for i in range(len(all_points) - 1):
        test_point = (all_points[i] + all_points[i+1]) / 2
        try:
            value = f_double_prime.subs(x, test_point)
            if value > 0:
                f_double_prime_signs.append(1)
            elif value < 0:
                f_double_prime_signs.append(-1)
            else:
                f_double_prime_signs.append(0)
        except:
            return "Không thể xác định dấu của f''(x)"
    
    # Kết luận
    f_prime_changes = any(f_prime_signs[i] * f_prime_signs[i+1] < 0 
                          for i in range(len(f_prime_signs)-1))
    f_double_prime_changes = any(f_double_prime_signs[i] * f_double_prime_signs[i+1] < 0 
                                 for i in range(len(f_double_prime_signs)-1))
    
    result = {
        "f'(x) đổi dấu": f_prime_changes,
        "f''(x) đổi dấu": f_double_prime_changes,
        "f'(x) signs": f_prime_signs,
        "f''(x) signs": f_double_prime_signs
    }
    return result

# Test
x = symbols("x")

# Ví dụ 1: f(x) = x^2 trên [0, 1]
print("f(x) = x^2 trên [0, 1]:")
print(check_sign_changes(x**2, 0, 1))

# Ví dụ 2: f(x) = x^3 trên [-1, 1]
print("\nf(x) = x^3 trên [-1, 1]:")
print(check_sign_changes(x**3, -1, 1))

# Ví dụ 3: f(x) = x^4 - 2x^2 trên [-2, 2]
print("\nf(x) = x^4 - 2x^2 trên [-2, 2]:")
print(check_sign_changes(x**4 - 2*x**2, -2, 2))