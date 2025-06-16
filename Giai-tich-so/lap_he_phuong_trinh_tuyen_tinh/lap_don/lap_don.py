import numpy as np
from tabulate import tabulate

DECIMALS = 7

def read_input(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    n = int(lines[0])
    B = np.array([list(map(float, lines[i + 1].split())) for i in range(n)])
    d = np.array(list(map(float, lines[n + 1].split())))
    x0 = np.array(list(map(float, lines[n + 2].split())))
    TOL = float(lines[n + 3])
    N = int(lines[n + 4])
    return B, d, x0, TOL, N

def simple_iteration_method(B, d, x0, TOL, N, stop_condition="absolute"):
    q = np.linalg.norm(B, ord=np.inf)
    if q >= 1:
        raise ValueError(f"Không thỏa mãn điều kiện hội tụ: ||B|| = {q:.{DECIMALS}f} ≥ 1")
    print(f"Thỏa mãn điều kiện hội tụ: ||B|| = {q:.{DECIMALS}f} <= 1")

    tol = TOL * (1 - q) / q
    k = 1
    x_prev = x0.copy()
    x_first = None
    iteration_log = []

    while k <= N:
        x_new = B @ x_prev + d
        iteration_log.append((k, *x_new))

        if k == 1:
            x_first = x_new.copy()

        # Tính các loại sai số
        abs_err = np.linalg.norm(x_new - x_prev, ord=np.inf)
        rel_err = abs_err / np.linalg.norm(x_new, ord=np.inf)
        post_err = (q / (1 - q)) * abs_err
        pre_err = (q ** k / (1 - q)) * np.linalg.norm(x_first - x0, ord=np.inf)

        print(f"Lặp lần {k}")
        print(f"sai số tuyệt đối:{abs_err}")
        print(f"sai số tương đối:{rel_err}")
        print(f"sai số hậu nghiệm:{post_err}")
        print(f"sai số tiên nghiệm:{pre_err}")
        print ("----------------------")

        # Điều kiện dừng dựa trên loại sai số
        if stop_condition == "absolute" and abs_err <= tol:
            break
        elif stop_condition == "relative" and rel_err <= tol:
            break
        elif stop_condition == "post" and post_err <= tol:
            break
        elif stop_condition == "pre" and pre_err <= tol:
            break

        x_prev = x_new
        k += 1

    # Tính các loại sai số cuối cùng
    abs_err_final = np.linalg.norm(x_new - x_prev, ord=np.inf)
    rel_err_final = abs_err_final / np.linalg.norm(x_new, ord=np.inf)
    post_err_final = (q / (1 - q)) * abs_err_final
    pre_err_final = (q ** k / (1 - q)) * np.linalg.norm(x_first - x0, ord=np.inf)

    return x_new, k, iteration_log, abs_err_final, rel_err_final, post_err_final, pre_err_final

def main():
    B, d, x0, TOL, N = read_input("input.txt")

    stop_condition = input("Chọn điều kiện dừng (absolute/relative/post/pre): ").lower()

    try:
        x_final, k, logs, abs_err, rel_err, post_err, pre_err = simple_iteration_method(B, d, x0, TOL, N, stop_condition)

        print("📘 Quá trình lặp:")
        print(tabulate(logs, headers=["Lần lặp"] + [f"x{i + 1}" for i in range(len(x0))], floatfmt=f".{DECIMALS}f", tablefmt="fancy_grid"))

        print("\n🔎 Nghiệm gần đúng cuối cùng:", x_final)
        print("🔁 Số lần lặp:", k)

        print("\n📊 Sai số cuối cùng:")
        print(tabulate([
            ["Sai số tuyệt đối", abs_err],
            ["Sai số tương đối", rel_err],
            ["Sai số hậu nghiệm (CT 3.3)", post_err],
            ["Sai số tiên nghiệm (CT 3.2)", pre_err],
        ], headers=["Loại sai số", "Giá trị"], floatfmt=f".{DECIMALS}e", tablefmt="fancy_grid"))

    except ValueError as e:
        print(e)

if __name__ == "__main__":
    main()
