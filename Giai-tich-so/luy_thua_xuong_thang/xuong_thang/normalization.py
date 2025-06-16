import numpy as np

# ✅ Chuẩn hóa chuẩn 2 (L2-norm): độ dài vector = 1
# Phù hợp với hầu hết bài toán trị riêng, PCA, giải tích số
def normalize_L2(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else v

# ✅ Chuẩn hóa chuẩn 1 (L1-norm): tổng trị tuyệt đối các phần tử = 1
# Dùng khi xử lý véc-tơ thưa, bài toán sparse hoặc khi cần tổng cố định
def normalize_L1(v):
    norm = np.sum(np.abs(v))
    return v / norm if norm != 0 else v

# ✅ Chuẩn hóa theo phần tử đầu tiên: v[0] = 1
# Dễ đọc, dễ so sánh, thường dùng trong hiển thị báo cáo hoặc trình bày
def normalize_first_element(v):
    first = v[0][0] if isinstance(v[0], np.ndarray) else v[0]
    return v / first if first != 0 else v

# ✅ Chuẩn hóa để phần tử lớn nhất (về trị tuyệt đối) là 1, và dương
# Dùng để giữ tỉ lệ giữa các phần tử nhưng đảm bảo dấu thống nhất
def normalize_with_positive_max(v):
    idx = np.argmax(np.abs(v))
    sign = np.sign(v[idx])
    return v * sign / np.max(np.abs(v))

# ✅ Không chuẩn hóa: giữ nguyên véc-tơ đầu vào
# Dùng khi chỉ cần đúng phương (direction), không cần độ dài cố định
def no_normalization(v):
    return v
