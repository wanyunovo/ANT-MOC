import matplotlib.pyplot as plt
import numpy as np

def hexagon(center, size):
    """生成一个正六边形的顶点
    Args:
        center (tuple): 六边形中心的坐标
        size (float): 六边形的边长
    Returns:
        numpy.ndarray: 六边形顶点的坐标
    """
    points = []
    for i in range(6):
        angle = 2 * np.pi / 6 * i
        x_i = center[0] + size * np.cos(angle)
        y_i = center[1] + size * np.sin(angle)
        points.append([x_i, y_i])
    return np.array(points)

# 设置绘图大小
plt.figure(figsize=(10, 10))
size = 0.401  # 六边形的边长

# 计算中心位置
cx, cy = 0, 0

# 绘制六边形
ax = plt.gca()
ax.set_aspect('equal')

# 五圈六边形
radius = 4  # 半径为5的六边形网格
for q in range(-radius, radius + 1):
    for r in range(-radius, radius + 1):
        if -radius <= (q + r) <= radius:
            hex_center_x = size * 3/2 * q
            hex_center_y = size * np.sqrt(3) * (r + q/2)
            hex_points = hexagon((cx + hex_center_x, cy + hex_center_y), size)
            ax.plot(hex_points[:, 0], hex_points[:, 1], 'k-')

# 隐藏坐标轴
plt.axis('off')

# 设置图片分辨率为2K以上
plt.savefig('Hex_2D.png', dpi=300, bbox_inches='tight')