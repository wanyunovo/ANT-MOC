import csv
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# # 源文件路径
# source_text_file = 'track_segment_in_hex'

# # 提取并保存 segment 数据
# segment_data = []
# segment_hex_data = []
# track_data = []

# # 正则表达式匹配规则
# segment_pattern = re.compile(r'segment start\(([^)]+)\), end\(([^)]+)\)')
# track_pattern = re.compile(r'2D Track Array id:(\d+), start\(([^)]+)\), end\(([^)]+)\)')

# with open(source_text_file, 'r') as file:
#     for line in file:
#         segment_match = segment_pattern.search(line)
#         track_match = track_pattern.search(line)
        
#         if segment_match:
#             start_coords, end_coords = segment_match.groups()
#             segment_data.append([start_coords, end_coords])
#         elif track_match:
#             track_id, start_coords, end_coords = track_match.groups()
#             track_data.append([track_id, start_coords, end_coords])

# # 保存 segment 数据到 CSV
# with open('segment_data_in_hex.csv', 'w', newline='') as csvfile:
#     csvwriter = csv.writer(csvfile)
#     csvwriter.writerow(['Start', 'End'])
#     csvwriter.writerows(segment_data)

# # 保存 track 数据到 CSV
# with open('track_data_in_hex.csv', 'w', newline='') as csvfile:
#     csvwriter = csv.writer(csvfile)
#     csvwriter.writerow(['ID', 'Start', 'End'])
#     csvwriter.writerows(track_data)

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

# 读取segment数据
segment_data = pd.read_csv('segment_data_in_hex.csv')

# 读取track数据
track_data = pd.read_csv('track_data_in_hex.csv')

# 创建一个新的绘图窗口
plt.figure()

size = 0.401  # 六边形的边长

# 计算中心位置
cx, cy = 0, 0

# 绘制六边形
ax = plt.gca()
ax.set_aspect('equal')

# 添加矩形边界
# 矩形边界的四条边
x_min, x_max = -3.0, 3.0
y_min, y_max = -3.5, 3.5

# 绘制4条边
ax.plot([x_min, x_max], [y_min, y_min], color='black') # 底边
ax.plot([x_min, x_max], [y_max, y_max], color='black') # 顶边
ax.plot([x_min, x_min], [y_min, y_max], color='black') # 左边
ax.plot([x_max, x_max], [y_min, y_max], color='black') # 右边

# 五圈六边形
radius = 4  # 半径为5的六边形网格
for q in range(-radius, radius + 1):
    for r in range(-radius, radius + 1):
        if -radius <= (q + r) <= radius:
            hex_center_x = size * 3/2 * q
            hex_center_y = size * np.sqrt(3) * (r + q/2)
            hex_points = hexagon((cx + hex_center_x, cy + hex_center_y), size)
            ax.plot(hex_points[:, 0], hex_points[:, 1], 'k-', linewidth=1, color = 'black')

# 绘制track线段
for i, row in track_data.iterrows():
    start_x, start_y = map(float, row['Start'].split(','))
    end_x, end_y = map(float, row['End'].split(','))
    plt.plot([start_x, end_x], [start_y, end_y], color='gray', linewidth=1, label=f'Track {i+1}')

# 绘制segment线段
colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'purple', 'orange', 'brown', 'pink']  # 自定义颜色列表
for i, row in segment_data.iterrows():
    start_x, start_y = map(float, row['Start'].split(','))
    end_x, end_y = map(float, row['End'].split(','))
    plt.plot([start_x, end_x], [start_y, end_y], color=colors[i % len(colors)], linewidth=1, label=f'Segment {i+1}')

# 设置等比例
plt.gca().set_aspect('equal', adjustable='box')

# 设置标题和轴标签
# plt.title('Track and Segments')
# plt.xlabel('X Coordinate')
# plt.ylabel('Y Coordinate')

# 设置图片分辨率为2K以上
plt.savefig('track_segment_in_hex.png', dpi=3000, bbox_inches='tight')