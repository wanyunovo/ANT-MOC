import csv
import re
import pandas as pd
import matplotlib.pyplot as plt

# 源文件路径
source_text_file = 'track_segment_without_hex'

# 提取并保存 segment 数据
segment_data = []
segment_hex_data = []
track_data = []

# 正则表达式匹配规则
segment_pattern = re.compile(r'segment not in CMFD Cell start\(([^)]+)\), end\(([^)]+)\)')
track_pattern = re.compile(r'2D Track Array id:(\d+), start\(([^)]+)\), end\(([^)]+)\)')
segment_hex = re.compile(r'segment not in Hex CMFD Cell start\(([^)]+)\), end\(([^)]+)\)')

with open(source_text_file, 'r') as file:
    for line in file:
        segment_match = segment_pattern.search(line)
        track_match = track_pattern.search(line)
        segment_hex_match = segment_hex.search(line)
        
        if segment_match:
            start_coords, end_coords = segment_match.groups()
            segment_data.append([start_coords, end_coords])
        elif track_match:
            track_id, start_coords, end_coords = track_match.groups()
            track_data.append([track_id, start_coords, end_coords])
        elif segment_hex_match:
            start_coords, end_coords = segment_hex_match.groups()
            segment_hex_data.append([start_coords, end_coords])

# 保存 segment 数据到 CSV
with open('segment_data_without_hex.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Start', 'End'])
    csvwriter.writerows(segment_data)

# 保存 track 数据到 CSV
with open('track_data_without_hex.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['ID', 'Start', 'End'])
    csvwriter.writerows(track_data)

with open('segment_hex_data_without_hex.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Start', 'End'])
    csvwriter.writerows(segment_hex_data)


# 读取segment数据
segment_data = pd.read_csv('segment_data_without_hex.csv')

# 读取segment数据
segment_hex_data = pd.read_csv('segment_hex_data_without_hex.csv')

# 读取track数据
track_data = pd.read_csv('track_data_without_hex.csv')

# 创建一个新的绘图窗口
plt.figure()

ax = plt.subplots()

# 绘制track线段
for i, row in track_data.iterrows():
    start_x, start_y = map(float, row['Start'].split(','))
    end_x, end_y = map(float, row['End'].split(','))
    plt.plot([start_x, end_x], [start_y, end_y], color='black', linewidth=1, label=f'Track {i+1}')

# 绘制segment线段
colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'purple', 'orange', 'brown', 'pink']  # 自定义颜色列表
for i, row in segment_data.iterrows():
    start_x, start_y = map(float, row['Start'].split(','))
    end_x, end_y = map(float, row['End'].split(','))
    plt.plot([start_x, end_x], [start_y, end_y], color=colors[i % len(colors)], linewidth=1, label=f'Segment {i+1}')

# 绘制segment_hex线段
colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'purple', 'orange', 'brown', 'pink']  # 自定义颜色列表
for i, row in segment_hex_data.iterrows():
    start_x, start_y = map(float, row['Start'].split(','))
    end_x, end_y = map(float, row['End'].split(','))
    plt.plot([start_x, end_x], [start_y, end_y], color=colors[i % len(colors)], linewidth=1, label=f'Segment {i+1}')

# 设置等比例
plt.gca().set_aspect('equal', adjustable='box')

# 设置标题和轴标签
plt.title('Track and Segments')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')

# 设置图片分辨率为2K以上
plt.savefig('track_segment_without_hex.png', dpi=300, bbox_inches='tight')