import re
import csv
from collections import Counter

# 假设你的txt文件名为data.txt
input_filename = '/home/gengxiaoxu/data/Hex-CMFD-v0.1/ant-moc-develop/Log-Hex-Modify/output'
output_filename = 'surfaces_count_sorted.csv'

# 使用正则表达式匹配surface数据
surface_pattern = re.compile(r'surface\[\d+\]:(\d+)')

# 初始化计数器
surface_counts = Counter()

# 读取txt文件
with open(input_filename, 'r') as file:
    content = file.readlines()

# 计算每个surface的出现次数
for line in content:
    match = surface_pattern.search(line)
    if match:
        surface_counts[match.group(1)] += 1

# 对surface_counts进行排序，键转换为整数以确保数字排序而非字符串排序
sorted_surfaces = sorted(surface_counts.items(), key=lambda x: int(x[0]))

# 写入CSV文件，已排序
with open(output_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    # 写入标题行
    writer.writerow(['surface', 'count'])
    # 按surface编号写入数据
    for surface, count in sorted_surfaces:
        writer.writerow([surface, count])