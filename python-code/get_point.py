import re
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random


# 假设你的txt文件名为data.txt
input_filename = '/home/gengxiaoxu/data/Hex-CMFD-v0.1/ant-moc-develop/Log-Hex-Modify/output'
output_filename = 'c5g7_xyz.csv'

# 使用正则表达式匹配point数据
point_pattern = re.compile(r'point:\((.*?),(.*?),(.*?)\)')

# 读取txt文件
with open(input_filename, 'r') as file:
    content = file.read()

# 使用正则表达式查找所有的point
points = point_pattern.findall(content)

# 写入CSV文件
with open(output_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    # 写入标题行
    writer.writerow(['x', 'y', 'z'])
    # 写入point数据
    for point in points:
        writer.writerow(point)

# 读取csv数据
data = pd.read_csv(output_filename)

fig, ax = plt.subplots()

# 随机颜色生成函数
def random_color():
    return (random.random(), random.random(), random.random())

# 使用DataFrame的每两行创建一个线段
for i in range(0, len(data), 2):
    plt.plot(data['x'].iloc[i:i+2], data['y'].iloc[i:i+2], data['z'].iloc[i:i+2], color=random_color(), linewidth=0.7) 

# 设置xy轴的尺度比例一致
ax.set_aspect('equal', adjustable='datalim')

# 设置图片分辨率为2K以上
plt.savefig('cefr_xyz.png', dpi=300, bbox_inches='tight')

# 如果要在屏幕上显示图形，可以取消注释下面这行
# plt.show()