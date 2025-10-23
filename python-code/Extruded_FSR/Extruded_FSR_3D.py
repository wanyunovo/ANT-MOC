import re
import csv
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 读取数据文件
with open('Extruded_FSR', 'r') as file:
    data = file.readlines()

# 正则表达式模式
pattern = r'FSR Track Start\((\-?\d+\.\d+),(\-?\d+\.\d+),(\-?\d+\.\d+)\), End\((\-?\d+\.\d+),(\-?\d+\.\d+),(\-?\d+\.\d+)\)'

# 创建 CSV 文件并写入数据
with open('Extruded_FSR.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Start_X', 'Start_Y', 'Start_Z', 'End_X', 'End_Y', 'End_Z'])
    
    for line in data:
        match = re.search(pattern, line)
        if match:
            csvwriter.writerow(match.groups())

# 读取 CSV 文件
df = pd.read_csv('Extruded_FSR.csv')

# 创建3D图
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 绘制每个线段
for index, row in df.iterrows():
    ax.plot([row['Start_X'], row['End_X']], 
            [row['Start_Y'], row['End_Y']], 
            [row['Start_Z'], row['End_Z']])

# 设置轴标签
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# 设置标题和轴标签
plt.title('Extruded FSR Segment')

# 设置图片分辨率为2K以上
plt.savefig('Extruded_FSR_3D.png', dpi=300, bbox_inches='tight')