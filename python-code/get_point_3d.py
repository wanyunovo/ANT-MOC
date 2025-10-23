import re
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter


# 假设你的txt文件名为data.txt
input_filename = '/home/gengxiaoxu/data/Hex-CMFD-v0.1/ant-moc-develop/Log-Hex-Modify/output'
output_filename = 'cefr_xyz.csv'

"""
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
"""

# 前面的代码不变，在此基础上添加动画的更新函数
def update(num, data, ax):
    ax.view_init(azim=num)
    return ax,

# 读取csv数据
data = pd.read_csv(output_filename)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')  # 注意这里变化，添加了三维投影

def random_color():
    return (random.random(), random.random(), random.random())

# 修改这一部分，以便绘制三维线段
for i in range(0, len(data), 2):
    ax.plot(data['x'].iloc[i:i+2], data['y'].iloc[i:i+2], data['z'].iloc[i:i+2], color=random_color(), linewidth=0.7) 

# 设置坐标轴标签
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

# 创建FuncAnimation对象，用于生成动画，这里设置动画每次旋转度数等于1
ani = FuncAnimation(fig, update, frames=np.arange(0, 360, 1), fargs=(data, ax))

# 设置动画的保存参数，定义动画的帧率和比特率等
writer = FFMpegWriter(fps=20, metadata=dict(artist='Me'), bitrate=1800)

# 保存这个动画到文件中（如果环境中没有安装ffmpeg，这一步会失败）
ani.save('cefr_xyz.mp4', writer=writer)

# 设置图片分辨率为2K以上
plt.savefig('cefr_xyz.png', dpi=300, bbox_inches='tight')

# 如果要在屏幕上显示图形，取消注释下面这行
plt.show()