import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

# 读取 CSV 文件
df = pd.read_csv('Extruded_FSR.csv')

# 创建图形
fig, ax = plt.subplots()

# 设置点的颜色和尺寸
num_points = len(df)
colors = np.random.rand(num_points)  # 为每个点生成不同的颜色
sizes = np.full(num_points, 10)  # 设置所有点的尺寸为10

# 添加矩形边界
# 矩形边界的四条边
x_min, x_max = -3.0, 3.0
y_min, y_max = -3.5, 3.5

# 绘制4条边
ax.plot([x_min, x_max], [y_min, y_min], color='black') # 底边
ax.plot([x_min, x_max], [y_max, y_max], color='black') # 顶边
ax.plot([x_min, x_min], [y_min, y_max], color='black') # 左边
ax.plot([x_max, x_max], [y_min, y_max], color='black') # 右边

# 绘制每个起点
scatter = ax.scatter(df['Start_X'], df['Start_Y'], c=colors, s=sizes, cmap='viridis', alpha=0.6, marker='.')

# 设置轴标签
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')

# 设置标题
plt.title('Start Points of Extruded FSR Segment (XY Plane)')

# 确保XY等比例
ax.set_aspect('equal', adjustable='box')

# 设置图片分辨率为2K以上
plt.savefig('Extruded_FSR_Start_2D.png', dpi=300, bbox_inches='tight')

# 显示图形（可选）
plt.show()