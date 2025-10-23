import re
import csv

# 读取文件内容
with open('FSR_CMFD', 'r') as file:
    data = file.readlines()

# 解析数据
entries = []
pattern = re.compile(r"FSR id: (\d+), Centroid:\((-?\d+\.\d+),(-?\d+\.\d+),(-?\d+\.\d+)\), CMFD Cell:(\d+)")
for line in data:
    match = pattern.search(line)
    if match:
        fsr_id = int(match.group(1))
        centroid = f"{match.group(2)},{match.group(3)},{match.group(4)}"
        cmfd_cell = int(match.group(5))
        entries.append((fsr_id, centroid, cmfd_cell))

# 按照CMFD Cell排序
entries.sort(key=lambda entry: entry[2])

# 写入CSV文件
with open('FSR_CMFD_sorted.csv', 'w', newline='') as csvfile:
    fieldnames = ['CMFD Cell', 'FSR id', 'Centroid']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for entry in entries:
        writer.writerow({'CMFD Cell': entry[2], 'FSR id': entry[0], 'Centroid': entry[1]})

print("数据已成功存储到 FSR_CMFD_sorted.csv 文件中。")