def cell_in_xy_boundary_with_stencil(cell, stencil_index):
    inboundary = False

    num_x, num_y, num_z, num_r = 9, 9, 10, 5
    x = (cell % (num_x * num_y)) % num_x
    y = (cell % (num_x * num_y)) // num_x
    z = cell // (num_x * num_y)

    length = 2 * (num_r - 1)
    mid = num_r - 1
    high = (3 * num_r) - 3

    if x == 0 and y >= mid:  # 左边界，包含上下两端点
        if y == length and stencil_index not in (2, 5, 6):
            inboundary = True
        elif y == mid and stencil_index not in (0, 2, 5):
            inboundary = True
        elif mid < y < length and stencil_index not in (2, 5):
            inboundary = True
    elif x == length and y <= mid:  # 右边界，包含上下两端点
        if y == 0 and stencil_index not in (0, 1, 4):
            inboundary = True
        elif y == mid and stencil_index not in (1, 4, 6):
            inboundary = True
        elif 0 < y < mid and stencil_index not in (1, 4):
            inboundary = True
    elif y == 0 and mid <= x < length:  # 下边界，只包含左边界端点
        if x == mid and stencil_index not in (0, 1, 2):
            inboundary = True
        elif x > mid and stencil_index not in (0, 1):
            inboundary = True
    elif y == length and 0 < x <= mid:  # 上边界，只包含右边界端点
        if x == mid and stencil_index not in (4, 5, 6):
            inboundary = True
        elif x < mid and stencil_index not in (5, 6):
            inboundary = True
    elif x < mid and y < mid and (x + y == mid):  # 左下边界，不含端点
        if stencil_index not in (0, 2):
            inboundary = True
    elif x > mid and y > mid and (x + y == high):  # 右上边界，不含端点
        if stencil_index not in (4, 6):
            inboundary = True

    return inboundary

# 测试案例
test_cases = [
    (4, 4),
]
    
# 顺序执行测试案例
for cell, stencil_index in test_cases:
    inboundary = cell_in_xy_boundary_with_stencil(cell, stencil_index)
    print(f"Cell: {cell}, Stencil Index: {stencil_index} -> In Boundary: {inboundary}")