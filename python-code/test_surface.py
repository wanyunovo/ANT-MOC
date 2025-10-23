NUM_FACES = 6
NUM_EDGES = 12
NUM_VERTICES = 8
NUM_SURFACES = NUM_FACES + NUM_EDGES + NUM_VERTICES

def calculate_direction(surface):
    direction = [0, 0, 0]
    
    if surface < NUM_FACES:
        ind = surface % 3
        dir = 2 * (surface // 3) - 1
        direction[ind] = dir
    elif surface < NUM_FACES + NUM_EDGES:
        surface -= NUM_FACES
        group = surface // 4
        skipped = 2 - group
        surface %= 4
        ind = [0, 0]
        ind[0] = surface % 2
        ind[1] = (surface - ind[0]) // 2
        n = 0
        for i in range(3):
            if i != skipped:
                direction[i] = 2 * ind[n] - 1
                n += 1
    elif surface < NUM_SURFACES:
        surface -= NUM_FACES + NUM_EDGES
        direction[0] = 2 * (surface // 4) - 1
        direction[1] = 2 * ((surface // 2) % 2) - 1
        direction[2] = 2 * (surface % 2) - 1
    else:
        raise ValueError(f"Invalid surface ID {surface}")

    return direction

# 测试函数
def test_calculate_direction():
    # 测试不同的surface值以确保函数的正确性
    test_surfaces = [0, NUM_FACES - 1, NUM_FACES, NUM_FACES + NUM_EDGES - 1, NUM_SURFACES - 1]
    for surf in test_surfaces:
        try:
            print(f"Surface {surf}: Direction {calculate_direction(surf)}")
        except ValueError as e:
            print(e)

# 运行测试
#test_calculate_direction()


def convert_direction_to_surface(direction):
    surface = 0
    # 计算方向向量中非零元素的数量
    num_crossings = abs(direction[0]) + abs(direction[1]) + abs(direction[2])
    
    if num_crossings == 1:
        for i in range(3):
            present = abs(direction[i])
            fwd = (direction[i] + 1) // 2
            surface += present * (3 * fwd + i)
    elif num_crossings == 2:
        surface += NUM_FACES
        ind1 = ind2 = 0
        if direction[0] == 0:
            ind1, ind2 = direction[1], direction[2]
            surface += 8
        elif direction[1] == 0:
            ind1, ind2 = direction[0], direction[2]
            surface += 4
        elif direction[2] == 0:
            ind1, ind2 = direction[0], direction[1]
        ind1 = (ind1 + 1) // 2
        ind2 = (ind2 + 1) // 2
        surface += 2 * ind2 + ind1
    elif num_crossings == 3:
        surface += NUM_FACES + NUM_EDGES
        fwd = [(direction[i] + 1) // 2 for i in range(3)]
        surface += 4 * fwd[0] + 2 * fwd[1] + fwd[2]
    else:
        raise ValueError("Invalid number of surface crossings")  # 代替ferror
    
    return surface


# 假设direction已经定义
direction = [1, 1, 1]  # 示例方向向量

# 初始化remainder_direction和partial_direction
remainder_direction = [0, 0, 0]
partial_direction = [0, 0, 0]

# 存储最终结果的变量
remainder_surfaces = [0, 0, 0]
partial_surfaces = [0, 0, 0]

for i in range(3):
    for j in range(3):
        if i == j:
            remainder_direction[j] = 0
            partial_direction[j] = direction[j]
        else:
            remainder_direction[j] = direction[j]
            partial_direction[j] = 0

    remainder_surfaces[i] = convert_direction_to_surface(remainder_direction)
    partial_surfaces[i] = convert_direction_to_surface(partial_direction)

# 打印结果进行检查
print("Remainder Surfaces: ", remainder_surfaces)
print("Partial Surfaces: ", partial_surfaces)