from math import sqrt

vertices = []
hairs = []
first_pos = []

def distance(p1, p2):
    dis = 0
    for i in range(3):
        dis += (p1[i] - p2[i]) * (p1[i] - p2[i])
    return sqrt(dis)

with open('C:/Users/v-shiyji/Desktop/creamystrand/assets/hairbundle.obj', 'r') as f:
    for line in f:
        items = line.split(' ')
        if items[0] == 'v':
            vertices.append([float(x) for x in items[1:]])
        elif items[0] == 'l':
            hairs.append([int(x) for x in items[1:]])
            first_pos.append(vertices[int(items[-1]) - 1])

num_vert = len(vertices)
num_hair = len(hairs)

center_pos = [0, 0, 0]
for pos in first_pos:
    for i in range(3):
        center_pos[i] += pos[i]
for i in range(3):
    center_pos[i] /= num_hair

new_hairs = []
new_vertices = []
globle_idx = 1
for i, hair in enumerate(hairs):
    if distance(center_pos, first_pos[i]) < 0.05:
        idx_list = []
        v_list = []
        for idx in hair:
            v_list.append(vertices[idx - 1])
            idx_list.append(globle_idx)
            globle_idx += 1
        new_hairs.append(idx_list)
        new_vertices.extend(v_list)

with open('C:/Users/v-shiyji/Desktop/creamystrand/assets/strands.obj', 'w') as f:
    for v in new_vertices:
        f.write('v {} {} {}\n'.format(*v))
    for h in new_hairs:
        f.write('l')
        for idx in h:
            f.write(' {}'.format(idx))
        f.write('\n')
        