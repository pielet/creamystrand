from random import random
from math import sqrt

def gen_random_dir(length, last_dir):
    dir = []
    while True:
        dir = []
        for i in range(3):
            dir.append(random() - 0.5)
        len_dir = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2])
        dir = list(map(lambda x: x * length / len_dir, dir))
        cos_angle = dir[0] * last_dir[0] + dir[1] * last_dir[1] + dir[2] * last_dir[2]
        if cos_angle > 0.8:
            break
    return dir


n_node = 40
length = 0.5

point_list = []
point_list.append([-1, 1, 1])
for i in range(n_node - 1):
    current_p = point_list[-1]
    dir = gen_random_dir(length, current_p)
    next_p = []
    for j in range(3):
        next_p.append(current_p[j] + dir[j])
    point_list.append(next_p)

with open("./random_init.ply", 'w') as f:
    f.write('''ply
format ascii 1.0
comment created by ADONIS
element vertex {}
property float x
property float y
property float z
property float theta
property int segment # strand index
property float ra
property float rb
property float ha
property float hb
property int group
property int actual
property float c0
property float c1
property float c2
property float c3
property float c4
element face 0
property list int int vertex_indices
end_header 
'''.format(n_node))
    for i in range(n_node):
        f.write("{} {} {} 0 0 0.056 0.005 0.005 0.005 1 1 1 1 1 1 1\n".format(*point_list[i]))
