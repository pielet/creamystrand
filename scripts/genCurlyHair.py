from math import sin, cos, pi

n_node = 40
radius = 0.5
n_circle = 8    # number of node per circle
height = 10  # total height

vert = []
delta_theta = 2 * pi / n_circle
delta_height = height / n_node
for i in range(n_node):
    x = sin(i * delta_theta)
    y = cos(i * delta_theta)
    z = i * delta_height
    vert.append([x, y, z])

with open("curly.txt", "w") as f:
    for pos in vert:
        f.write("<particle x=\"{} {} {}\" v=\"0.0 0.0 0.0\" fixed=\"0\"/>\n".format(*pos))