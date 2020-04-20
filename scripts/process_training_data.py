import os
import argparse
import numpy as np
import xml.etree.ElementTree as ET

parser = argparse.ArgumentParser(description="Genarate X and Y from .ply and .xml file")

parser.add_argument("--output_rate", type=int, required=True, help="output rate of .ply file")
parser.add_argument("--xml", type=str, required=True, help="path to xml file")
parser.add_argument("--ply", type=str, required=True, help="path to ply files")

args = parser.parse_args()

###################### read ply ############################
sim_data = []
ply_path = args.ply

for filename in os.listdir(ply_path):
    if filename.split('_')[0] == "rods":
        hairs = []
        verts = []
        n_vert = n_hair = 0
        with open(os.path.join(ply_path, filename), 'r') as f:
            for i, line in enumerate(f):
                if i == 3:
                    n_vert = int(line.split()[2])
                if i == 16:
                    n_hair = int(line.split()[2])
                if 19 <= i < 19 + n_vert:
                    nums = [float(x) for x in line.split()]
                    verts.append(nums[:3])
                if 19 + n_vert <= i < 19 + n_vert + n_hair:
                    nums = [int(x) for x in line.split()]
                    n_hair_vert = nums[0]
                    hair = []
                    for vi in nums[1:-1]:
                        hair.extend(verts[vi])
                    hairs.append(hair)
        sim_data.append(hairs)

np.save("X.npy", np.array(sim_data))

############################ read xml ########################
tree = ET.parse(args.xml)
root = tree.getroot()

duration = float(root.find("duration").get("time"))
dt = float(root.find("integrator").get("dt")) * args.output_rate
n_frame = len(sim_data)

# intialize group
n_group = 0
for particle in root.findall("particle"):
    gid = particle.get("group", "0")
    if int(gid) > n_group:
        n_group += 1
for hairobj in root.findall("hairobj"):
    oid = hairobj.get("group", "0")
    if int(oid) > n_group:
        n_group += 1

Y = np.zeros([n_frame, n_group * 6])    # [position, rotation] for every group

def process_translate(direction, start, end, position):
    for i in range(n_frame):
        t = duration * i / (n_frame - 1)
        if start <= t < end:
            position[i, :] += (t - start) / (end - start) * direction
        elif t >= end:
            position[i, :] += direction

def process_rotate(axis, angle, start, end, rotation):
    for i in range(n_frame):
        t = duration * i / (n_frame - 1)
        if start <= t < end:
            rotation[i, axis] += (t - start) / (end - start) * angle
        elif t >= end:
            rotation[i, axis] += angle

for script in root.findall("script"):
    start = float(script.get("start"))
    end = float(script.get("end"))
    group = int(script.get("group"))
    if script.get("type") == "translate":
        direction = np.zeros([1, 3])
        direction[0] = script.get("x")
        direction[1] = script.get("y")
        direction[2] = script.get("z")
        process_translate(direction, start, end, Y[:, 6 * group: 6 * group + 3])
    if script.get("type") == "rotate":
        if script.get("x") == "1":
            axis = 0
        if script.get("y") == "1":
            axis = 1
        if script.get("z") == "1":
            axis = 2
        angle = float(script.get("w"))
        process_rotate(axis, angle, start, end, Y[:, 6 * group + 3: 6 * group + 6])

np.save("Y.npy", Y)
