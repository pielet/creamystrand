import argparse
import xml.etree.ElementTree as ET
import numpy as np

parser = argparse.ArgumentParser(description="Generate training data from .xml and .ply")
# parser.add_argument("--ply_path", type=str, required=True, help="path to rods.ply")
parser.add_argument("--xml_file", type=str, required=True, help="xml file")

args = parser.parse_args()

# parse xml (parameters, Y)
tree = ET.parse(args.xml_file)
root = tree.getroot()

description = root.find("description").get("text")
duration = float(root.find("duration").get("time"))
dt = float(root.find("integrator").get("dt"))
n_frame = int(duration / dt)

for particle in root.findall("particle"):
    if int(particle.get("fixed")) == 3:
        pos = particle.get("x")
        fixed_point = np.array([float(x) for x in pos.split()])
        fixed_point.resize(1, 3)

Y = np.tile(fixed_point, [401, 1])
print(Y.shape)
np.save("Y.npy", Y)

# parse .ply (X)
