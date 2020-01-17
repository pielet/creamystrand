import argparse
import math
import random
import numpy as np


class HairGenerator:
    def __init__(self, output, radius, length, num_hair, num_vertex):
        self.output_file = output
        self.radius = radius
        self.length = length
        self.num_hair = num_hair
        self.num_vertex = num_vertex

        self.element_length = length / (num_vertex - 1)
        self.start_pos = []
        random.seed()

    def genHair(self):
        for i in range(self.num_hair):
            y = random.gauss(0, 1)
            while y <= 0:
                y = random.gauss(0, 1)
            x = random.gauss(0, 1)
            z = random.gauss(0, 1)
            vec = np.array([x, y, z])
            self.start_pos.append(vec / np.linalg.norm(vec))
    
    def output(self):
        with open(self.output_file, 'w') as f:
            for start_pos in self.start_pos:
                current_pos = start_pos * self.radius
                for i in range(self.num_vertex):
                    f.write("v {} {} {}\n".format(current_pos[0], current_pos[1], current_pos[2]))
                    current_pos += self.element_length * start_pos
            for i in range(self.num_hair):
                f.write("l ")
                f.write(" ".join(str(idx + 1) for idx in range(i * self.num_vertex, (i + 1)*self.num_vertex)))
                f.write('\n')


parser = argparse.ArgumentParser(description="Generate simple hair model")
parser.add_argument("-o", required=True, type=str, help="output file name (.obj file)")
parser.add_argument("-r", type=float, default=15, help="radius of base ball")
parser.add_argument("-l", type=float, default=20, help="hair length")
parser.add_argument("-ns", type=int, default=100, help="number of hair strand")
parser.add_argument("-nv", type=int, default=40, help="number of vertex in one hair")

args = parser.parse_args()

def main():
    hair_gen = HairGenerator(args.o, args.r, args.l, args.ns, args.nv)
    hair_gen.genHair()
    hair_gen.output()

if __name__ == "__main__":
    main()