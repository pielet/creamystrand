import math

n_node = 4

with open('../assets/strands.obj', 'r') as f:
    lines = f.readlines()

v_lines = filter(lambda x: x.startswith('v'), lines)
l_lines = filter(lambda x: x.startswith('l'), lines)

verteces = list(map(lambda x: [float(s) for s in x.split()[1:]], v_lines))
idx = list(map(lambda x: [int(s) for s in x.split()[1:]], l_lines))

new_verteces = []
new_idx = []
node_count = 1
for i in range(len(idx)):
    strand_n_node = len(idx[i])
    strip = int(strand_n_node / (n_node - 1))
    start_idx = idx[i][0]
    new_strand_idx = []
    for j in range(n_node):
        new_verteces.append(verteces[start_idx + j * strip])
        new_strand_idx.append(node_count)
        node_count += 1
    new_idx.append(new_strand_idx)
    
with open('../assets/{}_node_strands.obj'.format(n_node), 'w') as f:
    for vtx in new_verteces:
        f.write('v {} {} {}\n'.format(*vtx))
    for line in new_idx:
        f.write(('l' + ' {}' * n_node + '\n').format(*line))
