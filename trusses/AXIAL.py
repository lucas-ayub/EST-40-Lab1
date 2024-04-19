from node import *
from bar import *
from structure import *


E1 = 0.8e9
A1 = 1000e-6
E2 = 2.7e9
A2 = 1875e-6
P = 25e3
q = 10e3
L1 = 150e-3
L2 = 262.5e-3

q1 = -2*q
q2 = q
F1 = 0
F2 = -2*P
F3 = 0

node_1 = Node(x = 0, y = 0, fx = F1, fy = 0, fixed_in_x = True, fixed_in_y = False)
node_2 = Node(x = L1, y = 0, fx = F2, fy = 0, fixed_in_x = False, fixed_in_y = False)
node_3 = Node(x = L1 + L2, y = 0, fx = F3, fy = 0, fixed_in_x = True, fixed_in_y = False)

list_of_nodes = []

list_of_nodes.append(node_1)

# Creating elements between node_1 and node_2
n1 = 2
bars_between_1_2, nodes_between_1_2 = discretizateBar(node_1, node_2, q1, E1, A1, n1)

list_of_nodes += nodes_between_1_2
list_of_nodes.append(node_2)

list_of_bars = bars_between_1_2

# Creating elements between node_2 and node_3
n2 = 4
bars_between_2_3, nodes_between_2_3 = discretizateBar(node_2, node_3, q2, E2, A2, n2)

list_of_nodes += nodes_between_2_3
list_of_nodes.append(node_3)

list_of_bars += bars_between_2_3

axial_bar = Structure(list_of_nodes = list_of_nodes, list_of_bars = list_of_bars)

displacements_and_forces_of_nodes = axial_bar.getSolution()
bars_forces_and_stresses = axial_bar.getBarsStressesAndNormals()

print(f'System solution: \n{displacements_and_forces_of_nodes}')

print(f'\n Normal forces and stresses: \n {bars_forces_and_stresses}')

axial_bar.plotStructure(displacement_scale = 1)