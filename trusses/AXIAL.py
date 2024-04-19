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



node_1 = Node(x = 0, y = 0, fx = 0, fy = 0, fixed_in_x = True, fixed_in_y = False)
node_2 = Node(x = L1, y = 0, fx = -2*P, fy = 0, fixed_in_x = False, fixed_in_y = False)
node_3 = Node(x = L1 + L2, y = 0, fx = 0, fy = 0, fixed_in_x = True, fixed_in_y = False)

list_of_nodes = []
list_of_nodes.append(node_1)

n1 = 10
bars_between_1_2, nodes_between_1_2 = discretizateBar(node_1, node_2, q1, E1, A1, n1)

list_of_nodes += nodes_between_1_2
list_of_nodes.append(node_2)

list_of_bars = bars_between_1_2

n2 = 10
bars_between_2_3, nodes_between_2_3 = discretizateBar(node_2, node_3, q2, E2, A2, n2)

list_of_nodes += nodes_between_2_3
list_of_nodes.append(node_3)

list_of_bars += bars_between_2_3

# list_of_nodes = [node_1, node_2, node_3]

# bar_1 = Bar(left_node=node_1, right_node=node_2, E=E1, A=A1, q=q1)
# bar_2 = Bar(left_node=node_2, right_node=node_3, E=E2, A=A2, q=q2)

# list_of_bars = [bar_1, bar_2]

axial_bar = Structure(list_of_nodes = list_of_nodes, list_of_bars = list_of_bars)

print(f"\n\n {axial_bar.getSolution()}")
# axial_bar.plotStructure(displacement_scale = 1000)

# print(-2*P + q1*L1 + q2*L2)
