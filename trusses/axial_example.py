from node import *
from bar import *
from truss import *


E1 = 0.8e9
A1 = 1000e-6
E2 = 2.7e9
A2 = 1875e-6
P = 25e3
q = 10e3
L1 = 150e-3
L2 = 262.5e-3

node_1 = Node(x = 0, y = 0, fx = 0, fy = 0, fixed_in_x = True, fixed_in_y = False)
node_2 = Node(x = L1, y = 0, fx = -2*P, fy = 0, fixed_in_x = False, fixed_in_y = False)
node_3 = Node(x = L2, y = 0, fx = 0, fy = 0, fixed_in_x = True, fixed_in_y = False)

list_of_nodes = [node_1, node_2, node_3]

bar_1 = Bar(left_node = node_1, right_node = node_2, q = -2*q, E = E1, A = A1)
bar_2 = Bar(left_node = node_2, right_node = node_3, q = q, E = E2, A = A2)

list_of_bars = [bar_1, bar_2]

truss = Truss(list_of_nodes = list_of_nodes, list_of_bars = list_of_bars)

print(f"\n\n {truss.getSolution()}")
truss.plot_truss(displacement_scale = 1000)
