from node import Node
from bar import Bar
from truss import Truss

node_1 = Node(x = 0, y = 0, fy = 0, fx = 0, fixed_in_x = True, fixed_in_y = True)
node_2 = Node(x = 40e-2, y = -30e-2, fy = -20e3, fx = 0, fixed_in_x = False, fixed_in_y = False)
node_3 = Node(x = 60e-2, y = 0, fy = 0, fx = 0, fixed_in_x = True, fixed_in_y = True)

list_of_nodes = [node_1, node_2, node_3]

bar_1 = Bar(left_node = node_1, right_node = node_2,  E = 1000e3, A = 1)
bar_2 = Bar(left_node = node_2, right_node = node_3,  E = 1000e3, A = 1)

list_of_bars = [bar_1, bar_2]

truss = Truss(list_of_nodes = list_of_nodes, list_of_bars = list_of_bars)

print(f"R: \n {truss.r}")
print(f"u: \n {truss.u}")

print(truss.solve())

print(f"u: \n {truss.u}")