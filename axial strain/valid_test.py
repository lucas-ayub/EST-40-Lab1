from node import *
from bar import *
from structure import *

node_1 = Node(x = 0, fe = 0, fixed = True)
node_2 = Node(x = 1, fe = 0, fixed = False)
node_3 = Node(x = 2, fe = 150, fixed = False)
node_4 = Node(x = 3, fe = 0, fixed = False)


nodes = [node_1, node_2, node_3, node_4]

bar_1 = Bar(left_node = node_1, right_node = node_2, q = 45, E = 2*28000, A = 1)
bar_2 = Bar(left_node = node_2, right_node = node_3, q = 22.5, E = 28000, A = 1)
bar_3 = Bar(left_node = node_3, right_node = node_4, q = 22.5, E = 2*28000, A = 1)


bars = [bar_1, bar_2, bar_3]

structure = Structure(nodes, bars)
print(structure.num_nodes)
print(structure.symbols)
print(structure.K)
print(structure.solve())    