from node import *
from bar import *
from structure import *



node_1 = Node(x = 0, y = 0, fx = 0, fy = 0, fixed_in_x = True, fixed_in_y = False)
node_2 = Node(x = 1, y = 0, fx = 0, fy = 0, fixed_in_x = False, fixed_in_y = False)
node_3 = Node(x = 2, y = 0, fx = 150, fy = 0, fixed_in_x = False, fixed_in_y = False)
node_4 = Node(x = 3, y = 0, fx = 0, fy = 0, fixed_in_x = False, fixed_in_y = False)


list_of_nodes = [node_1, node_2, node_3, node_4]

E = 28000
A = 1
q = 45

bar_1 = Bar(left_node = node_1, right_node = node_2, q = 2*q, E = 2*E, A = A)
bar_2 = Bar(left_node = node_2, right_node = node_3, q = q, E = E, A = A)
bar_3 = Bar(left_node = node_3, right_node = node_4, q = q, E = 2*E, A = A)


list_of_bars = [bar_1, bar_2, bar_3]

structure = Structure(list_of_nodes, list_of_bars)
# print(structure.solve())    
print(structure.getStiffnessMatrix()) 
print(f"\n\n {structure.getSolution()}")
structure.plotStructure(displacement_scale = 1000)