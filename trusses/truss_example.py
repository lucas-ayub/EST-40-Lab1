from node import Node
from bar import Bar
from truss import *

"""
TODO:
- Calculate the stress at all bars: check if the current implementation is correct for the method
- Implement a method in Truss class to insert at all nodes theirs displacements
"""

node_1 = Node(x = 0, y = 0, fy = 0, fx = 0, fixed_in_x = True, fixed_in_y = True)
node_2 = Node(x = 10, y = 0, fy = 0, fx = 0, fixed_in_x = False, fixed_in_y = True)
node_3 = Node(x = 5, y = 8, fx = 20e3, fy = 30e3, fixed_in_x = False, fixed_in_y = False)
node_4 = Node(x = 5, y = 4, fy = 0, fx = 0, fixed_in_x = False, fixed_in_y = False) 

list_of_nodes = [node_1, node_2, node_3, node_4]

E = 200e9
A = 900e-6

bar_1 = Bar(left_node = node_1, right_node = node_2,  E = E, A = A)
bar_2 = Bar(left_node = node_2, right_node = node_3, E = E, A = A)
bar_3 = Bar(left_node = node_3, right_node = node_1, E = E, A = A)
bar_4 = Bar(left_node = node_1, right_node = node_4, E = E, A = A)
bar_5 = Bar(left_node = node_2, right_node = node_4, E = E, A = A)
bar_6 = Bar(left_node = node_3, right_node = node_4, E = E, A = A)

list_of_bars = [bar_1, bar_2, bar_3, bar_4, bar_5, bar_6]

truss = Truss(list_of_nodes = list_of_nodes, list_of_bars = list_of_bars)

print(truss.getSolution())