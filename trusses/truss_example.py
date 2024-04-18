from node import Node
from bar import Bar
from structure import *

"""
TODO:
- Calculate the stress at all bars: check if the current implementation is correct for the method
"""

node_1 = Node(x = 0, y = 0, fy = 0, fx = 0, fixed_in_x = True, fixed_in_y = True)
node_2 = Node(x = 10, y = 0, fy = 0, fx = 0, fixed_in_x = False, fixed_in_y = True)
node_3 = Node(x = 5, y = 8, fx = 20e3, fy = 30e3, fixed_in_x = False, fixed_in_y = False)
node_4 = Node(x = 5, y = 4, fy = 0, fx = 0, fixed_in_x = False, fixed_in_y = False) 

list_of_nodes = [node_1, node_2, node_3, node_4]

E = 200e9
A = 900e-6

bar_1 = Bar(left_node = node_1, right_node = node_2, q = 0, E = E, A = A)
bar_2 = Bar(left_node = node_3, right_node = node_2, q = 0, E = E, A = A)
bar_3 = Bar(left_node = node_1, right_node = node_3, q = 0, E = E, A = A)
bar_4 = Bar(left_node = node_1, right_node = node_4, q = 0, E = E, A = A)
bar_5 = Bar(left_node = node_4, right_node = node_2, q = 0, E = E, A = A)
bar_6 = Bar(left_node = node_3, right_node = node_4, q = 0, E = E, A = A)

list_of_bars = [bar_1, bar_2, bar_3, bar_4, bar_5, bar_6]

truss = Structure(list_of_nodes = list_of_nodes, list_of_bars = list_of_bars)

print(f"\n\n {truss.getSolution()}")
# print(f"\n\n {truss.getBarsStressesAndNormals()}")


# for i, node in enumerate(list_of_nodes):      
#     print(f'node_{i+1}:\n')
#     print(f'u_{i+1} = {node.getDisplacement()[0]}, v_{i+1} = {node.getDisplacement()[1]}\n')
#     print(f'H_{i+1} = {node.getTotalForces()[0]}, V_{i+1} = {node.getTotalForces()[1]}\n')

# for i, bar in enumerate(list_of_bars):
#     print(f'bar_{i+1}:\n')
#     print(f'N_{i+1} = {bar.getBarNormal()}\n')
    # print(f'sigma_{i+1} = {bar.getBarStress()}\n')

# for i, node in enumerate(list_of_nodes):
#     print(f'node_{i+1}:')
#     print(f'u_{i+1} = {node.getDisplacement()[0]}, v_{i+1} = {node.getDisplacement()[1]}')
#     print(f'H_{i+1} = {node.getTotalForces()[0]}, V_{i+1} = {node.getTotalForces()[1]}')    

truss.plotStructure(displacement_scale = 1000)
# print(truss.nodes_initial_positions == truss.nodes_final_positions)
# print(truss.nodes_final_positions)
