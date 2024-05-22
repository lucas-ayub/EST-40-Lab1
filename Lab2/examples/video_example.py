import sys
import os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../classes')))

from bar import Bar
from node import Node
from structure import Structure

node_1 = Node(index=1, x=0, y=0, support_type='fixed', prescribed_displacement_x=0, prescribed_displacement_y=0, prescribed_rotation=0)

node_2 = Node(index=2, x=3, y=2.25, support_type='free', global_f_x=50e3, global_f_y=25e3)

node_3 = Node(index=3, x=4, y=3, support_type='double_roller', prescribed_displacement_x=2e-3, prescribed_displacement_y=0)

nodes = [node_1, node_2, node_3]

bar_1 = Bar(left_node=node_1, right_node=node_2, E=25e9, A=0.036, I=2.7e-4, local_q_x=40e3)
bar_2 = Bar(left_node=node_2, right_node=node_3, E=25e9, A=0.036, I=2.7e-4, local_q_y=20e3)

bars = [bar_1, bar_2]

structure = Structure(list_of_nodes=nodes, list_of_bars=bars)

print(structure.solution)
