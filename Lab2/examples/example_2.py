import sys
import os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../classes')))

from bar import Bar
from node import Node
from structure import Structure

E = 70e9
A1 = (160e-3) * (160e-3)
A2 = (180e-3) * (180e-3)
I1 = (160e-3) * (160e-3) ** 3 / 12
I2 = (180e-3) * (180e-3) ** 3 / 12

node_1 = Node(index=1, x=0, y=0, support_type='fixed', prescribed_displacement_x=0, prescribed_displacement_y=0, prescribed_rotation=0)

node_2 = Node(index=2, x=4/3, y=1, support_type='free')

node_3 = Node(index=3, x=2*4/3, y=2, support_type='free')

node_4 = Node(index=4, x=4, y=3, support_type='free', global_f_x=15e3, global_f_y=-30e3, external_momentum=25e3)

node_5 = Node(index=5, x=7, y=3, support_type='double_roller', prescribed_displacement_x=0, prescribed_displacement_y=5e-3)

node_6 = Node(index=6, x=7, y=6, support_type='free', global_f_x=-12e3)

nodes = [node_1, node_2, node_3, node_4, node_5, node_6]

bar_1 = Bar(index=1, left_node=node_1, right_node=node_2, E=E, A=A1, I=I1, local_q_x=-5e3, local_q_y=-15e3)

bar_2 = Bar(index=2, left_node=node_2, right_node=node_3, E=E, A=A1, I=I1, local_q_x=-5e3, local_q_y=-15e3)

bar_3 = Bar(index=3, left_node=node_3, right_node=node_4, E=E, A=A1, I=I1, local_q_x=-5e3, local_q_y=-15e3)

bar_4 = Bar(index=4, left_node=node_4, right_node=node_5, E=E, A=A2, I=I2)

bar_5 = Bar(index=5, left_node=node_5, right_node=node_6, E=E, A=A2, I=I2, local_q_y=-10e3)

bars = [bar_1, bar_2, bar_3, bar_4, bar_5]

structure = Structure(list_of_nodes=nodes, list_of_bars=bars)

solution = structure.solution


print('========================================================================')
print('Exemplo 2:\n')
structure.printNodalParameters()
structure.printBarParameters()
