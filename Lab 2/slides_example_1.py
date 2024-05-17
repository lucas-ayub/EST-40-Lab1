from node import *
from bar import *
from structure import *

node_1 = Node(x=0, y=0, support_type='fixed')

node_2 = Node(x=3, y=2.25, support_type='free', global_f_x=50e3, global_f_y=25e3)

node_3 = Node(x=4, y=3, support_type='double_roller', prescribed_displacement_x=2e-3, prescribed_displacement_y=0)

nodes = [node_1, node_2, node_3]

bar_1 = Bar(left_node=node_1, right_node=node_2, E=25e9, A=0.036, I=2.7e-4, local_q_x=lambda x: 40e3)
bar_2 = Bar(left_node=node_2, right_node=node_3, E=25e9, A=0.036, I=2.7e-4, local_q_y=lambda y: 20e3)

bars = [bar_1, bar_2]

structure = Structure(list_of_nodes=nodes, list_of_bars=bars)