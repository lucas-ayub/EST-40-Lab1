import sys
import os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../classes')))

from bar import Bar
from node import Node
from structure import Structure

xi = 500
xf = 700
yf = -150

M = 8000

node_1 = Node(x=0, y=0, support_type='fixed', prescribed_displacement_x=0, prescribed_displacement_y=0, prescribed_rotation=0)
node_2 = Node(x=xi, y=0, support_type='free', external_momentum=M)
node_3 = Node(x=xf, y=yf, support_type='fixed', prescribed_displacement_x=0, prescribed_displacement_y=0, prescribed_rotation=0)

nodes = [node_1, node_2, node_3]

E = 20000
A = 60 
I = 500

q1_y = -0.4

bar_1 = Bar(left_node=node_1, right_node=node_2, E=E, A=A, I=I, local_q_y=q1_y)
bar_2 = Bar(left_node=node_2, right_node=node_3, E=E, A=A, I=I)

bars = [bar_1, bar_2]

structure = Structure(list_of_nodes=nodes, list_of_bars=bars)

print(structure.solution)
