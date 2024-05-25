import sys
import os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../classes')))

from bar import Bar
from node import Node
from structure import Structure

L1 = 5.5
L2 = 2.5

E = 200e9
A = 4812e-6
I = 1.24204e-4

def getLinearLoad(start_value, final_value, L):
    """
    Initializes a function to get a linear load function.

    :param start_value: The initial value of the load.
    :type start_value: float
    :param final_value: The final value of the load.
    :type final_value: float
    :param L: Bar length.
    
    :return: A function that computes the load at a given point x, where 0 <= x <= L.
    :rtype: function
    """
    return lambda x: start_value + (final_value - start_value) * (x / L)

load = getLinearLoad(90e3, 30e3, L1 + L2)

node_1 = Node(index=1, x=0, y=0, support_type='fixed', prescribed_displacement_x=0, prescribed_displacement_y=0, prescribed_rotation=0)

node_2 = Node(index=2, x=L1/3, y=0, support_type='free')
node_3 = Node(index=3, x=2*L1/3, y=0, support_type='free')

node_4 = Node(index=4, x=L1, y=0, support_type='horizontal_roller', prescribed_displacement_y=0)

node_5 = Node(index=5, x=L1 + 0.5*L2, y=0, support_type='free')
node_6 = Node(index=6, x=L1 + L2, y=0, support_type='free')

nodes = [node_1, node_2, node_3, node_4, node_5, node_6]

q0 = load(0)
qf = load(L1/3)
load_1 = getLinearLoad(q0, qf, L1/3)

bar_1 = Bar(index = 1, left_node = node_1, right_node = node_2, E = E, A = A, I = I, local_q_y = lambda x: -load_1(x))

q0 = load(L1/3)
qf = load(2*L1/3)
load_2 = getLinearLoad(q0, qf, L1/3)

bar_2 = Bar(index = 2, left_node = node_2, right_node = node_3, E = E, A = A, I = I, local_q_y = lambda x: -load_2(x))

q0 = load(2*L1/3)
qf = load(L1)
load_3 = getLinearLoad(q0, qf, L1/3)

bar_3 = Bar(index = 3, left_node = node_3, right_node = node_4, E = E, A = A, I = I, local_q_y = lambda x: -load_3(x))

q0 = load(L1)
qf = load(L1 + 0.5*L2)
load_4 = getLinearLoad(q0, qf, 0.5*L2)

bar_4 = Bar(index = 4, left_node = node_4, right_node = node_5, E = E, A = A, I = I, local_q_y = lambda x: -load_4(x))

q0 = load(L1 + 0.5*L2)
qf = load(L1 + L2)
load_5 = getLinearLoad(q0, qf, 0.5*L2)

bar_5 = Bar(index = 5, left_node = node_5, right_node = node_6, E = E, A = A, I = I, local_q_y = lambda x: -load_5(x))

bars = [bar_1, bar_2, bar_3, bar_4, bar_5] 

structure = Structure(list_of_nodes = nodes, list_of_bars = bars)

solution = structure.solution

print('========================================================================')
print('Exemplo 1:\n')
structure.printNodalParameters()
structure.printBarParameters()
