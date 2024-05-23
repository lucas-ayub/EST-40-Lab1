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

b1 = np.linspace(0, L1, 5)
b2 = np.linspace(L1, L1 + L2, 5)[1:]

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

load = getLinearLoad(-90e3, -30e3, L1 + L2)

node_e_1 = Node(index=1, x=b1[0], y=0, support_type='fixed', prescribed_displacement_x=0, prescribed_displacement_y=0, prescribed_rotation=0)
node_2 = Node(index=2, x=b1[1], y=0, support_type='free')
node_3 = Node(index=3, x=b1[2], y=0, support_type='free')
node_4 = Node(index=4, x=b1[3], y=0, support_type='free')
node_e_5 = Node(index=5, x=b1[4], y=0, support_type='horizontal_roller', prescribed_displacement_y=0)

node_6 = Node(index=6, x=b2[0], y=0, support_type='free')
node_7 = Node(index=7, x=b2[1], y=0, support_type='free')
node_8 = Node(index=8, x=b2[2], y=0, support_type='free')
node_e_9 = Node(index=9, x=b2[3], y=0, support_type='free')

nodes = [node_e_1, node_2, node_3, node_4, node_e_5, node_6, node_7, node_8, node_e_9]


bars = []
# Loop para a primeira parte da estrutura (b1)
for i in range(4):
    q0 = load(b1[i])
    qf = load(b1[i + 1])
    bars.append(Bar(index=i + 1, left_node=nodes[i], right_node=nodes[i + 1], E=E, A=A, I=I, local_q_y=getLinearLoad(q0, qf, b1[i + 1] - b1[i])))

# Loop para a segunda parte da estrutura (b2)
for i in range(3):  # Ajustado para evitar o acesso fora dos limites
    q0 = load(b2[i])
    qf = load(b2[i + 1])
    bars.append(Bar(index=i + 5, left_node=nodes[i + 4], right_node=nodes[i + 5], E=E, A=A, I=I, local_q_y=getLinearLoad(q0, qf, b2[i + 1] - b2[i])))

# Adicionando a última barra manualmente para garantir que todas as barras sejam incluídas
q0 = load(b2[3])
qf = load(b2[3])  # O último ponto não tem um ponto final seguinte, então a carga final será igual à carga inicial
bars.append(Bar(index=8, left_node=nodes[7], right_node=nodes[8], E=E, A=A, I=I, local_q_y=getLinearLoad(q0, qf, b2[3] - b2[2])))

structure = Structure(list_of_nodes=nodes, list_of_bars=bars)

solution = structure.solution

structure.printNodalParameters()
structure.printBarParameters()
