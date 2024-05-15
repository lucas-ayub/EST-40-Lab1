from node import *
from bar import *
from structure import *

E1 = 0.8e9
A1 = 1000e-6
E2 = 2.7e9
A2 = 1875e-6
P = 25e3
q = 10e3
L1 = 150e-3
L2 = 262.5e-3

q1 = -2*q
q2 = q

node_1 = Node(x=0, fe=0, fixed=True)
node_2 = Node(x=L1, fe=-2*P, fixed=False)
node_3 = Node(x=L1 + L2, fe=0, fixed=True)

list_of_nodes = [node_1, node_2, node_3]

bar_1 = Bar(left_node=node_1, right_node=node_2, E=E1, A=A1, q=q1)
bar_2 = Bar(left_node=node_2, right_node=node_3, E=E2, A=A2, q=q2)

list_of_bars = [bar_1, bar_2]

axial_bar = Structure(list_of_nodes=list_of_nodes, list_of_bars=list_of_bars)


solution = axial_bar.solve() 
print(f"\n\n {solution}")
values = list(solution.values())
print(values[0]+values[1])
