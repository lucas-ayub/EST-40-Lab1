from node import Node
from bar import Bar
import numpy as np

def getLinearLoad(start_value, final_value, left_node, right_node):
    """
    Initializes a function to get a linear load function.

    :param start_value: The initial value of the load.
    :type start_value: float
    :param final_value: The final value of the load.
    :type final_value: float
    :param left_node: The left node of the bar.
    :type left_node: Node
    :param right_node: The right node of the bar.
    :type right_node: Node

    :return: A function that computes the load at a given point x, where 0 <= x <= L.
    :rtype: function
    """
    L = np.linalg.norm(left_node.position - right_node.position)
    return lambda x: start_value + (final_value - start_value) * (x / L)

def discretizateBar(left_node, right_node, q, E, A, n_elements):
    """
    Initializes a function to discretize a bar between two given nodes.

    :param left_node: The left node of the bar.
    :type left_node: Node
    :param right_node: The right node of the bar.
    :type right_node: Node
    :param q: Distributed load along the bar.
    :type q: float
    :param E: Young's modulus of the material.
    :type E: float
    :param A: Cross-sectional area of the bar.
    :type A: float
    :param n_elements: Number of elements to discretize the bar into (number of nodes + 1 between left and right node).
    :type n_elements: int

    :return: List of discretized bars and list of intermediate nodes between left_node and right_node.
    :rtype: tuple
    """
    left_node_position, right_node_position = left_node.getPosition(), right_node.getPosition()
    left_node_x, left_node_y = left_node_position[0], left_node_position[1]
    right_node_x, right_node_y = right_node_position[0], right_node_position[1]
    dx = (right_node_x - left_node_x) / n_elements
    dy = (right_node_y - left_node_y) / n_elements
    bars = []
    nodes = []
    for i in range(1, n_elements):
        node = Node(x = left_node_x + i*dx, y = left_node_y + i*dy, fx = 0, fy = 0, fixed_in_x = False, fixed_in_y = False)
        nodes.append(node)
    for i in range(n_elements):
        node_i = left_node if i == 0 else nodes[i-1]
        node_j = right_node if i == n_elements-1 else nodes[i]
        bar = Bar(left_node = node_i, right_node = node_j, q = q, E = E, A = A)
        bars.append(bar)
    return bars, nodes