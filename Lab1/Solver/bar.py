import numpy as np
from node import *

def discretizateBar(left_node, right_node, q, E, A, n_elements):
    """
    Function to discretize a bar between two given nodes.

    Parameters:
    - left_node (Node): The left node of the bar.
    - right_node (Node): The right node of the bar.
    - q (float): Distributed load along the bar.
    - E (float): Young's modulus of the material.
    - A (float): Cross-sectional area of the bar.
    - n_elements (int): Number of elements to discretize the bar into (number of nodes + 1 between left and right node).

    Returns:                  
    - bars (list): List of discretized bars.
    - nodes (list): List of intermediate nodes between left_node and right_node.
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

class Bar:
    def __init__(self, left_node, right_node, q, E, A):
        """
        Initializes a Bar object.

        :param left_node: The left node of the bar.
        :type left_node: _Node_
        :param right_node: The right node of the bar.
        :type right_node: _Node_
        :param q: The distributed load on the bar.
        :type q: float
        :param E: The Young's modulus of the bar material.
        :type E: float
        :param A: The cross-sectional area of the bar.
        :type A: float
        :param N: Normal
        :type N: float
        :param sigma: Stress
        :type sigma: float
        """
        self.left_node = left_node
        self.right_node = right_node
        self.E = E
        self.A = A
        self.L = self.calculateBarLength()
        self.Li = self.calculateBarLength()
        self.angle = self.getBarAngle()
        self.stiffness_matrix = self.calculateStiffnessMatrix()
        self.load = q
        self.N = 0
        self.sigma = 0
        
    def getBarLength(self):
        """
        Get the length of the bar.
        """
        return self.L
    
    def calculateBarLength(self):
        """
        Calculate the length of the bar.
        
        :return: The length of the bar.
        :rtype: float
        """
        return np.linalg.norm(self.left_node.position - self.right_node.position)
        
    def calculateStiffnessMatrix(self):
        """
        Calculates the stiffness matrix of the bar.

        :return: The stiffness matrix.
        :rtype: numpy.ndarray
        """
        c = np.cos(self.angle)
        s = np.sin(self.angle)
        k = self.E * self.A / self.L
        K = k * np.array([[c**2, c * s, -c**2, -c * s],
                          [c * s, s**2, -c * s, -s**2],
                          [-c**2, -c * s, c**2, c * s],
                          [-c * s, -s**2, c * s, s**2]])
        
        return K
    
    def getStiffnessMatrix(self):
        """
        Gets the stiffness matrix of the bar.

        :return: The stiffness matrix.
        :rtype: numpy.ndarray
        """
        return self.stiffness_matrix
    
    def getBarAngle(self):
        """
        Gets the angle of the bar.
        
        :return: The angle of the bar
        :rtype: float
        """
        dx = self.right_node.position[0] - self.left_node.position[0]
        dy = self.right_node.position[1] - self.left_node.position[1]
        
        if dx == 0:  
            if dy > 0:
                return np.pi / 2  
            else:   
                return -np.pi / 2  
        else:
            return np.arctan(dy / dx)
        
    def updateNodeForces(self):
        """
        Computes the load forces in the bar nodes.
        """
        factor = self.load * self.L / 2
        delta_f_x, delta_f_y = factor * np.cos(self.angle), factor * np.sin(self.angle)
        self.left_node.addNewForce(delta_f_x, delta_f_y)
        self.right_node.addNewForce(delta_f_x, delta_f_y)

    def setBarNormalAndStress(self):
        """
        Calculates the stress and normal force of the bar.
        """
        self.L = self.calculateBarLength() 
        sigma = self.E * (self.L - self.Li) / self.Li 
        self.sigma = sigma 
        self.N = sigma * self.A  


    def getBarNormal(self):
        """
        Gets the bar normal.
        
        :return: The normal of the bar.
        :rtype: float
        """
        return self.N
    
    def getBarStress(self):
        """
        Gets the bar stress.
        
        :return: The stress of the bar.
        :rtype: float
        """
        return self.sigma