import numpy as np
from node import Node

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
        """
        self.left_node = left_node
        self.right_node = right_node
        self.E = E
        self.A = A
        self.L = np.linalg.norm(self.left_node.position - self.right_node.position)
        self.load = q
        self.k = self.E * self.A / self.L
        self.stiffness_matrix = self.k * np.array([[1, -1], [-1, 1]])
        
    def getStiffnessMatrix(self):
        """
        Calculates the stiffness matrix of the bar.

        :return: The stiffness matrix.
        :rtype: numpy.ndarray
        """
        return self.stiffness_matrix
    
    def updateNodeForces(self):
        """
        Computes the load forces in the bar nodes.
        """
        self.left_node.addNewForce(self.load*self.L)
        self.right_node.addNewForce(self.load*self.L)