import numpy as np
from node import *

class Bar:
    def __init__(self, left_node, right_node, E, A):
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
        self.L = self.calculateLength()
        self.angle = self.getBarAngle()
        self.stiffness_matrix = self.calculateStiffnessMatrix()
        self.N = 0
        
    def calculateLength(self):
        """
        Calculates the length of the bar.
        
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

    def getBarStress(self):
        """
        Gets the stress of the bar.

        :return: The stress of the bar.
        :rtype: float
        """
        Li = self.L
        Lf = self.calculateLength()
        N = self.E * self.A * (Lf - Li) / Li
        self.L = self.calculateLength()
        
        return N
    
    def setBarStress(self):
        """
        Sets the stress of the bar.
        """
        self.N = self.getBarStress()