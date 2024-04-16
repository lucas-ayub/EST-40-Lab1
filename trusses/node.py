import numpy as np

class Node:
    def __init__(self, x, y, fy, fx, fixed_in_x, fixed_in_y):
        self.position = np.array([x, y])
        self.displacement = np.array([0, 0])
        self.neighbors = []
        self.constraints = {'isFixedInX': fixed_in_x, 'isFixedInY': fixed_in_y}
        self.external_forces = np.array([fx, fy])
   
        
    def newNeighbor(self, neighbor):
        """Add a new neighbor to the node."""
        self.neighbors.append(neighbor)
        
    def setDisplacement(self, delta_x, delta_y):
        """
        Sets the displacement of the node.

        :param delta_x: The horizontal displacement of the node.
        :type delta_x: float
        :param delta_y: The vertical displacement of the node.
        :type delta_y: float
        """
        self.displacement = np.array([delta_x, delta_y])
    
    def getFixedState(self):
        """
        Gets the fixed state of the node.

        :return: The fixed state of the node.
        :rtype: bool
        """
        return (self.constraints['isFixedInX'], self.constraints['isFixedInY'])