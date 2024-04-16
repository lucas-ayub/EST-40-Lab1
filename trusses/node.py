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