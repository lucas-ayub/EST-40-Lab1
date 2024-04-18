import numpy as np

class Node:
    def __init__(self, x, y, fy, fx, fixed_in_x, fixed_in_y):
        """
        Initializes a Node object.

        :param x: The x-coordinate of the node.
        :type x: float
        :param y: The y-coordinate of the node.
        :type y: float
        :param fy: The vertical external force applied to the node.
        :type fy: float
        :param fx: The horizontal external force applied to the node.
        :type fx: float
        :param fixed_in_x: Indicates whether the node is fixed in the x-direction.
        :type fixed_in_x: bool
        :param fixed_in_y: Indicates whether the node is fixed in the y-direction.
        :type fixed_in_y: bool
        """
        self.position = np.array([x, y])
        self.displacement = np.array([0, 0])
        self.external_forces = np.array([fx, fy])
        self.constraints = {'isFixedInX': fixed_in_x, 'isFixedInY': fixed_in_y}
   
    def getPosition(self):
        """
        Gets the position of the node.
        
        :return: The position of the node.
        :rtype: np.array
        """    
        return self.position
        
    def setDisplacement(self, delta_x, delta_y):
        """
        Sets the displacement of the node.

        :param delta_x: The horizontal displacement of the node.
        :type delta_x: float
        :param delta_y: The vertical displacement of the node.
        :type delta_y: float
        """
        self.displacement = np.array([delta_x, delta_y])
        
    def getDisplacement(self):
        """
        Gets the displacement of the node.

        :return: The displacement of the node.
        :rtype: numpy.ndarray
        """
        return self.displacement
        
    def setTotalForces(self, r_x, r_y):
        """
        Sets the displacement of the node.

        :param r_x: The total horizontal force of the node.
        :type r_x: float
        :param r_y: The total vertical force of the node.
        :type r_y: float
        """
        self.external_forces = np.array([r_x, r_y])
        
    def getTotalForces(self):
        """
        Gets the total forces of the node.

        :return: The total forces of the node.
        :rtype: numpy.ndarray
        """
        return self.external_forces
    
    def getFixedState(self):
        """
        Gets the fixed state of the node.

        :return: The fixed state of the node.
        :rtype: bool
        """
        return (self.constraints['isFixedInX'], self.constraints['isFixedInY'])