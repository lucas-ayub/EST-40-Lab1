import numpy as np

class Node:
    def __init__(self, x, fe, fixed):
        """
        Creates a node.

        :param x: node position.
        :type x: float
        :param fe: resultant of concentrated external forces that are not from a wall acting on the node.
        :type fe: float
        :param fixed: indicates if the node is fixed.
        :type fixed: bool
        """
        self.position = x
        self.displacement = None
        self.is_fixed = fixed
        self.f = fe 

    def addNewForce(self, delta_f):
        """
        Adds a new force to the node.

        :param delta_f: The force to be added to the node.
        :type delta_f: float

        This method is used to compute the load forces on the node.
        """
        self.f += delta_f
        
    def setDisplacement(self, delta):
        """
        Sets the displacement of the node.

        :param delta: The displacement of the node.
        :type delta: float
        """
        self.displacement = delta
        
    def getFixedState(self):
        """"""
        return self.is_fixed