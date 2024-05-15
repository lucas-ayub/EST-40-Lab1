import numpy as np

class Node:
    def __init__(self, x, y, global_f_x, local_f_x, global_f_y, local_f_y, support_type, prescribed_displacement_x, prescribed_displacement_y, prescribed_rotation):
        """
        Initializes a Node object.
        :param x: The x-coordinate of the node.
        :type x: float
        :param y: The y-coordinate of the node.
        :type y: float
        :param global_f_x: The horizontal external force applied to the node.
        :type global_f_x: float
        :param local_f_x: The horizontal external force applied to the node in local coordinates.
        :type local_f_x: float
        :param global_f_y: The vertical external force applied to the node.
        :type global_f_y: float
        :param local_f_y: The vertical external force applied to the node in local coordinates.
        :type local_f_y: float
        :param support_type: Indicates the support type.
        :type support_type: str
        :param prescribed_displacement_x: The prescribed displacement in x of the node.
        :type prescribed_displacement_x: float
        :param prescribed_displacement_y: The prescribed displacement in y of the node.
        :type prescribed_displacement_y: float
        :param prescribed_rotation: The prescribed rotation of the node.
        :type prescribed_rotation: float
        """
        self.position = np.array([x, y], dtype=np.float64)
        self.displacement = np.array([0, 0], dtype=np.float64)
        self.support = support_type
        self.global_forces = np.array([global_f_x, global_f_y], dtype=np.float64)
        self.local_forces = np.array([local_f_x, global_f_y], dtype=np.float64)
        # TODO: Prescribed as a list instead of a array to modify it later and to differ if the user knows a value of displacement or if it is actually zero
        self.prescribed_displacement = [prescribed_displacement_x,prescribed_displacement_y]
        self.prescribed_rotation = prescribed_rotation
        
        
    def updatePosition(self):
        """
        Updates the position of the node.
        """
        for i in range(len(self.position)):
            self.position[i] = self.position[i] + self.displacement[i]   
   
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
        
    def addNewForce(self, delta_f_x, delta_f_y):
        """
        Updates the forces of the node.

        :param f_x: The horizontal force to be added to the node.
        :type f_x: float
        :param f_y: The vertical force to be added to the node.
        :type f_y: float
        """
        self.external_forces += np.array([delta_f_x, delta_f_y], dtype = np.float64)
        
    def getTotalForces(self):
        """
        Gets the total forces of the node.

        :return: The total forces of the node.
        :rtype: numpy.ndarray
        """
        return self.external_forces
