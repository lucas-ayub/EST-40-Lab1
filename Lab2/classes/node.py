import numpy as np

class Node:
    def __init__(self, index, x, y, support_type, prescribed_displacement_x=None, prescribed_displacement_y=None, prescribed_rotation=None, 
                 global_f_x=0, global_f_y=0, external_momentum=0):
        """
        Initializes a Node object.
        :param index: The index of the node.
        :type index: int
        :param x: The x-coordinate of the node.
        :type x: float
        :param y: The y-coordinate of the node.
        :type y: float
        :param support_type: Indicates the support type:
        horizontal_roller, vertical_roller, double_roller, pinned, fixed, or free.
        :type support_type: str
        :param prescribed_displacement_x: The prescribed displacement in x of the node.
        :type prescribed_displacement_x: float or None
        :param prescribed_displacement_y: The prescribed displacement in y of the node.
        :type prescribed_displacement_y: float or None
        :param prescribed_rotation: The prescribed rotation of the node.
        :type prescribed_rotation: float or None
        :param global_f_x: The horizontal external force applied to the node.
        :type global_f_x: float
        :param global_f_y: The vertical external force applied to the node.
        :type global_f_y: float
        :param external_momentum: The external momentum applied to the node
        :type external_momentum: float
        """
        self.index = index - 1
        self.initial_position = np.array([x, y], dtype=np.float64)
        self.position = np.array([x, y], dtype=np.float64)
        self.displacement = np.array([0, 0, 0], dtype=np.float64)
        self.support = support_type
        self.global_forces = np.array([global_f_x, global_f_y, external_momentum], dtype=np.float64)
        self.prescribed_displacements = [prescribed_displacement_x, prescribed_displacement_y, prescribed_rotation]

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

    def setDisplacement(self, delta_x, delta_y, delta_theta):
        """
        Sets the displacement of the node.

        :param delta_x: The horizontal displacement of the node.
        :type delta_x: float
        :param delta_y: The vertical displacement of the node.
        :type delta_y: float
        :param delta_theta: The rotational displacement of the node.
        :type delta_theta: float
        """
        if self.prescribed_displacements[0] is not None:
            d_x = self.prescribed_displacements[0]
        else:
            d_x = delta_x
        
        if self.prescribed_displacements[1] is not None:
            d_y = self.prescribed_displacements[1]
        else:
            d_y = delta_y
        
        if self.prescribed_displacements[2] is not None:
            d_theta = self.prescribed_displacements[2]
        else:
            d_theta = delta_theta
        
        self.displacement = np.array([d_x, d_y, d_theta], dtype=np.float64)

    def getDisplacement(self):
        """
        Gets the displacement of the node.

        :return: The displacement of the node.
        :rtype: numpy.ndarray
        """
        return self.displacement

    def setFinalForces(self, f_x, f_y, momentum):
        """
        Sets the final forces acting on the node.
        
        :param f_x: The horizontal force.
        :type f_x: float
        :param f_y: The vertical force.
        :type f_y: float
        :param momentum: The moment applied to the node.
        :type momentum: float
        """
        self.global_forces = np.array([f_x, f_y, momentum])

    def getGlobalForces(self):
        """
        Gets the forces acting on the node.
        
        :return: The forces acting on the node.
        :rtype: numpy.ndarray
        """
        return self.global_forces

    def getSupportType(self):
        """
        Gets the support type of the node.
        
        :return: The support type of the node.
        :rtype: str
        """
        return self.support

    def getPrescribedDisplacements(self):
        """
        Gets the prescribed displacements of the node.
        
        :return: The prescribed displacements of the node.
        :rtype: list
        """
        return self.prescribed_displacements
