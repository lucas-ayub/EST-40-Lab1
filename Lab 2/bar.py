import numpy as np
from node import *
from scipy.integrate import quad


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
    def __init__(self, left_node, right_node, E, A, I, global_q_x, local_q_x, global_q_y, local_q_y):
        """
        Initializes a Bar object.

        :param left_node: The left node of the bar.
        :type left_node: _Node_
        :param right_node: The right node of the bar.
        :type right_node: _Node_
        :param E: Young's modulus of the material.
        :type E: float
        :param A: Cross-sectional area of the bar.
        :type A: float
        :param I: Inertia Momentum of the bar.
        :type I: float
        :param global_q_x: Distributed load in the x-direction (Global).
        :type global_q_x: float
        :param local_q_x: Distributed load in the x-direction (Local).
        :type local_q_x: float
        :param global_q_y: Distributed load in the y-direction (Global).
        :type global_q_y: float
        :param local_q_y: Distributed load in the y-direction (Local).
        :type local_q_y: float
        """
        self.left_node = left_node
        self.right_node = right_node
        self.E = E
        self.A = A
        self.I = I
        self.L = self.calculateBarLength()
        self.Li = self.calculateBarLength()
        self.angle = self.getBarAngle()
        self.stiffness_matrix = self.calculateStiffnessMatrix()
        self.global_loads = [global_q_x, global_q_y]
        self.local_loads = [local_q_x, local_q_y]
        self.rotation_matrix = self.calculateRotationMatrix()
        self.force_vector = self.calculateForceVector()
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
    
    def calculateForceVector(self):
        """
        Calculates the force vector of the bar.

        :return: The force vector.
        :rtype: numpy.ndarray
        """
        force_vector = np.zeros(6)
        c = np.cos(self.angle)
        s = np.sin(self.angle)
        
        phi_1 = lambda x: 2 * (x**3) / self.L**3 - 3 * (x**2) / self.L**2 + 1
        phi_2 = lambda x: x - 2 * (x**2) / self.L + (x**3) / self.L**2
        phi_3 = lambda x: -2 * (x**3) / self.L**3 + 3 * (x**2) / self.L**2
        phi_4 = lambda x: -(x**3) / self.L**2 + (x**2) / self.L
        
        
        phi_5 = lambda x: 1 - x / self.L
        phi_6 = lambda x: x / self.L

        decomposed_q_x = lambda x: self.global_q_x(x) * c - self.global_q_y(x) * s
        decomposed_q_y = lambda x: self.global_q_x(x) * s + self.global_q_y(x) * c

        total_q_x = lambda x: self.local_q_x(x) + decomposed_q_x(x)
        total_q_y = lambda x: self.local_q_y(x) + decomposed_q_y(x)

        self.force_vector[0] = quad(lambda x: total_q_x(x) * phi_5(x), 0, self.L)[0]
        self.force_vector[1] = quad(lambda x: total_q_y(x) * phi_1(x), 0, self.L)[0]
        self.force_vector[2] = quad(lambda x: total_q_y(x) * phi_2(x), 0, self.L)[0]
        self.force_vector[3] = quad(lambda x: total_q_x(x) * phi_6(x), 0, self.L)[0]
        self.force_vector[4] = quad(lambda x: total_q_y(x) * phi_3(x), 0, self.L)[0]
        self.force_vector[5] = quad(lambda x: total_q_y(x) * phi_4(x), 0, self.L)[0]

        global_force_vector = self.rotation_matrix @ force_vector
        
        return global_force_vector
                   
    def calculateStiffnessMatrix(self):
        """
        Calculates the stiffness matrix of the bar.

        :return: The stiffness matrix.
        :rtype: numpy.ndarray
        """
        c = np.cos(self.angle)
        s = np.sin(self.angle)
        mu = self.A * self.L**2 / 2 * self.I
        k = 2 * self.E * self.I / self.L**3
        L = self.L
        K = k * np.array([
            [mu * c**2 + 6 * s**2, (mu - 6) * c * s, -3 * L * s, -mu * c**2 - 6 * s**2, -(mu - 6) * c * s, -3 * L * s],
            [(mu - 6) * c * s, mu * s**2 + 6 * c**2, 3 * L * c, -(mu - 6) * c * s, -mu * s**2 - 6 * c**2, 3 * L * c],
            [-3 * L * s, 3 * L * c, 2 * L**2, 3 * L * s, -3 * L * c, L**2],
            [-mu * c**2 - 6 * s**2, -(mu - 6) * c * s, 3 * L * s, mu * c**2 + 6 * s**2, (mu - 6) * c * s, 3 * L * s],
            [-(mu - 6) * c * s, -mu * s**2 - 6 * c**2, -3 * L * c, (mu - 6) * c * s, mu * s**2 + 6 * c**2, -3 * L * c],
            [-3 * L * s, 3 * L * c, L**2, 3 * L * s, -3 * L * c, 2 * L**2]
        ])
        
        return K
    
    def calculateRotationMatrix(self):
        """
        Calculates the rotation matrix as shown in the image.

        :return: The rotation matrix.
        :rtype: numpy.ndarray
        """
        c = np.cos(self.angle)
        s = np.sin(self.angle)

        R = np.array([
            [c, -s, 0, 0, 0, 0],
            [s,  c, 0, 0, 0, 0],
            [0,  0, 1, 0, 0, 0],
            [0,  0, 0, c, -s, 0],
            [0,  0, 0, s,  c, 0],
            [0,  0, 0, 0,  0, 1]
        ])

        return R
    
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