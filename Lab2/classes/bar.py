import numpy as np
import sympy as sp
from scipy.integrate import quad
from node import Node

class Bar:
    def __init__(self,index, left_node, right_node, E, A, I, global_q_x=0, local_q_x=0, global_q_y=0, local_q_y=0):
        """
        Initializes a Bar object.

        :param index: The index of the bar.
        :type index: int
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
        :type global_q_x: function
        :param local_q_x: Distributed load in the x-direction (Local).
        :type local_q_x: function
        :param global_q_y: Distributed load in the y-direction (Global).
        :type global_q_y: function
        :param local_q_y: Distributed load in the y-direction (Local).
        :type local_q_y: function
        """
        self.index = index 
        self.left_node = left_node
        self.right_node = right_node
        self.E = E
        self.A = A
        self.I = I
        self.L = self.calculateBarLength()
        self.Li = self.calculateBarLength()
        self.angle = self.getBarAngle()
        self.rotation_matrix_6x6, self.rotation_matrix_3x3 = self.calculateRotationMatrix()
        self.stiffness_matrix = self.calculateStiffnessMatrix()
        self.global_loads = [
            (lambda x, val=global_q_x: val if isinstance(global_q_x, (int, float)) else global_q_x(x)),
            (lambda y, val=global_q_y: val if isinstance(global_q_y, (int, float)) else global_q_y(y))
        ]
        self.local_loads = [
            (lambda x, val=local_q_x: val if isinstance(local_q_x, (int, float)) else local_q_x(x)),
            (lambda y, val=local_q_y: val if isinstance(local_q_y, (int, float)) else local_q_y(y))
        ]
        self.force_vector = self.calculateForceVector()
        self.N = 0
        self.sigma = 0
        self.V = 0
        self.M = 0
        
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

    def calculateRotationMatrix(self):
        """
        Calculates the 6x6 and 3x3 rotation matrices for the given angle.

        The rotation matrices are used to rotate a vector or a set of points in 3D space.
        The 6x6 rotation matrix is used to rotate a 6-dimensional vector, while the 3x3
        rotation matrix is used to rotate a 3-dimensional vector.

        :return: A tuple containing the 6x6 and 3x3 rotation matrices.
        :rtype: tuple(numpy.ndarray, numpy.ndarray)
        """
        c = np.cos(self.angle)
        s = np.sin(self.angle)

        R_6 = np.array([
            [c, -s, 0, 0, 0, 0],
            [s,  c, 0, 0, 0, 0],
            [0,  0, 1, 0, 0, 0],
            [0,  0, 0, c, -s, 0],
            [0,  0, 0, s,  c, 0],
            [0,  0, 0, 0,  0, 1]
        ])
        
        R_3 = np.array([
            [c, -s, 0],
            [s,  c, 0],
            [0,  0, 1]
        ])


        return R_6, R_3
    
    def calculateForceVector(self):
        """
        Calculates the force vector of the bar.

        :return: The force vector.
        :rtype: numpy.ndarray
        """
        force_vector = np.zeros(6)
        c = np.cos(self.angle)
        s = np.sin(self.angle)
        
        phi_1 = lambda x: x / self.L
        phi_2 = lambda x: 2 * (x**3) / self.L**3 - 3 * (x**2) / self.L**2 + 1
        phi_3 = lambda x: (x**3) / (self.L**2) - 2 * (x**2) / (self.L) + x
        
        phi_4 = lambda x: 1 - x / self.L
        phi_5 = lambda x: -2 * (x**3) / self.L**3 + 3 * (x**2) / self.L**2
        phi_6 = lambda x: (x**3) / self.L**2 - (x**2) / self.L
        
        decomposed_q_x = lambda x: self.global_loads[0](x) * c - self.global_loads[1](x) * s
        decomposed_q_y = lambda x: self.global_loads[0](x) * s + self.global_loads[1](x) * c
        total_q_x = lambda x: self.local_loads[0](x) + decomposed_q_x(x)
        total_q_y = lambda x: self.local_loads[1](x) + decomposed_q_y(x)

        force_vector[0] = quad(lambda x: total_q_x(x) * phi_1(x), 0, self.L)[0]
        force_vector[1] = quad(lambda x: total_q_y(x) * phi_2(x), 0, self.L)[0]
        force_vector[2] = quad(lambda x: total_q_y(x) * phi_3(x), 0, self.L)[0]
        force_vector[3] = quad(lambda x: total_q_x(x) * phi_4(x), 0, self.L)[0]
        force_vector[4] = quad(lambda x: total_q_y(x) * phi_5(x), 0, self.L)[0]
        force_vector[5] = quad(lambda x: total_q_y(x) * phi_6(x), 0, self.L)[0]

        global_force_vector = self.rotation_matrix_6x6 @ force_vector

        return global_force_vector

    
    def getForceVector(self):
        """
        Get the force vector of the bar.
        
        :return: The force vector.
        :rtype: numpy.ndarray
        """
        return self.force_vector
                   
    def calculateStiffnessMatrix(self):
        """
        Calculates the stiffness matrix of the bar.

        :return: The stiffness matrix.
        :rtype: numpy.ndarray
        """
        c = np.cos(self.angle)
        s = np.sin(self.angle)
        mu = self.A * self.L**2 / (2 * self.I)
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
        
    def getStiffnessMatrix(self):
        """
        Gets the stiffness matrix of the bar.

        :return: The stiffness matrix.
        :rtype: numpy.ndarray
        """
        return self.stiffness_matrix

    def calculateBarNormalAndStress(self):
        """
        Calculates the stress and normal force of the bar.
        """
        
        left_displacements = self.left_node.getDisplacement()
        right_displacements = self.right_node.getDisplacement()
    
        R_T = self.rotation_matrix_3x3.T
        
        local_left_displacements = R_T @ left_displacements
        local_right_displacements = R_T @ right_displacements
        
        u_l1 = local_left_displacements[0]
        u_l2 = local_right_displacements[0]
        
        sigma = self.E * (u_l1 - u_l2) / self.Li 
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

    def calculateBarShearForceAndBendingMoment(self):
        """
        Calculates the shear force and bending moment of the bar.
        """
        
        left_displacements = self.left_node.getDisplacement()
        right_displacements = self.right_node.getDisplacement()
    
        R_T = self.rotation_matrix_3x3.T
        
        local_left_displacements = R_T @ left_displacements
        local_right_displacements = R_T @ right_displacements
        
        v_l1 = local_left_displacements[1]  
        theta_l2 = local_left_displacements[2]
        v_l3 = local_right_displacements[1]
        theta_l4 = local_right_displacements[2]

        E = self.E
        I = self.I
        L = self.Li
        
        self.V = E * I * (v_l1 * (12 / L**3) + theta_l2 * (6 / L**2) + v_l3 * (-12 / L**3) + theta_l4 * (6 / L**2) )
        
        self.M = lambda x: E * I * (v_l1 * (12 * x / L**3 - 6 / L**2) + theta_l2 * (6 * x / L**2 - 4 / L) + v_l3 * (-12 * x / L**3 + 6 / L**2) + theta_l4 * (6 * x / L**2 - 2 / L))
    
    def getBarShear(self):
        """
        Gets the bar shear force.
        
        :return: The shear force of the bar.
        :rtype: float
        """
        return self.V
    
    def getBarBendingMomentum(self):
        """
        Gets the bar bending moment.
        
        :return: The bending moment of the bar.
        :rtype: function
        """
        return self.M
        
    def getBendingMomentumExpression(self):
        """
        Gets the bending moment expression of the bar.
        
        :return: The symbolic bending moment expression of the bar.
        :rtype: function
        """
        x = sp.symbols('x')

        M_symbolic = self.M(x)
        
        return M_symbolic
    
    