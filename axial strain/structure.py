from node import Node
from bar import Bar

import numpy as np
import sympy as sp


class Structure:
    def __init__(self, list_of_nodes, list_of_bars):
        """
        Initializes a System object.

        :param list_of_nodes: List of nodes in the system.
        :type list_of_nodes: list of _Node_
        :param list_of_bars: List of bars in the system.
        :type list_of_bars: list of _Bar_
        """
        self.nodes = list_of_nodes
        self.bars = list_of_bars
        self.num_nodes = len(self.nodes)
        self.K = self.calculateStiffnessMatrix()
        self.f = self.calculateForceVector()
        self.r = self.calculateReactionForceVector()
        self.symbols = ['Re', 'Rd']
        
    def calculateStiffnessMatrix(self):
        """
        Calculates the global stiffness matrix of the system.

        :return: The global stiffness matrix.
        :rtype: numpy.ndarray
        """
        global_matrix = np.zeros((self.num_nodes, self.num_nodes))
        
        for i, bar in enumerate(self.bars):
            local_matrix = bar.getStiffnessMatrix()
            global_matrix[i:i+2, i:i+2] += local_matrix
            
        return global_matrix

    def calculateForceVector(self):
        """
        Calculates the force vector of the system.

        :return: The force vector.
        :rtype: numpy.ndarray
        """
        force_vector = np.zeros(self.num_nodes)
        
        for bar in self.bars:
            bar.updateNodeForces()
        for i, bar in enumerate(self.bars):
            if i == 0:
                force_vector[i] += bar.left_node.f
            force_vector[i+1] += bar.right_node.f
            
        return force_vector
    
    def calculateReactionForceVector(self):
        """
        Calculates the reaction force vector of the system.

        :return: The reaction force vector.
        :rtype: sympy.Matrix
        """
        Re, Rd = sp.symbols('Re Rd')
        reaction_vector = sp.zeros(self.num_nodes, 1)
        
        reaction_vector[0] = Re if self.nodes[0].getFixedState() else 0
        reaction_vector[-1] = Rd if self.nodes[-1].getFixedState() else 0
        
        return reaction_vector
    
    def solve(self):
        """
        Solves the system.
        
        :return: The solution vector.
        :rtype: dict
        """
        f = sp.Matrix(self.f) + self.r
        K = sp.Matrix(self.K)
        for i in range(self.num_nodes):
            symbol = 'u_' + str(i+1)
            self.symbols.append(symbol)
        self.symbols = sp.symbols(self.symbols)
        u = sp.Matrix(self.symbols[2:])
        print(u)
        for i in range(len(u)):
            if self.nodes[i].getFixedState():
                u[i] = 0
        solution = sp.solve(K*u - f, self.symbols)
        
        return solution
