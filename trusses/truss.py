from node import Node
from bar import Bar

import numpy as np
import sympy as sp


class Truss:
    def __init__(self, list_of_nodes, list_of_bars):
        """
        Initializes a Truss object.

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
        self.symbols = []
        self.r, self.u = self.calculateReactionAndDisplacementVectors()
        self.solution = self.calculateSolution()
        self.setNodesDisplacementsAndForces()
    
    def setNodesDisplacementsAndForces(self):
        """
        Sets the displacements and forces of the nodes based on the system solution.
        """
        solution = self.getSolution()
        values = list(solution.values())  
        print(values)
        for i in range(1, self.num_nodes + 1):
            index = (i - 1) * 4
            H, V, u, v = values[index:index+4]
            self.nodes[i - 1].setDisplacement(u, v)
            self.nodes[i - 1].setTotalForces(H, V)
 
        
    def calculateStiffnessMatrix(self):
        """
        Calculates the global stiffness matrix of the system.

        :return: The global stiffness matrix.
        :rtype: numpy.ndarray
        """
        global_matrix = np.zeros((2 * self.num_nodes, 2 * self.num_nodes))

        for bar in self.bars:
            i = self.nodes.index(bar.left_node)
            j = self.nodes.index(bar.right_node)
            
            local_matrix = bar.getStiffnessMatrix()
            
            global_matrix[2*i:2*(i+1), 2*i:2*(i+1)] += local_matrix[:2, :2]
            global_matrix[2*i:2*(i+1), 2*j:2*(j+1)] += local_matrix[:2, 2:]
            global_matrix[2*j:2*(j+1), 2*i:2*(i+1)] += local_matrix[2:, :2]
            global_matrix[2*j:2*(j+1), 2*j:2*(j+1)] += local_matrix[2:, 2:]

        return global_matrix

            
    def calculateForceVector(self):
        """
        Calculates the force vector of the system.

        :return: The force vector.
        :rtype: numpy.ndarray
        """
        force_vector = np.zeros(2 * len(self.nodes))

        for i, node in enumerate(self.nodes):
            force_vector[2*i] += node.external_forces[0]  
            force_vector[2*i + 1] += node.external_forces[1]  

        return force_vector
            
    def calculateReactionAndDisplacementVectors(self):
        """
        Calculates the reaction force vector and the displacement vector of the system.

        :return: The reaction force vector and the displacement vector.
        :rtype: sympy.Matrix, sympy.Matrix
        """
        reaction_vector = sp.zeros(2*self.num_nodes, 1)
        displacement_vector = sp.zeros(2*self.num_nodes, 1)
        
        for i, node in enumerate(self.nodes):
            H_symbol, V_symbol = 'H_' + str(i+1), 'V_' + str(i+1)
            u_symbol, v_symbol = 'u_' + str(i+1), 'v_' + str(i+1)
            
            self.symbols.extend([H_symbol, V_symbol, u_symbol, v_symbol])
            
            is_fixed_x, is_fixed_y = node.getFixedState()
            
            reaction_vector[2*i] = H_symbol if is_fixed_x else 0
            displacement_vector[2*i] = 0 if is_fixed_x else u_symbol
            reaction_vector[2*i + 1] = V_symbol if is_fixed_y else 0
            displacement_vector[2*i + 1] = 0 if is_fixed_y else v_symbol
                         
        self.symbols = sp.symbols(self.symbols) 
            
        return reaction_vector, displacement_vector

    def sympySolution(self):
        """
        Returns the sympy solution to the system.
        
        :return: Symply dict with the variables and their values
        :rtype: dict
        """
        f = sp.Matrix(self.f) + self.r
        K = sp.Matrix(self.K)
        u = sp.Matrix(self.u)
        variables = sp.solve(K*u - f, self.symbols)
        
        return variables

    def calculateSolution(self):
        """
        Calculates the system solution.
        
        :return: The solution vector.
        :rtype: dict
        """    
        variables = self.sympySolution()
        solution = {key: 0 for key in self.symbols}        
        for key in solution.keys():
            if key in variables.keys():
                solution[key] = variables[key]
        return solution
    
    def getSolution(self):
        """
        Returns the system solution.
        
        :return: The solution vector.
        :rtype: dict
        """
        
        return self.solution