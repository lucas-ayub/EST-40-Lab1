from node import *
from bar import *

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


class Structure:
    def __init__(self, list_of_nodes, list_of_bars):
        """
        Initializes a structure object.

        :param list_of_nodes: List of nodes in the system.
        :type list_of_nodes: list of _Node_
        :param list_of_bars: List of bars in the system.
        :type list_of_bars: list of _Bar_
        """
        self.nodes = list_of_nodes
        self.bars = list_of_bars
        self.num_nodes = len(self.nodes)
        self.dof_per_node = 3
        self.K = self.calculateStiffnessMatrix()
        self.f = self.calculateForceVector()
        self.symbols = []
        self.r, self.u = self.calculateReactionAndDisplacementVectors()
        self.solution = self.calculateSolution()
        self.nodes_initial_positions = [node.getPosition() for node in self.nodes]
        self.setNodesDisplacementsAndForces()
        self.updateNodesPositions()
        self.nodes_final_positions = [node.getPosition() for node in self.nodes]
        self.setBarsStressesAndNormals()


    def calculateStiffnessMatrix(self):
        """
        Calculates the global stiffness matrix of the system.

        :return: The global stiffness matrix.
        :rtype: numpy.ndarray
        """
        global_size = self.num_nodes * self.dof_per_node
        global_matrix = np.zeros((global_size, global_size))

        for i, bar_matrix in enumerate(self.bars):
            start_index = i * self.dof_per_node
            end_index = start_index + 2 * self.dof_per_node
            global_matrix[start_index:end_index, start_index:end_index] += bar_matrix

        return global_matrix

    def getStiffnessMatrix(self):
        """
        Returns the global stiffness matrix of the system.

        :return: The global stiffness matrix.
        :rtype: numpy.ndarray
        """
        return self.K
 
    def calculateForceVector(self):
        """
        Calculates the force vector of the system.

        :return: The force vector.
        :rtype: numpy.ndarray
        """
        force_vector = np.zeros(self.dof_per_node * len(self.nodes))

        for i, bar in enumerate(self.bars):
            bar_global_force_vector = bar.getForceVector()
            
            left_node_index = i * self.dof_per_node
            right_node_index = (i + 1) * self.dof_per_node
            
            force_vector[left_node_index:left_node_index+self.dof_per_node] += bar_global_force_vector[:self.dof_per_node]
            force_vector[right_node_index:right_node_index+self.dof_per_node] += bar_global_force_vector[self.dof_per_node:]

        return force_vector
            
    def calculateReactionAndDisplacementVectors(self):
        """
        Calculates the reaction force vector and the displacement vector of the system.

        :return: The reaction force vector and the displacement vector.
        :rtype: sympy.Matrix, sympy.Matrix
        """
        reaction_vector = sp.zeros(self.dof_per_node*self.num_nodes, 1)
        displacement_vector = sp.zeros(self.dof_per_node*self.num_nodes, 1)
        
        for i, node in enumerate(self.nodes):
            H_symbol, V_symbol, M_symbol = 'H_' + str(i+1), 'V_' + str(i+1), 'M_' + str(i+1)
            u_symbol, v_symbol, o_symbol = 'u_' + str(i+1), 'v_' + str(i+1), 'o_' + str(i+1)
            
            self.symbols.extend([H_symbol, V_symbol, M_symbol, u_symbol, v_symbol, o_symbol])
            
            support_type = node.getSupportType()
            
            if support_type == 'horizontal_roller':
                reaction_vector[self.dof_per_node * i + 1] = V_symbol  

            elif support_type == 'vertical_roller':
                reaction_vector[self.dof_per_node * i] = H_symbol  

            elif support_type == 'pinned' or support_type == 'double_roller':
                reaction_vector[self.dof_per_node * i] = H_symbol  
                reaction_vector[self.dof_per_node * i + 1] = V_symbol  

            elif support_type == 'fixed':
                reaction_vector[self.dof_per_node * i] = H_symbol  
                reaction_vector[self.dof_per_node * i + 1] = V_symbol  
                reaction_vector[self.dof_per_node * i + 2] = M_symbol  

            prescripted_displacements = node.getPrescribedDisplacements()
            
            if prescripted_displacements[0] is not None:
                displacement_vector[self.dof_per_node * i] = prescripted_displacements[0]
            else:
                displacement_vector[self.dof_per_node * i] = u_symbol  
   
            if prescripted_displacements[1] is not None:
                displacement_vector[self.dof_per_node * i + 1] = prescripted_displacements[1]        
            else:
                displacement_vector[self.dof_per_node * i + 1] = v_symbol
                
            if prescripted_displacements[2] is not None:
                displacement_vector[self.dof_per_node * i + 2] = prescripted_displacements[2]
            else:
                displacement_vector[self.dof_per_node * i + 2] = o_symbol

                         
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
    
    def setNodesDisplacementsAndForces(self):
        """
        Sets the displacements and forces of the nodes based on the system solution.
        """
        solution = self.getSolution()
        values = list(solution.values())  
        for i in range(1, self.num_nodes + 1):
            index = (i - 1) * 4
            H, V, u, v = values[index:index+4]
            self.nodes[i - 1].setDisplacement(u, v)
            self.nodes[i - 1].setTotalForces(H, V)
            
    def updateNodesPositions(self):
        """
        Updates all nodes positions based after their displacements.
        """
        for node in self.nodes:
            node.updatePosition()
    
    def setBarsStressesAndNormals(self):
        """
        Sets the normals and stresses of the bars based on the system solution.
        """
        for bar in self.bars:
            bar.setBarNormalAndStress()
            
    def getBarsStressesAndNormals(self):
        """
        Gets the stresses and normals of the bars.
        """
        keys, values = [], []
        for i, bar in enumerate(self.bars):
            keys.append(f'N_{i+1}')
            keys.append(f'sigma_{i+1}')
            values.append(bar.getBarNormal())
            values.append(bar.getBarStress())
        infos = dict(zip(keys, values))
        return infos
