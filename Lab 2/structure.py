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
        global_matrix = np.zeros((3 * self.num_nodes, 3 * self.num_nodes))

        for bar in self.bars:
            i = self.nodes.index(bar.left_node)
            j = self.nodes.index(bar.right_node)
            
            local_matrix = bar.getStiffnessMatrix()
            
            global_matrix[2*i:2*(i+1), 2*i:2*(i+1)] += local_matrix[:2, :2]
            global_matrix[2*i:2*(i+1), 2*j:2*(j+1)] += local_matrix[:2, 2:]
            global_matrix[2*j:2*(j+1), 2*i:2*(i+1)] += local_matrix[2:, :2]
            global_matrix[2*j:2*(j+1), 2*j:2*(j+1)] += local_matrix[2:, 2:]

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
        force_vector = np.zeros(3 * len(self.nodes))

        for bar in self.bars:
            bar.updateNodeForces()
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

    def plotStructure(self, displacement_scale = 1.0):
        """
        Plots the initial and final configurations of the structure.
        
        :param displacement_scale: Scale factor for displacements, default is 1.0.
        """
        _, ax = plt.subplots()
        
        initial_x = [position[0] for position in self.nodes_initial_positions]
        initial_y = [position[1] for position in self.nodes_initial_positions]
        ax.plot(initial_x, initial_y, 'ko', label='Initial')

        for bar in self.bars:
            left_initial_pos = bar.left_node.getPosition()
            right_initial_pos = bar.right_node.getPosition()
            ax.plot([left_initial_pos[0], right_initial_pos[0]], [left_initial_pos[1], right_initial_pos[1]], 'k-')

        for bar in self.bars:
            left_initial_pos = bar.left_node.getPosition()
            right_initial_pos = bar.right_node.getPosition()
            left_final_pos = None
            right_final_pos = None
            for idx, pos in enumerate(self.nodes_final_positions):
                if np.array_equal(pos, left_initial_pos):
                    left_final_pos = pos + displacement_scale * bar.left_node.getDisplacement()
                if np.array_equal(pos, right_initial_pos):
                    right_final_pos = pos + displacement_scale * bar.right_node.getDisplacement()
                if left_final_pos is not None and right_final_pos is not None:
                    break
            ax.plot([left_final_pos[0], right_final_pos[0]], [left_final_pos[1], right_final_pos[1]], 'r--')
        
        final_x = [position[0] + displacement_scale * node.getDisplacement()[0] for position, node in zip(self.nodes_final_positions, self.nodes)]
        final_y = [position[1] + displacement_scale * node.getDisplacement()[1] for position, node in zip(self.nodes_final_positions, self.nodes)]
        ax.plot(final_x, final_y, 'ro', label='Final')

        ax.legend()
        ax.set_aspect('equal', 'box')
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        plt.title(f'Structure  (Displacement Factor: {displacement_scale})')
        plt.grid(True)
        plt.show()