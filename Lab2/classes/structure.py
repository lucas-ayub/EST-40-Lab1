import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from node import Node
from bar import Bar

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
        self.f_star = sp.Matrix(self.f) + self.r
        self.solution = self.calculateSolution()
        self.calculateNodalParams()
        self.calculateBarParams()
        self.updateNodesPositions()
        

    def calculateStiffnessMatrix(self):
        """
        Calculates the global stiffness matrix of the system.

        :return: The global stiffness matrix.
        :rtype: numpy.ndarray
        """
        global_size = self.num_nodes * self.dof_per_node
        global_matrix = np.zeros((global_size, global_size))

        for bar in self.bars:
            bar_matrix = bar.getStiffnessMatrix()
            
            left_node_index = bar.left_node.index * self.dof_per_node
            right_node_index = bar.right_node.index * self.dof_per_node
            
            global_matrix[left_node_index:left_node_index+self.dof_per_node, left_node_index:left_node_index+self.dof_per_node] += bar_matrix[:self.dof_per_node, :self.dof_per_node]
            global_matrix[left_node_index:left_node_index+self.dof_per_node, right_node_index:right_node_index+self.dof_per_node] += bar_matrix[:self.dof_per_node, self.dof_per_node:]
            global_matrix[right_node_index:right_node_index+self.dof_per_node, left_node_index:left_node_index+self.dof_per_node] += bar_matrix[self.dof_per_node:, :self.dof_per_node]
            global_matrix[right_node_index:right_node_index+self.dof_per_node, right_node_index:right_node_index+self.dof_per_node] += bar_matrix[self.dof_per_node:, self.dof_per_node:]
            
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

        for bar in self.bars:
            bar_global_force_vector = bar.getForceVector()
            
            left_node_index = bar.left_node.index * self.dof_per_node
            right_node_index = bar.right_node.index * self.dof_per_node
            
            force_vector[left_node_index:left_node_index+self.dof_per_node] += bar_global_force_vector[:self.dof_per_node]
            force_vector[right_node_index:right_node_index+self.dof_per_node] += bar_global_force_vector[self.dof_per_node:]
            
        global_concentrated_forces = np.zeros(self.dof_per_node * len(self.nodes))

        for node in self.nodes:
            node_global_concentrated_forces = node.getGlobalForces()
            i = node.index * 3
            global_concentrated_forces[i:i+3] = node_global_concentrated_forces

        total_force_vector = force_vector + global_concentrated_forces
        
        return total_force_vector
            
    def calculateReactionAndDisplacementVectors(self):
        """
        Calculates the reaction force vector and the displacement vector of the system.

        :return: The reaction force vector and the displacement vector.
        :rtype: sympy.Matrix, sympy.Matrix
        """
        reaction_vector = sp.zeros(self.dof_per_node*self.num_nodes, 1)
        displacement_vector = sp.zeros(self.dof_per_node*self.num_nodes, 1)
        
        for node in self.nodes:
            index = node.index * self.dof_per_node
            N_symbol, V_symbol, M_symbol = f'N_{index+1}', f'V_{index+2}', f'M_{index+3}'
            u_symbol, v_symbol, theta_symbol = f'u_{index+1}', f'v_{index+2}', f'theta_{index+3}'
            
            self.symbols.extend([N_symbol, V_symbol, M_symbol, u_symbol, v_symbol, theta_symbol])
            
            support_type = node.getSupportType()
            prescribed_displacements = node.getPrescribedDisplacements()
            
            if support_type == 'horizontal_roller':
                reaction_vector[index + 1] = V_symbol  
                self.symbols.extend([V_symbol])

            elif support_type == 'vertical_roller':
                reaction_vector[index] = N_symbol  
                self.symbols.extend([N_symbol])

            elif support_type == 'double_roller':
                reaction_vector[index] = N_symbol  
                reaction_vector[index + 1] = V_symbol  
                self.symbols.extend([N_symbol, V_symbol])

            elif support_type == 'pinned':
                reaction_vector[index] = N_symbol  
                reaction_vector[index + 1] = V_symbol  
                self.symbols.extend([N_symbol, V_symbol])

            elif support_type == 'fixed':
                reaction_vector[index] = N_symbol  
                reaction_vector[index + 1] = V_symbol  
                reaction_vector[index + 2] = M_symbol 
                self.symbols.extend([N_symbol, V_symbol, M_symbol])

            
            if prescribed_displacements[0] is not None:
                displacement_vector[index] = prescribed_displacements[0]
            else:
                displacement_vector[index] = u_symbol  
                self.symbols.append(u_symbol)
   
            if prescribed_displacements[1] is not None:
                displacement_vector[index + 1] = prescribed_displacements[1]        
            else:
                displacement_vector[index + 1] = v_symbol
                self.symbols.append(v_symbol)
                
            if prescribed_displacements[2] is not None:
                displacement_vector[index + 2] = prescribed_displacements[2]
            else:
                displacement_vector[index + 2] = theta_symbol
                self.symbols.append(theta_symbol)

                         
        self.symbols = sp.symbols(self.symbols) 
            
        return reaction_vector, displacement_vector

    def sympySolution(self):
        """
        Returns the sympy solution to the system.
        
        :return: Symply dict with the variables and their values
        :rtype: dict
        """
        f = sp.Matrix(self.f_star)
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

        solution = {str(k): v for k, v in solution.items()}
        
        return solution
    
    def getSolution(self):
        """
        Returns the system solution.
        
        :return: The solution vector.
        :rtype: dict
        """
        return self.solution
    
    def calculateNodalParams(self):
        """
        Sets the displacements and forces of the nodes based on the system solution.
        """
        solution = self.getSolution()
        
        solution = {str(k): v for k, v in solution.items()}
        
        for node in self.nodes:
            index = node.index * self.dof_per_node

            N = solution.get(f'N_{index+1}', 0)
            V = solution.get(f'V_{index+2}', 0)
            M = solution.get(f'M_{index+3}', 0)

            node.setFinalForces(N, V, M)

            u = solution.get(f'u_{index+1}', 0)
            v = solution.get(f'v_{index+2}', 0)
            theta = solution.get(f'theta_{index+3}', 0)

            node.setDisplacement(u, v, theta)

                           
    def updateNodesPositions(self):
        """
        Updates all nodes positions based after their displacements.
        """
        for node in self.nodes:
            node.updatePosition()
    
    def calculateBarParams(self):
        """
        Sets the normals, stresses, shear forces and bending moments of the bars based on the system solution.
        """
        for bar in self.bars:
            bar.calculateBarNormalAndStress()
            bar.calculateBarShearForceAndBendingMoment()
            
    def printBarParameters(self):
        """
        Prints the bar parameters.
        """
        for bar in self.bars:
            bar.calculateBarShearForceAndBendingMoment() 
            print(f'Bar {bar.index}:')
            print(f'N = {bar.N}, V = {bar.V}, M = {bar.getBendingMomentumExpression()}')
            print('\n')
    

    def printNodalParameters(self):
        """
        Prints the nodal parameters.
        """
        for node in self.nodes:
            positions = node.getDisplacement()
            forces = node.getGlobalForces()
            i = 3*node.index + 1  
            print(f'node {node.index+1}:')
            print(f'u_{i} = {positions[0]}, v_{i+1} = {positions[1]}, theta_{i+2} = {positions[2]}')
            print('\n')
            