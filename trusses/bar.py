import numpy as np
import node.py

class Bar:
    def __init__(self, left_node, right_node, E, A):
        self.left_node = left_node
        self.right_node = right_node
        self.E = E
        self.A = A
        self.length = np.linalg.norm(self.right_node.position - self.left_node.position)    
        

        