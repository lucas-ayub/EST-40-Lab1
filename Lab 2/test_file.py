import numpy as np

def create_global_stiffness_matrix(bars, num_nodes):
    """
    Cria uma matriz de rigidez global para uma estrutura com várias barras.

    :param bars: Lista de matrizes de rigidez locais de cada barra.
    :param num_nodes: Número total de nós na estrutura.
    :return: Matriz de rigidez global.
    """
    dof_per_node = 3  # Graus de liberdade por nó
    global_size = num_nodes * dof_per_node
    K_global = np.zeros((global_size, global_size))

    # Adicionar cada barra à matriz de rigidez global
    for i, bar_matrix in enumerate(bars):
        start_index = i * dof_per_node
        end_index = start_index + 2 * dof_per_node
        K_global[start_index:end_index, start_index:end_index] += bar_matrix

    return K_global

# Matrizes de rigidez locais das barras fornecidas
K1_global = np.array([
    [2400, 0, 0, -2400, 0, 0],
    [0, 0.96, 240, 0, -0.96, 240],
    [0, 240, 80000, 0, -240, 40000],
    [-2400, 0, 0, 2400, 0, 0],
    [0, -0.96, -240, 0, 0.96, -240],
    [0, 240, 40000, 0, -240, 80000]
])

K2_global = np.array([
    [3074.7, -2300.3, 576, -3074.7, 2300.3, 576],
    [-2300.3, 1733, 768, 2300.3, -1733, 768],
    [576, 768, 160000, -576, -768, 80000],
    [-3074.7, 2300.3, -576, 3074.7, -2300.3, -576],
    [2300.3, -1733, -768, -2300.3, 1733, -768],
    [576, 768, 80000, -576, -768, 160000]
])

# Número total de nós na estrutura
num_nodes = 3  # Assume que temos 3 nós para 2 barras

# Calcula a matriz de rigidez global
K_global = create_global_stiffness_matrix([K1_global, K2_global], num_nodes)

# Imprime a matriz de rigidez global calculada
print("Matriz de rigidez global calculada para a estrutura:")
print(K_global)
