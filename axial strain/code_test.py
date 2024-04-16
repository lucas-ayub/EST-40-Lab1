import numpy as np

def matriz_a_partir_de_matrizes(matrizes):
    # Determinar o tamanho da matriz de saída
    tamanho_saida = len(matrizes) + 1
    
    # Inicializar a matriz de saída com zeros
    matriz_saida = np.zeros((tamanho_saida, tamanho_saida))
    
    # Preencher a diagonal principal com as matrizes de entrada
    for i, matriz in enumerate(matrizes):
        matriz_saida[i:i+2, i:i+2] += matriz
    
    return matriz_saida

# Exemplo de uso
input_matrices = [np.array([[1, -1], [-1, 1]]), np.array([[1, -1], [-1, 1]])]


a = np.array( [ [56, -56], [-56, 56]] )
b = np.array( [ [28, -28], [-28, 28]] )
c = np.array( [ [56, -56], [-56, 56]] )
input_matrices = [a, b, c]
output_matrix = matriz_a_partir_de_matrizes(input_matrices)
print(output_matrix)
