import sympy as sp

# Definindo a função lambda usando funções de sympy
f_lambda = lambda x: sp.sin(x) + sp.cos(x)

# Criando a variável simbólica
x = sp.symbols('x')

# Aplicando a função lambda à variável simbólica
f_applied = f_lambda(x)

# A expressão simbólica resultante
print("Expressão simbólica resultante:")
print(f_applied)
