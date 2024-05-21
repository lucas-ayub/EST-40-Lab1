
for bar in bars:
    print(bar.getForceVector())
    print(bar.rotation_matrix_6x6)
    print('\n\n\n\n')
    print(bar.rotation_matrix_6x6 @ bar.getForceVector())
