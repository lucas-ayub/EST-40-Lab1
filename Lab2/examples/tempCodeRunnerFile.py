for node in nodes:
    positions = node.getDisplacement()
    forces = node.getGlobalForces()
    i = node.index + 1  # Ajustar para exibição com base em 1
    print(f'u_{i} = {positions[0]}, v_{i} = {positions[1]}, theta_{i} = {positions[2]}')
    print(f'N_{i} = {forces[0]}, V_{i} = {forces[1]}, M_{i} = {forces[2]}')
