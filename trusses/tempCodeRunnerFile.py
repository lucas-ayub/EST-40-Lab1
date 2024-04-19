for i, bar in enumerate(list_of_bars):
    print(f'bar_{i+1}:\n')
    print(f'N_{i+1} = {bar.getBarNormal()}\n')
    print(f'sigma_{i+1} = {bar.getBarStress()}\n')