# solution = {
#     'H_1': -20000.0000000000, 'V_1': -31000.0000000000, 'u_1': 0, 'v_1': 0,
#     'H_2': 0, 'V_2': 999.999999999942, 'u_2': -0.000140066315204849, 'v_2': 0,
#     'H_3': 0, 'V_3': 0, 'u_3': 0.00179579866629321, 'v_3': 0.000770139790283814,
#     'H_4': 0, 'V_4': 0, 'u_4': -7.00331576024243e-5, 'v_4': 0.000546410462377163
# }

# values = list(solution.values())

# for i in range(1, 5):  # Iterating from 1 to 4
#     index = (i - 1) * 4
#     print(f'H_{i}: {values[index]}, V_{i}: {values[index + 1]}, u_{i}: {values[index + 2]}, v_{i}: {values[index + 3]}')


class test:
	def __init__(self, a):
		self.a = self.calculate(a)
	
	def calculate(self, value):
		return value/2

a = 2
object_test = test(a)

object_test.calculate
object_test.calculate
object_test.calculate
object_test.calculate

print(object_test.a)