import numpy as np

# читаем бинарный файл как комплексные числа с float64
data = np.fromfile("./data/complex/1/jg", dtype=np.int32)

for i, val in enumerate(data):
    print(i, val)
