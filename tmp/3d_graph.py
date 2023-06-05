import matplotlib.pyplot as plt
import os.path
import numpy as np
from sys import exit
from sys import argv
import json

if (not(os.path.exists("result.txt"))):
    print("The file does not exist")
    exit()

if (len(argv) != 3):
    print("Incorrect launch")
    exit()

file = open("result.txt")
data = json.loads(file.read())

var = int(argv[2])
label = ""
zlabel = ""
arr = []
a = float(data["a"])
b = float(data["b"])
c = float(data["c"])
d = float(data["d"])
n = int(data["n"])
m = int(data["m"])
test = bool(data["test"])

h = 0.0
k = 0.0

if var == 1:
	arr = data["arr_u"][0];
	if test:
		label = "Численное решение тестовой задачи"
	else:
		label = "Решение основной задачи на основной сетке"
	zlabel = "u1"
	h = (b-a)/n
	k = (d-c)/m
	print("var == 1")
elif var == 2:
	arr = data["arr_u"][1];
	if test:
		label = "Точное решение тестовой задачи"
	else:
		label = "Решение основной задачи на контрольной сетке"
	zlabel = "u2"
	if test:
		h = (b-a)/n
		k = (d-c)/m
	else:
		h = (b-a)/(2*n)
		k = (d-c)/(2*m)
	print("var == 2")
elif var == 3:
	arr = data["arr_err"]
	if test:
		label = "Погрешность тестовой задачи"
	else:
		label = "Точность основной задачи"
	zlabel = "|u1-u2|"
	h = (b-a)/n
	k = (d-c)/m
	print("var == 3")
else:
	print("Incorrect launch")
	exit()

#with open("result.txt") as res:
    #all_result = [row.strip() for row in res]

#result = []
#for i, row in enumerate(all_result):
    #result.append(row.split(";"))
    #result[i].pop()

#numberVariant = 1
m = len(arr)
n = len(arr[0])

new_result = np.array(arr)

x = []
for i in range(n):
	x.append(a + i * h)

y = []
for j in range(m):
	y.append(c + j * k)

xgrid, ygrid = np.meshgrid(x, y)

fig = plt.figure("")
ax = plt.axes(projection ='3d')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel(zlabel)
ax.set_title(label)

#ax.plot_wireframe(xgrid, ygrid, new_result)
ax.plot_surface(ygrid, xgrid, new_result, cmap='viridis', 
                       edgecolor='black', shade = False, linewidth = 0.3)

def y_fmt(x, y):
    return '${:2.1e}'.format(x).replace('e', '\\cdot 10^{') + '}$'

ax.zaxis.set_major_formatter(y_fmt)
plt.show()
