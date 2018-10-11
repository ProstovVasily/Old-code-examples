from random import *
import matplotlib.pyplot as plt
import numpy as np
from tree import *
from boosting import *


def generate_set (center, radius, amount):
	return [[center[0] + randint(-radius, +radius), center[1] + randint(-radius, +radius)] for x in range(0, amount)]


def generate_noise (amount):
	return [[randint(0,100), randint(0, 100)] for x in range(0, amount)]

def show_data (X, y):
	x1 = []
	x2 = []
	for i in range(0, len(y)):
		if y[i] == 0:
			x1.append(X[i])
		else:
			x2.append(X[i])
			
	v1 = [t[0] for t in x1]
	v2 = [t[1] for t in x1]
	plt.plot(np.array(v1), np.array(v2), 'ro')

	v1 = [t[0] for t in x2]
	v2 = [t[1] for t in x2]
	plt.plot(np.array(v1), np.array(v2), 'bo')
	
def show_tree(node, xmin, xmax, ymin, ymax):
	if node.leaf == True: 
		return
	elif node.feature == 0:
		plt.axvline (node.border, ymin/100.0, ymax/100.0)
		show_tree(node.left, xmin, node.border, ymin, ymax)
		show_tree(node.right, node.border, xmax, ymin, ymax)
		
	else:
		plt.axhline (node.border, xmin/100.0,  xmax/100.0, color = "red")	
		show_tree(node.left, xmin, xmax, ymin, node.border)
		show_tree(node.right, xmin, xmax, node.border, ymax)
		
	

tst = generate_set ([50, 50], 10, 30);
X = tst
y = [0 for x in tst]
tst = generate_set ([20, 50], 10, 30);
X = X + tst
y = y + [1 for x in tst]
tst = generate_set ([70, 80], 10, 30);
X = X + tst
y = y + [1 for x in tst]
tst = generate_set ([80, 20], 10, 30);
X = X + tst
y = y + [0 for x in tst]
tst = generate_noise(10)
X = X + tst
y = y + [randint(0,1) for x in tst]

obj = decision_tree(4)
obj.fit(X, y)

test = boosting(4)
test.fit(X, y)

print test.classify([[20, 20], [20, 80], [80, 20], [80, 80]])
print obj.classify([[20, 20], [20, 80], [80, 20], [80, 80]])

show_data(X, y)
plt.plot ([0, 0, 100, 100], [0, 100, 0, 100], 'o')
show_tree(obj.root, 0, 100, 0, 100)


plt.show()

