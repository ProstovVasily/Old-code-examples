from tree import decision_tree
from math import exp, log


class boosting:
	def __init__(self, max_depth = 20, number_of_trees = 1, learning_rate = 1):
		self.max_depth = max_depth
		self.n = number_of_trees
		self.learning_rate = learning_rate		
		self.trees = []

	def fit(self, X, y):
		f = [0.0 for x in X]
		for i in range(0, self.n):
			tree = decision_tree (self.max_depth)
			res = residual (f, y)
			tree.fit (X, res)
			self.trees.append(tree)
			f = new_f(f, tree, X, self.learning_rate)
            
	def classify(self, X):
		labels = []
		for x in X:
			tmp = self.f_x(x)
			if tmp > 0.5:
				labels.append (1)
			else:
				labels.append (0)
			
		return labels
    
	def f_x(self, x):
		res = 0
		for tree in self.trees:
			res += self.learning_rate * tree.classify([x])[0]
		return res


def residual(f, y):
	res = []
	for i in range(len(y)):
		res.append(2.0*y[i]/(1+exp(2*y[i]*f[i])))
	return res
		

def new_f(f, tree, X, learning_rate):
	nf = []
	for i in range(len(X)):
		nf.append(f[i] + learning_rate * tree.classify([X[i]])[0])
		
	return nf
