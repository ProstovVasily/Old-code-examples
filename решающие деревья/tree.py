from math import log


def log2(x):
	if x <=0 : return 0
	return log(x, 2)
	
def entropy (y):	
	n = len(y)
	if n == 0: return 0
	p1 = sum(y)/(n*1.0)
	p2 = 1 - p1
	H = -(p1*log2(p1) + p2*log2(p2))
	return H

class tree:
	def __init__ (self, feature = 0, val = 0, leaf = False, label = None):
		self.border = val
		self.feature = feature
		self.left = None
		self.right = None
		self.leaf = leaf
		self.label = label
		
		def is_leaf(self): 
			return self.leaf
			
class decision_tree:
	def __init__ (self, max_depth = 20):
		self.max_depth = max_depth
		self.root = None
		
	def fit (self, X, y):
		self.root = self.build_tree (X, y, self.max_depth)
				
		
	def classify (self, X):
		labels = []
		
		for x in X:
			label = None
			node = self.root
			while label == None:
				if node.leaf:
					label = node.label
				else:
					feat = node.feature
					val = node.border
					if x[feat] >= val:
						node = node.right 
					else:
						node = node.left						
			labels.append(label)
		return labels
		
	def information_gain(self, X, y, f, val):
		set1 = []
		set2 = []
		n = len(y)
		for i in range(0, n):
			if X[i][f] >= val:
				set2.append(y[i])
			else:
				set1.append(y[i])
		if n <= 0: print n
		return entropy(y) - entropy(set1)* len(set1)/(n*1.0) - entropy(set2)* len(set2)/(n*1.0)
		
		
	def best_divide (self, X, y):
		maxi = -10000
		bval = 0
		bf = 0
		n = len(X[0])
		for f in range(0, n):
			for val in range(0, 100):
				curr = self.information_gain (X, y, f, val)
				if curr > maxi:
					maxi = curr 
					bf = f
					bval = val
		return [bf, bval] 
		
	def build_tree (self, X, y, depth):
		if (depth <= 0) or (len(y) == sum(y)) or (sum(y) == 0):
			node = tree(0, 0, True, 0)
			if sum(y)*2 >= len(y):
				node.label = 1
			return node
		[f, val] = self.best_divide (X, y)
		X1 = []
		X2 = []
		y1 = []
		y2 = []
		n = len(y)
		for i in range(0, n):
			if X[i][f] >= val:
				X2.append (X[i])
				y2.append (y[i])
			else:
				X1.append (X[i])
				y1.append (y[i])
		node = tree(f, val)
		node.left = self.build_tree (X1, y1, depth-1)
		node.right = self.build_tree (X2, y2, depth-1)
		return node
		
		
		
		
		
		
		
			
		
