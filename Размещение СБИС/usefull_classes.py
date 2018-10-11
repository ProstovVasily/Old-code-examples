


class node:
    def __init__(self):
        self.pos = [0,0]
        self.size = [0,0]
        self.fixed = False
        self.placed = False
        self.edges = []
        self.square = 0
        self.in_block = False
        
    def __str__(self):
        print self.pos, self.size, self.fixed, len(self.edges)
        return ''

class row:
    def __init__(self):
        self.coordx = 0
        self.coordy = 0
        self.height = 0
        self.shift  = 0
        self.numsites = 0
