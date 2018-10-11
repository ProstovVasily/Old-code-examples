from usefull_classes import *

from math import sqrt
from scipy import ndimage
from scipy import misc
import matplotlib.pyplot as plt
import numpy as np


class scheme:
    
    def __init__(self):
        self.nodes = {}
        self.min_x = 0
        self.min_y = 0
        self.max_x = 0
        self.max_y = 0
        self.max_size = 0
        self.NumRows = 0
        self.rows = []

        
    def show_scheme(self):
        #here would be visualisation
        #for x in self.nodes:
         #   print x, self.nodes[x], type(x), type(self.nodes[x])
        #GOGOGOGOGOGO
        print len(self.rows), self.NumRows, len(self.nodes)
        size = (self.max_x - self.min_x + 1, self.max_y - self.min_y + 1, 3)
        print self.min_x, self.min_y, self.max_x, self.max_y
        scale = 1.0
        if size[0] > 2000:
            scale = size[0]/1000
            size = (1000, 1000, 3)
        print size, scale
        im = np.zeros(size)
        shift_x = -self.min_x
        shift_y = -self.min_y
        FUcount  = 0
        for p in self.nodes:
            t = self.nodes[p]
            if t.placed == False:
                a = [0,1,0]
            else:
                a = [1,1,1]
            if (t.in_block == False) and (t.placed == False):
                a = [1,0,0]
                
            x0 = min(max(0,(t.pos[0] + shift_x)/scale), size[0]-2)
            x1 = min(max(0,(t.pos[0] + t.size[0] + shift_x)/scale), size[0]-2)
            y0 = min(max(0,(t.pos[1] + shift_y)/scale), size[1]-2)
            y1 = min(max(0,(t.pos[1] + t.size[1] + shift_y)/scale), size[1]-2)
            if (x1 <= x0 or y1 <= y0) and x1 < 990 and y1 < 990 and x0 > 10 and y0 > 10:
                print "FFFF", x0, x1, y0, y1 
                #x1 = x0 + 1
                #y1 = y0 + 1
                FUcount += 1
            if t.fixed:
                im[x0:min(size[0], x1+1), y0:min(size[1], y1+1)] = [0.5, 0, 0.5]
            else:
                im[x0:min(size[0], x1+1), y0:min(size[1], y1+1)] = a
        #im[900:999, 900:999] = [1,1,1]
        #im[0:100, 900:999] = [0,1,0]
        #im[900:999, 0:100] = [1,0,1]
        plt.imshow(im)
        print FUcount
        plt.show()
                
             
    def parsing (self, file_name):
        self.parse_nodes (file_name + ".nodes")
        self.parse_pl (file_name + ".pl")
        self.parse_nets (file_name + ".nets")
        self.parse_scl (file_name + ".scl")
        print len(self.nodes)
    
        
    def parse_nodes(self, file_name):
        f = open(file_name, "r")
        for string in f.readlines():
            if string == '\n':
                continue   
            s = (' '.join(string.split('\t'))).split(' ')
            s = [x for x in s if x != '']
            if s[0] in ['#', 'UCLA', 'NumNodes', 'NumTerminals']:
                continue
            tmp = node()
            tmp.size = [float(s[1]), float(s[2])]
            tmp.square = float(s[1])*float(s[2])
            self.max_size = max (self.max_size, float(s[1]), float(s[2]))
            if len(s) > 3:
                tmp.fixed = True
                tmp.placed = True
            self.nodes[s[0]] = tmp
        f.close()
            
    def parse_pl(self, file_name):
        f = open(file_name, "r")
        for string in f.readlines():
            if string == '\n':
                continue   
            s = (' '.join(string.split('\t'))).split(' ')
            s = [x for x in s if x != '']
            if s[0] in ['#', 'UCLA']:
                continue
            self.nodes[s[0]].pos = [float(s[1]), float(s[2])]
            self.max_x = max(float(s[1]) + self.nodes[s[0]].size[0], self.max_x)
            self.max_y = max(float(s[2]) + self.nodes[s[0]].size[1], self.max_y)
            self.min_x = min(float(s[1]), self.min_x)
            self.min_y = min(float(s[2]), self.min_y)
        f.close()
        
    def parse_nets(self, file_name):
        f = open(file_name, "r")
        string = ' '
        while string != '':
            string = f.readline()
            if string == '':
                break
            if string == '\n':
                continue   
            s = (' '.join(string.split('\t'))).split(' ')
            s = [x for x in s if x != '']
            if s[0] in ['#', 'UCLA', 'NumNets', 'NumPins']:
                continue
            inp = []
            out = []
            if s[0] == 'NetDegree':
                for x in xrange(int(s[2])):
                    tmp = f.readline()
                    s = (' '.join(tmp.split('\t'))).split(' ')
                    s = [x for x in s if x != '']
                    if s[1] == 'I':
                        inp.append([s[0],float(s[3]),float(s[4])])
                    else:
                        out.append([s[0],float(s[3]),float(s[4])])
                for x in inp:
                    for t in out:
                        self.nodes[x[0]].edges.append(t + x[1:3])
                        self.nodes[t[0]].edges.append(x + t[1:3])
        f.close()

    def parse_scl (self, file_name):
        f = open(file_name, "r")
        string = ' '
        while string != '':
            string = f.readline()
            if string == '':
                break
            if string == '\n':
                continue   
            s = (' '.join(string.split('\t'))).split(' ')
            s = [x for x in s if x != '']
            if s[0] == 'NumRows':
                self.NumRows = int(float(s[2]))
                continue
            if s[0] in ['#', 'UCLA']:
                continue
            if s[0] == 'CoreRow':
                crow = row()
                for i in xrange(8):
                    string = f.readline()
                    s = (' '.join(string.split('\t'))).split(' ')
                    s = [x for x in s if x != '']
                    if s[0] == 'Coordinate':
                        crow.coordx = int(float(s[2]))
                    if s[0] == 'Height':
                        crow.height = int(float(s[2]))
                    if s[0] == 'Sitewidth':
                        crow.unused1 = int(float(s[2]))
                    if s[0] == 'Sitespacing':
                        crow.shift = int(float(s[2]))
                    if s[0] == 'Siteorient':
                        crow.unused2 = int(float(s[2]))
                    if s[0] == 'Sitesymmetry':
                        crow.unused3 = int(float(s[2]))
                    if s[0] == 'SubrowOrigin':
                        crow.coordy = int(float(s[2]))
                        crow.numsites = int(float(s[5]))
                    if s[0] == 'End':
                        break
                self.rows.append(crow)
                    
        
        
    def global_placement (self, num_it):
        for i in xrange(num_it):
            print i
            cnt = 0
            for p in self.nodes:
                nod = self.nodes[p]
                if nod.fixed:
                    continue
                Fx = 0
                Fy = 0
                ntsh = False
                for edg in nod.edges:
                    if (abs(self.nodes[edg[0]].pos[0] - nod.pos[0]) < 1) and (abs(self.nodes[edg[0]].pos[1] - nod.pos[1]) < 1):
                        ntsh = True
    
                    Fx += (self.nodes[edg[0]].pos[0] + edg[1] - (nod.pos[0] + edg[3]))*c # + self.nodes[edg[0]].size[0]/2
                    Fy += (self.nodes[edg[0]].pos[1] + edg[2] - (nod.pos[1] + edg[4]))*c # + self.nodes[edg[0]].size[1]/2
                #Fx = Fx - (nod.pos[0])*len(nod.edges)
                #Fy = Fy - (nod.pos[1])*len(nod.edges)
                nod.pos = [nod.pos[0] + Fx/20, nod.pos[1] + Fy/20]
                self.nodes[p] = nod
                if ntsh:
                    cnt += 1
            print cnt
            
    def legalisation(self):
        n = 10 #number of rows in block 
        h = self.rows[0].height
        minx = self.rows[0].coordx
        miny = self.rows[0].coordy
        maxy = self.rows[0].coordy + self.rows[0].numsites 
        maxx = self.rows[-1].coordx + h
        sz = int(self.NumRows/n)
        blocks = [[[] for x in xrange(self.NumRows/n)] for x in xrange(self.NumRows/n)] 
        mx = n*h
        my = int((maxy-miny)/sz)
        for w in self.nodes:
            nd = self.nodes[w]
            x = nd.pos[0]
            y = nd.pos[1]
            if x < minx:
                if not nd.fixed:
                    i = 0
                else:
                    if x + nd.size[0] < minx:
                        continue
                    else:
                        i = 0
            elif x >= maxx:
                if nd.fixed:
                    continue
                else:
                    i = sz - 1
            else:
                i = int(sz*(x - minx)/(maxx - minx))
                
            if y < miny:
                if not nd.fixed:
                    j = 0
                else:
                    if y + nd.size[1] < miny:
                        continue
                    else:
                        j = 0
            elif y > maxy:
                if nd.fixed:
                    continue
                else:
                    j = sz - 1
            else:
                j = int(sz*(y - miny)/(maxy - miny))
                
            if nd.fixed:
                blocks[i][j].insert(0, w)
            else:
                blocks[i][j].append(w)
                self.nodes[w].in_block = True
            #add check for terminal
            if nd.fixed:
                i1 = min(int(sz*(x + nd.size[0] - minx)/(maxx - minx)), sz - 1)
                j1 = min(int(sz*(y + nd.size[1] - miny)/(maxy - miny)), sz - 1)
                i0 = i
                j0 = j
                for i in xrange(i0, i1+1):
                    for j in xrange (j0, j1+1):
                        if (i == i0) and (j == j0):
                            continue
                        blocks[i][j].insert(0, w)                    
        mconc = 0
        cnt = 0
        omfg = 0
        concentrations = np.zeros((sz,sz))
        
        sm = 0
        for ttt in self.nodes:
            sm += self.nodes[ttt].square
        print 'Conc', sm*1.0/((maxx-minx)*(maxy-miny))
        
        
        for i in xrange(sz):
            for j in xrange(sz):
                sm = 0
                for k in blocks[i][j]:
                    if self.nodes[k].fixed:
                        continue
                    sm += self.nodes[k].square
                conc = float(sm)/float(mx*my) #actually here would be sun squares/total zone square
                concentrations[i][j] = conc
                if conc > mconc:
                    mconc = conc
                    print i, j
                if conc > 1:
                    cnt += 1
                if conc > 10:
                    omfg += 1
        i = 0
        j = 0                
        while (i < sz):
            j = 0
            redirect = False
            while (j < sz):
                redirect = False
                conc = concentrations[i][j]
                #here goes placement
                global_x = self.rows[i*n].coordx
                global_y = j*my + self.rows[i*n].coordy
                pl_mat = np.zeros((mx, my))#change to more precice values
                cur_x = 0
                cur_y = 0
                #place terminals
                k = 0
                if len(blocks[i][j]) == 0:
                    j += 1
                    continue
                if self.nodes[blocks[i][j][-1]].placed:
                    j += 1
                    continue
                for w in blocks[i][j]:
                    nd = self.nodes[w]
                    if not nd.placed:
                        break
                    k += 1
                    pl_mat[max(nd.pos[0]-global_x, 0):min(nd.pos[0]+nd.size[0]+1 - global_x, mx), 
                            max(nd.pos[1]-global_y, 0):min(nd.pos[1]+nd.size[1]+1 - global_y, my)] = 1
                
                #blocks[i][j] = blocks[i][j][k:]
                #place other
                #for w in blocks[i][j]:
                for k in xrange(k, len(blocks[i][j])):
                    w = blocks[i][j][k]
                    nd = self.nodes[w]
                    if nd.fixed:
                        print "WWWWWWWWWWW"
                    x = nd.size[0]
                    y = nd.size[1]
                    stop = False
                    if (cur_x > mx) or (cur_y > my):
                        break
                    for cur_x in xrange(cur_x, mx, h):
                        for cur_y in xrange(cur_y, my):
                            if cur_x + x >= mx:
                                break
                            if cur_y + y >= my:
                                break
                            if sum(sum(pl_mat[cur_x:cur_x+x+1, cur_y:cur_y+y+1])) < 0.5:
                                stop = True
                                pl_mat[cur_x:cur_x+x+1, cur_y:cur_y+y+1] = 1
                                nd.pos = [cur_x + global_x, cur_y + global_y]
                                nd.placed = True
                                cur_y += int(y) + self.rows[i].shift
                                self.nodes[w] = nd
                                break
                        if stop:
                            break
                        cur_y = 0
                    if not stop:
                        #change concentration status
                        concentrations[i][j] = 1 #or smth else
                        #move extra nodes to other blocks
                        first_guess = max(int((sqrt(conc) + 1)/2), 1)
                        p0 = max(0, i - first_guess)
                        p1 = min(sz, i + first_guess + 1)
                        q0 = max(0, j - first_guess)
                        q1 = min(sz, j + first_guess + 1)
                        tmp = []
                        while (len(tmp) == 0):
                            p0 = max(0, i - first_guess)
                            p1 = min(sz, i + first_guess + 1)
                            q0 = max(0, j - first_guess)
                            q1 = min(sz, j + first_guess + 1)
                            if first_guess > 2*sz:
                                print "redistribution impossible"
                                break
                            for p in xrange(p0,p1):
                                for q in xrange(q0,q1):
                                    if concentrations[p][q] > 0.99:
                                        continue
                                    if (i == p) and (j == q):
                                        print "WTF", concentrations[i][j]
                                        continue 
                                    tmp.append([p,q])
                            first_guess *= 2
                        
                        if len(tmp) == 0:
                            break
                        extra = blocks[i][j][k:]
                        lt = len(tmp)
                        le = len(extra)
                        s = le/lt + 1 * (le % lt != 0)
                        for t in xrange(0, lt):
                            p = tmp[t][0]
                            q = tmp[t][1]
                            if t*s >= le:
                                break
                            blocks[p][q] += extra[t*s:min((t+1)*s, le)]
                        #remove extra from current block
                        blocks[i][j] = blocks[i][j][:k]
                        if (i > tmp[0][0]) or ((i == tmp[0][0]) and (j > tmp[0][1])):
                            i = tmp[0][0]
                            j = tmp[0][1]
                        #i = 0
                        #j = 0
                            redirect = True
                        break
                    #add spaces between
                if not redirect:
                    j += 1
            if not redirect:
                i += 1
        sm = 0
        for ttt in self.nodes:
            sm += self.nodes[ttt].square
        print 'Conc', sm*1.0/((maxx-minx)*(maxy-miny))
        print 'Res: ', mconc, cnt, omfg

            
a = scheme()
a.parsing('test1/bigblue1')  
print "Parsing Done!"    
a.global_placement(0)
print "Global placement Done"
a.legalisation()
a.show_scheme()

    
    
  
    

 
