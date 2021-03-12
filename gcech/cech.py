from itertools import combinations
import math

# code credit to https://stackoverflow.com/questions/55816902/finding-the-intersection-of-two-circles
def get_intersections(x0, y0, r0, x1, y1, r1):
    # circle 1: (x0, y0), radius r0
    # circle 2: (x1, y1), radius r1

    d=math.sqrt((x1-x0)**2 + (y1-y0)**2)
    # non intersecting
    if d > r0 + r1 :
        return []
    # One circle within other
    if d < abs(r0-r1):
        return []
    # coincident circles
    if d == 0 and r0 == r1:
        return []
    else:
        a=(r0**2-r1**2+d**2)/(2*d)
        h=math.sqrt(r0**2-a**2)
        x2=x0+a*(x1-x0)/d   
        y2=y0+a*(y1-y0)/d   
        x3=x2+h*(y1-y0)/d     
        y3=y2-h*(x1-x0)/d 

        x4=x2-h*(y1-y0)/d
        y4=y2+h*(x1-x0)/d
        
        return [[x3, y3], [x4, y4]]

class CechComplexBase(object):
    """
    Base class for Cech complex
    """
    def __init__(self) -> None:
        pass

    def cech(self, positions, radius, maxK=None):
        """
        Compute Cech complex given position and radius
        """
        self.positions = positions
        self.radius = radius
        self.maxK = maxK
        
        self.N = len(positions)
        self.neighbor = [[] for i in range(self.N)]
        self.S = [self.sim0(), self.sim1()]
        self.simk()
        return self.S

    def sim0(self):
        return [[i] for i in range(self.N)]

    def sim1(self):
        S1 = []
        for i in range(self.N):
            self.neighbor[i] = set()
        for i in range(self.N):
            for j in range(i+1, self.N):
                if self.intersects(i, j):
                    self.neighbor[i].add(j)
                    self.neighbor[j].add(i)
                    S1.append([i, j])
        return S1

    def simk(self):
        pass

    def isCellInsideCell(self, i, j):
        return self.dist(i, j) <= abs(self.radius[i]-self.radius[j])


    def isPointInsideCell(self, p, cellid):
        x1, y1 = self.positions[cellid]
        x2, y2 = p
        return (x1-x2)**2 + (y1-y2)**2 <= self.radius[cellid]**2


    def intersects(self, i, j, d=None):
        if d is None:
            d = self.dist(i, j)
        
        x1, y1 = self.positions[i]
        x2, y2 = self.positions[j]
        
        if d > self.radius[i] + self.radius[j]:
            return False
        
        return True

    def dist(self, i, j):
        x1, y1 = self.positions[i]
        x2, y2 = self.positions[j]
        return math.sqrt((x1-x2)**2+(y1-y2)**2)
    
    def intersection(self, i,j):
        x1, y1 = self.positions[i]
        x2, y2 = self.positions[j]
        return get_intersections(x1, y1, self.radius[i], x2, y2, self.radius[j])

class CechComplexLE(CechComplexBase):
    """
    Class for algorithm in "Construction of the generalized Cech complex" 
    url: https://arxiv.org/abs/1409.8225
    """
    def simk(self):
        def _simk():
            Sk = []
            for cand in self.getCandidates(k):
                if self.verify(cand):
                    Sk.append(cand)
            return Sk
            
        verified = {}
        k = 2
        while True:
            if (self.maxK and k > self.maxK):
                break
            Sk = _simk()
            if len(Sk) == 0:
                break
            else:
                self.S.append(Sk)
                k += 1
            
                
    def getCandidates(self, k):
        for i in range(self.N):
            for cand in combinations(self.neighbor[i], k):
                ns = [i, *cand]
                yield ns

    
    def verify(self, cells):
        smallest = cells[0]
        for c in cells:
            if self.radius[c] < self.radius[smallest]:
                smallest = c
        
        inside = True
        for c in cells:
            if c == smallest or self.isCellInsideCell(c, smallest):
                continue
            else:
                inside = False
        
        if(inside):
            return True
        else:
            inter = []
            for i,j in combinations(cells, 2):
                x = self.intersection(i,j)
                if x:
                    inter.extend([(p,i,j) for p in x])
            
            for p,i,j in inter:
                exists = True
                for cellid in cells:
                    if cellid  == i or cellid == j:
                        continue
                    elif not self.isPointInsideCell(p, cellid):
                        exists = False
                        break
                if exists:
                    return True
            
            return False

class CechComplex(CechComplexBase):
    """
    Class for improved algorithm to construct Cech complex.
    """
    def simk(self):
        self.intersections = {}
        verified = {}
        k = 2
        
        def _simk():
            Sk = []
            for cand in self.getCandidates(k):
                ns_sig = frozenset(cand)
                if ns_sig not in verified:
                    verified[ns_sig], inter = self.verify(cand)
                    if verified[ns_sig]:
                        self.intersections[ns_sig] = inter
                        Sk.append(cand)
            return Sk
        
        while True:
            if (self.maxK and k > self.maxK):
                break
            #print(f"{k-1}:{len(S[k-1])}")
            Sk = _simk()
            if len(Sk) != 0:
                self.S.append(Sk)
                k += 1
            else:
                break

    def getCandidates(self, k):
        for s in self.S[k-1]:
            cand = set.intersection(*[self.neighbor[i] for i in s]).difference(set(s))
            for c in cand:
                ns = [*s, c]
                yield ns
    
    def verify(self, cells):
        smallest = cells[0]
        for c in cells:
            if self.radius[c] < self.radius[smallest]:
                smallest = c
        
        inside = True
        for c in cells:
            if c == smallest or self.isCellInsideCell(c, smallest):
                continue
            else:
                inside = False
        
        if(inside):
            return True, self.positions[smallest]
        else:
            new_cand = cells[-1]
            pre_set = cells[:-1]
            if len(pre_set) == 2:
                x = self.intersection(cells[0], cells[1])
                for p in x:
                    if self.isPointInsideCell(p, new_cand):
                        return True, p
            elif len(pre_set) > 2 and self.isPointInsideCell(self.intersections[frozenset(pre_set)], new_cand):
                return True, self.intersections[frozenset(pre_set)]
            inter = []
            i = len(cells)-1
            for j in range(0, i):
                x = self.intersection(cells[i], cells[j])
                if x:
                    for p in x:
                        exists = True
                        for k, cellid in enumerate(cells):
                            if k == i or k == j:
                                continue
                            if not self.isPointInsideCell(p, cellid):
                                exists = False
                                break
                        if exists:
                            return True, p
            return False, []