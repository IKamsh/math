from fractions import Fraction
import copy

class RationalMatrix(object):
    '''
    matrices over the field of rational numbers nxm
    
    '''
    
    def __init__(self, matrix=None, n=1, m=1):
        if type(matrix) == tuple or type(matrix) == list:
            self.n = len(matrix)
            if type(matrix[0]) != list and type(matrix[0]) != tuple:
                self.m = 1
                self.matrix = [[matrix[i]] for i in range(len(x))]
            else:
                self.m = len(matrix[0])
                self.matrix = [[Fraction(matrix[i][j]) for j in range(len(matrix[i]))] for i in range(len(matrix))]
        if matrix == None:
            self.n = n
            self.m = m
            self.matrix = [[Fraction(0) for j in range(m)] for i in range(n)]
            
        self.det = None
        
    def __iadd__(self, a):
        return self + a
    
    def __add__(self, a):
        assert type(a) == RationalMatrix
        assert self.n == a.n and self.m == self.m
        return RationalMatrix([[self.matrix[i][j]  + a.matrix[i][j] for j in range(a.m)]
                      for i in range(a.n)])
    
    def __mul__(self, a):
        assert type(a) == RationalMatrix
        assert self.m == a.n
        res = [[0]*self.n for i in range(a.m)]
        for i in range(self.n):
            for j in range(a.m):
                for k in range(a.m):
                    res[i][j] += self.matrix[i][k] * a.matrix[k][j]
        return RationalMatrix(res)
        
    def __imul__(self, a):
        return self * a
        
    def determinant(self):
        '''
        return determinant of the matrix
        '''
        assert self.m == self.n
        if self.det == None:
            self.upper_triang()
        return self.det
    
    def upper_triang(self):
        '''
        returns a copy of the matrix reduced to the upper triangular form
        '''
        triang = copy.deepcopy(self)
        det_sign = 1
        num_row = 0
        for j in range(self.m): #column
            for i in range(num_row, self.n):
                if triang.matrix[i][j] != 0: #choose row with non-zero element
                    if i != num_row:
                        det_sign += -1
                    triang.matrix[num_row], triang.matrix[i] = triang.matrix[i], triang.matrix[num_row]
                    for k in range(num_row + 1, self.n):
                        factor = triang.matrix[k][j] / triang.matrix[num_row][j]
                        triang.add_row(k, num_row, -factor)
                    num_row += 1
                    break
                    
        if self.n == self.m: # determinant computation
            for i in range(self.n):
                det_sign *= triang.matrix[i][i]
            self.det = det_sign
            
        return triang
                
    def rank(self):
        '''
        return rank of the matrix
        '''
        c = 0
        mat = self.upper_triang()
        for i in range(mat.n - 1, -1, -1):
            for j in range(mat.m):
                if mat.matrix[i][j] != 0:
                    return mat.n - c
            c += 1
        return mat.n
    
    def add_row(self, toAddInd, addedInd, factor):
        '''
        add row with index 'addedInd' multiplied by 'factor' to row with index 'toAddInd'
        '''
        for i in range(self.m):
            self.matrix[toAddInd][i] += self.matrix[addedInd][i] * factor
        return self
            
    def __print__(self):
        st = ''
        for i in range(self.n):
            for j in range(self.m):
                st += str(self.matrix[i][j]) + ' '
            st += '\n'
        return st
    
    def __repr__(self):
        st = ''
        for i in range(self.n):
            for j in range(self.m):
                st += str(self.matrix[i][j]) + ' '
            st += '\n'
        st += 'class: Rational Matrix'
        return st

    def solve(self, b):
        '''
        solve system of linear exuations
        Ax = b
        A - is RationalMatrix object
        b - expected to be a list or a tuple of matrix rows number size
        '''
        assert len(b) == self.n
        extended_mat = copy.deepcopy(self)
        for i in range(self.n):
            extended_mat.matrix[i].append(Fraction(b[i]))
        
        extended_mat.m += 1
        
        if extended_mat.rank() == self.rank():
            if extended_mat.rank() < self.m:
                print('System has an infinite number of solutions')
                return None
            extended_mat = extended_mat.upper_triang()
            x = [0] * self.m
            
            for i in range(self.m - 1, -1, -1):
                x[i] = extended_mat.matrix[i][extended_mat.m - 1]
                for j in range(extended_mat.m - 2, i, -1):
                    x[i] += -x[j] * extended_mat.matrix[i][j]
                x[i] = x[i] / extended_mat.matrix[i][i]
                  
        else:
            print('System has no solution')
            return None
        
        return x