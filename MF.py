#!/usr/bin/env python

import math

class MF:
    def __init__ (self,term, mf):
        self.term = term    
        self.mf = mf
    def eval(self, x):
        return self.mf.eval(x)


class Sigmoid():
    def __init__ (self,a,c):
        self.a = float(a)
        self.c = float(c)

    def eval(self,x):
        try:
            exp_result = math.exp( -self.a * (x-self.c) )
        except:
            exp_result =  math.exp(709)


        result =  1.0 / exp_result
        return result
                

class Triangular():
    def __init__ (self,a,b,c):
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)
        
    def eval(self,x):
        if x <= self.a:
            result = 0.0
        elif self.a <= x <= self.b:
            result = (x-self.a) / (self.b-self.a)
        elif self.b <= x  <= self.c:
            result = (self.c-x) / (self.c- self.b)
        elif self.c <= x:
            result = 0.0
        return result       


class Bell():
    def __init__ (self,a,b,c):
        if a == 0:
            return "Error: invalid parameter"
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)
        
    def eval(self,x):
        return 1/(1+ math.pow( math.fabs((x-self.c)/self.a),2*self.b))
 
           
class Trapezoidal():
    def __init__ (self,a,b,c,d):

        self.a = float(a)
        self.b = float(b)
        self.c = float(c)
        self.d = float(d)
   
    def eval(self,x):
        if self.a == self.b == x:
            return 1.0
        if self.c == self.d == x :
            return 1.0

        if x <= self.a:
           result = 0.0
        if self.a <= x <= self.b:
           result = (x-self.a) / (self.b-self.a)
        if self.b <= x <= self.c:
            result = 1.0
        if self.c <= x  <= self.d:
            result = (self.d-x) / (self.d- self.c)
        if self.d <= x:
            result = 0.0
        return result

class Gaussian():
    def __init__ (self,sigma,c ):
        self.c = float(c)
        self.sigma = float(sigma)
   
    def eval(self,x):
        return math.exp(-1*(x-self.c)**2 / (2*self.sigma**2) )


 
if __name__ == "__main__":
 
 params = (20.0,25.0,30)
 tri = Triangular(*params)
 print "Triangular"
 print tri.eval(0.5)
    
 params = (0.0, 1.0, 3.0, 4.0)
 tri = Trapezoidal(*params)
 print "Trapezoidal"
 print tri.eval(1.0)

 params = (0.0, 0.0, 4.0, 4.0)
 tri = Trapezoidal(*params)
 print "Trapezoidal"
 print tri.eval(4.0)



 print "Gaussian"
 gauss = Gaussian(1.5,10.0)
 print gauss.eval(3.0)

