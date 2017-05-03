#!/usr/bin/env python

import MF

class LinguisticVariable:
    def __init__(self, name, type="in", mfs=None, current_value = 0.0, range = None ):
        self.name = name
        if mfs is None:
            self.mfs = {}
        else:
            self.mfs = mfs
        self.type =type
        self.current_value = current_value
        
        ## Las variables de salida deben tener un Rango
        if type=="out" and range is None:
            self.range = (0,100)
        else:
            self.range = range
        
    def addMF(self,name,mf):
        if name in self.mfs:
            self.mfs[name].append(mf)
        else:
            self.mfs[name] = mf
        #self.mfs.append(MF.MF(name,mf))

class FuzzyProposition:
    def __init__(self, lv, mf, weight=1, negated = False):
        self.lv = lv
        self.mf = mf
        self.weight = weight
        self.negated = negated
    def eval(self, x=None):
        if x is None:
            return self.mf.eval(self.lv.current_value)
        else:
            return self.mf.eval(x)

class FuzzyOperator:
    def __init__(self,type,l,r):
        self.type = type
        self.l = l
        self.r =r
    def eval(self):
        if self.type == 'or':
            return max(self.l.eval(),self.r.eval())
        elif self.type == 'and':
            return min(self.l.eval(),self.r.eval())
        else:
            return "Error"
    
class FuzzyRule:
    def __init__(self, name="",antecedent = None,consequent = None):
        self.name = name
        if antecedent is None:
            self.antecedent = []
        else:
            self.antecedent = antecedent  
        
        if consequent is None:
            self.consequent = []
        else:
            self.consequent = consequent
    
    def getQualConsequence(self, out_var = 0 ):
        degreeOfSupport = self.antecedent[0].eval()
        return QualConsequence(degreeOfSupport, self.consequent[out_var].mf)
   
   
class QualConsequence:
    def __init__(self, degreeOfSupport, mf):
        self.degreeOfSupport = degreeOfSupport
        self.mf = mf
    
    def eval(self, x):
        return min(self.degreeOfSupport, self.mf.eval(x))

class FIS():
    def __init__(self, rules = None):
        if rules == None:
            self.rules = []
        else:
            self.rules = rules
    def composition(self , out_var):
		return [r.getQualConsequence(out_var) for r in self.rules]
		
    def eval(self, out_var):
        Composition = [r.getQualConsequence(out_var) for r in self.rules]

        smallest,largest = self.rules[0].consequent[out_var].lv.range
        num = 0.0
        denom = 0.0
        
        delta = (largest - smallest) / 60.0
        #delta = 0.1

        x =float( smallest)
        
        while x <= largest:
            values = [q.eval(x) for q in Composition]
            m = max(values)
            num += x * m
            denom += m;
            x+=delta
        
        if denom == 0:
           result = 0
        else:
            result = num / denom   
        return round(result,3)

if __name__ == "__main__":
    pass