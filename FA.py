import math
#import matplotlib.pyplot as plt
import numpy as np
import time

#!/usr/bin/env python
from FIS import * 



def difuso(ns):
    
    ns = ns.sort()

    EV = LinguisticVariable('ev')
    EV.addMF('N',MF.Trapezoidal(ns[0], ns[4], ns[6], ns[14]))
    EV.addMF('Z',MF.Triangular(ns[10], ns[21], ns[29]))
    EV.addMF('P',MF.Trapezoidal(ns[22], ns[30], ns[34], ns[39]))
    
    EW = LinguisticVariable('ew')
    EW.addMF('N',MF.Trapezoidal(ns[1], ns[5], ns[7], ns[15]))
    EW.addMF('Z',MF.Triangular(ns[11], ns[20], ns[28]))
    EW.addMF('P',MF.Trapezoidal(ns[23], ns[31], ns[35], ns[38]))
    
    t1 = LinguisticVariable('T1' , type = 'out', range = (-1,1))
    t1.addMF('N',MF.Triangular(ns[2], ns[8], ns[16]))
    t1.addMF('Z',MF.Triangular(ns[12], ns[19], ns[27]))
    t1.addMF('P',MF.Triangular(ns[24], ns[32], ns[37]))
    
    t2 = LinguisticVariable('T2' , type = 'out', range = (-1,1))
    t2.addMF('N',MF.Triangular(ns[3], ns[9], ns[17]))
    t2.addMF('Z',MF.Triangular(ns[13], ns[18], ns[26]))
    t2.addMF('P',MF.Triangular(ns[25], ns[33], ns[36]))
    
    # Rules
    
    r1 = FuzzyRule()
    r1.antecedent.append(FuzzyOperator('and',FuzzyProposition(EV,EV.mfs['N']),FuzzyProposition(EW,EW.mfs['N'])))
    r1.consequent.append(FuzzyProposition(T1,T1.mfs['N']),FuzzyProposition(T2,T2.mfs['N']))
    
    r2 = FuzzyRule()
    r2.antecedent.append(FuzzyOperator('and',FuzzyProposition(EV,EV.mfs['Z']),FuzzyProposition(EW,EW.mfs['Z'])))
    r2.consequent.append(FuzzyProposition(T1,T1.mfs['Z']),FuzzyProposition(T2,T2.mfs['Z']))
    
    r3 = FuzzyRule()
    r3.antecedent.append(FuzzyOperator('and',FuzzyProposition(EV,EV.mfs['P']),FuzzyProposition(EW,EW.mfs['P'])))
    r3.consequent.append(FuzzyProposition(T1,T1.mfs['P']),FuzzyProposition(T2,T2.mfs['P']))
    
    
    reglas = [r1,r2,r3]
    
    fis = FIS(reglas)
        
    def eval(ver, vi):
        EV.current_value = ErrorV
        EW.current_value = ErrorW
        return fis.eval()
    
    if __name__ == '__main__':
        print eval(12.5 , 6.31)


para = np.array([20, 500, 0.5, 0.2, 1])
convergence = []
n=para[0]  
MaxGeneration=para[1]
alpha=para[2]
betamin=para[3]
gamma=para[4]


# Numero total de Evaluaciones
NumEval=n*MaxGeneration

# Simple bounds/limits for d-dimensional problems
d = 40
Lb = np.zeros(d)
Ub = np.ones(d)

#Lb = np.ones(d)*-5.12
#Ub = np.ones(d)*5.12


zn=np.ones(n)
zn.fill(float("inf")) 

#ns(i,:)=Lb+(Ub-Lb).*rand(1,d);
ns=np.random.uniform(0,1,(n,d)) *(Ub-Lb)+Lb

Lightn=np.ones(n)
Lightn.fill(float("inf")) 
    
#[ns,Lightn]=init_ffa(n,d,Lb,Ub,u0)

def alpha_new(alpha, NGen):
    delta=1-(10**(-4)/0.9)**(1/NGen);
    alpha=(1-delta)*alpha
    return alpha

def objf(x):
    dim=len(x);
    z=np.sum(x**2-10*np.cos(2*math.pi*x))+10*dim
    return z


print("CS is optimizing  \""+objf.__name__+"\"")    
    
timerStart=time.time() 
startTime=time.strftime("%Y-%m-%d-%H-%M-%S")

# parameters [n N_iteration alpha betamin gamma]

for k in range(int(MaxGeneration)):
    #% This line of reducing alpha is optional
    alpha = alpha_new(alpha, MaxGeneration)# generando el vector de aleatoriedad alpha
    
    #% Evaluate new solutions (for all n fireflies)
    for i in range(int(n)):
        zn[i]=objf(ns[i,:])
        Lightn[i]=zn[i]
   
    # Ranking fireflies by their light intensity/objectives
    Lightn=np.sort(zn)
    Index=np.argsort(zn)
    ns=ns[Index,:]
    
    #Find the current best
    nso=ns
    Lighto=Lightn
    nbest=ns[0,:] 
    Lightbest=Lightn[0]
        
    #% For output only
    fbest=Lightbest;
  
    
   #% Move all fireflies to the better locations
#    [ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,nbest,...
#          Lightbest,alpha,betamin,gamma,Lb,Ub);   
        
          
    scale=np.ones(d)*abs(Ub-Lb)
    for i in range (int(n)):
        # The attractiveness parameter beta=exp(-gamma*r)
        for j in range(int(n)):
            r=np.sqrt(np.sum((ns[i,:]-ns[j,:])**2));
            #r=1
            # Update moves
            if Lightn[i]>Lighto[j]: # Brighter and more attractive
               beta0=1
               beta=(beta0-betamin)*math.exp(-gamma*r**2)+betamin
               tmpf=alpha*(np.random.rand(d)-0.5)*scale
               ns[i,:]=ns[i,:]*(1-beta)+nso[j,:]*beta+tmpf        
              
                
                  
    #ns=numpy.clip(ns, lb, ub)
        
    convergence.append(fbest)
        	
    IterationNumber=k
    BestQuality=fbest
        
    if (k%1==0):
           print(['At iteration '+ str(k)+ ' the best fitness is '+ str(BestQuality)])
    #    
   ####################### End main loop 
   
timerEnd=time.time()  
endTime=time.strftime("%Y-%m-%d-%H-%M-%S")
executionTime=timerEnd-timerStart
convergence=convergence

#print optimizer="FFA"
print objf.__name__      
print "Star", startTime
print "End", endTime
print "Total Time", executionTime             
                  
                       
                                
#    ############### End of Firefly Algorithm implementation ################