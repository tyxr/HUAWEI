from dealdate import *
import copy

def getBoostTrainSet(data,L):
    LL=len(data)
    X=[]
    Y=[]
    for i in range(LL-L):
        X.append(splitDateList(data,i,L))
        Y.append(data[i+L])
    return X,Y
