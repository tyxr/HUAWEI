import dealdate
import math
import copy

#init data weight
def initWeight(num,para):
    weight=[]
    for i in range(num):
        weight.append(para)
    return weight

#mat x * mat y
def mul(x,y):
    L=len(x)
    Sum=0
    for i in range(L):
        Sum+=x[i]*y[i]
    return Sum

def initAllStump(xL,stepLen,border):
    stump=initWeight(xL,0)
    allStump=[]
    loc=0
    while loc<xL:
        if stump[loc]>(border-stepLen):
            loc+=1
            if loc>=xL:
                break
            stump[loc]+=stepLen
        elif stump[loc]<border:
            stump[loc]+=stepLen
        allStump.append(copy.deepcopy(stump))
    return allStump


#init stump
def createStump(X,Y,stepLen,border):
    allStump=[]
    Ares=[]
    L=len(Y)
    xL=len(X[0])
    #predict
    allStump=initAllStump(xL,stepLen,border)      
    sL=len(allStump)
    for i in range(sL):
        Are=[]
        for j in range(L):
            pre=mul(allStump[i],X[j])
            if Y[j]==0:
                Are.append(abs(float(pre-Y[j])/3.0))
            else:
                Are.append(abs(float(pre-Y[j])/float(Y[j])))   
        Ares.append(Are)
    return Ares,allStump


def selectBestStump(Ares,dataWeight,errThreshold,paras):
    #init
    L=len(Ares)
    L1=len(Ares[0])
    minErr=100000
    minLoc=-1
    minPara=[]
    for i in range(L):
        #get each sum weight
        err=0
        r=0
        for j in range(L1):
            if Ares[i][j]>errThreshold:
                err+=dataWeight[j]
            else:
                r+=1
        #updata
        if err<minErr:
            minErr=err
            minLoc=i
            minPara=paras[i]
    return minPara,minErr,minLoc

def normalize(x):
    Sum=float(sum(x))
    L=len(x)
    for i in range(L):
        x[i]=x[i]/Sum

def updateDataWeight(dataWeight,Are,errThreshold,beta):
    D=copy.deepcopy(dataWeight)
    L=len(D)
    for i in range(L):
        if Are[i]<=errThreshold:
            D[i]=D[i]*beta
    normalize(D)
    return D

#use adaboost to train data
  #X:trainset
  #Y:label
  #num:the num of stump
def adaboost(X,Y,num,stepLen,border,errThreshold,h):
    m=len(X)
    #init para
    dataWeight=initWeight(m,1.0/float(m))
    classTree=[]
    treeWeight=[]
    #get all trainSet
    Ares,paras=createStump(X,Y,stepLen,border)
    for i in range(num):
        #are:residual
        para,err,loc=selectBestStump(Ares,dataWeight,errThreshold,paras)
        beta=err**h
        #update dateWeight
        dataWeight=updateDataWeight(dataWeight,Ares[loc],errThreshold,beta)  #
        #update para
        treeWeight.append(math.log(1.0/float(beta)))
        classTree.append(para)
    normalize(treeWeight)
    return classTree,treeWeight


#predict
def adaPre(x,classTree,treeWeight):
    L=len(classTree)
    result=0
    for i in range(L):
        pre=mul(classTree[i],x)
        result+=pre*treeWeight[i]
    if result<0:
        result=0
    result=int(result+0.5)
    return result
        
