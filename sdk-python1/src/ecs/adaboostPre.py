from adaboost import *
import time
import datetime

#product new data
def createNewData(X,Y,pre):
    oneList=[]
    LL=len(X[0])
    for i in range(1,LL):
        oneList.append(X[-1][i])
    oneList.append(Y[-1])
    X.append(oneList)
    Y.append(pre)

#get new trainset
def getNewTrainSet(X,Y):
    oneList=[]
    LL=len(X[0])
    for i in range(1,LL):
        oneList.append(X[-1][i])
    oneList.append(Y[-1])
    return oneList

#adaboost predict
def predict(lastD,beginD,endD,X,Y,num,stepLen,border,errThreshold,h):
    #dealtime
    lastT=time.strptime(lastD,'%Y-%m-%d')
    beginT=time.strptime(beginD[0],'%Y-%m-%d')
    endT=time.strptime(endD[0],'%Y-%m-%d')
    dateL=datetime.datetime(lastT[0],lastT[1],lastT[2])
    dateB=datetime.datetime(beginT[0],beginT[1],beginT[2])
    dateE=datetime.datetime(endT[0],endT[1],endT[2])
    if int(beginD[1].split(':')[0])>12:
        dateB=dateB+datetime.timedelta(days=1)
    if int(endD[1].split(':')[0])>12:
        dateE=dateE+datetime.timedelta(days=1)
    #cal TL,DL
    TL=int(str(dateB-dateL).split()[0])-1
    DL=int(str(dateE-dateB).split()[0])
    dateL=dateL+datetime.timedelta(days=1)
    #train
    classTree,treeWeight=adaboost(X,Y,num,stepLen,border,errThreshold,h)
    #init
    copyX=copy.deepcopy(X)
    copyY=copy.deepcopy(Y)
    p=[]
    dateL=dateB
    #begin
    while dateL<dateE:
        #get train set
        oneList=getNewTrainSet(copyX,copyY)
        #predict
        pre=adaPre(oneList,classTree,treeWeight)
        if dateL>=dateB:
            p.append(pre)
        #create new data
        createNewData(copyX,copyY,pre)
        #move time
        dateL=dateL+datetime.timedelta(days=1)
    return p,TL,DL
            
