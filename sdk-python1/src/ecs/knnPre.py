import time
import datetime
import copy

#get test set
def getKnnPredictSet(lastDate,lastNum,beforeT,beginT,endT,X,Y):
    #analyze time
    last=time.strptime(lastDate,'%Y-%m-%d')
    before=time.strptime(beforeT,'%Y-%m-%d')
    begin=time.strptime(beginT[0],'%Y-%m-%d')
    end=time.strptime(endT[0],'%Y-%m-%d')
    lastD=datetime.datetime(last[0],last[1],last[2])
    beforeDate=datetime.datetime(before[0],before[1],before[2])
    beginDate=datetime.datetime(begin[0],begin[1],begin[2])
    endDate=datetime.datetime(end[0],end[1],end[2])
    if int(beginT[1].split(':')[0])>12:
        beginDate=beginDate+datetime.timedelta(days=1)
    if int(endT[1].split(':')[0])>12:
        endDate=endDate+datetime.timedelta(days=1)
    beforeDate=beforeDate+datetime.timedelta(days=1)
    #count data
    Sum=0
    #copy the data
    copyX=copy.deepcopy(X)
    copyY=copy.deepcopy(Y)
    while beforeDate<endDate:
        #get test data
        lastT=int(str(beforeDate-lastD).split()[0])
        oneList=[beforeDate.weekday(),lastT,lastNum]
        #predict
        onePre=KnnPredict(copyX,copyY,oneList,3)
        if beforeDate>=beginDate:
            Sum+=onePre
        if onePre>=1:
            lastD=beforeDate
            lastNum=onePre
        #move time
        beforeDate=beforeDate+datetime.timedelta(days=1)
    return int(Sum)

def KnnPredict(X,Y,testData,n):
    #make the data nomal
    xf1=[x[0] for x in X]
    xf1Max=max(xf1)
    xf1Max=max(xf1Max,testData[0])
    xf1Min=min(xf1)
    xf1Min=min(xf1Min,testData[0])
    x1=[float(x-xf1Min)/float(xf1Max-xf1Min) for x in xf1]
    testData[0]=float(testData[0]-xf1Min)/float(xf1Max-xf1Min)

    xf2=[x[1] for x in X]
    xf2Max=max(xf2)
    xf2Max=max(xf2Max,testData[1])
    xf2Min=min(xf2)
    xf2Min=min(xf2Min,testData[1])
    x2=[float(x-xf2Min)/float(xf2Max-xf2Min) for x in xf2]
    testData[1]=float(testData[1]-xf2Min)/float(xf2Max-xf2Min)

    xf3=[x[2] for x in X]
    xf3Max=max(xf3)
    xf3Max=max(xf3Max,testData[2])
    xf3Min=min(xf3)
    xf3Min=min(xf3Min,testData[2])
    x3=[float(x-xf3Min)/float(xf3Max-xf3Min) for x in xf3]
    testData[2]=float(testData[2]-xf3Min)/float(xf3Max-xf3Min)

    loc=0
    numDic={}
    #begin KNN
    for a1,a2,a3 in zip(x1,x2,x3):
        numDic[loc]=(a1-testData[0])**2+(a2-testData[1])**2+(a3-testData[2])**2
        loc+=1
    numList=sorted(numDic.items(),key=lambda e:e[1])
    #get the best n
    num=0
    for i in range(0,n):
        num+=Y[numList[i][0]]
    onePre=float(num)/float(n)
    if onePre>0:
        onePre=int(onePre)+1
    return onePre
