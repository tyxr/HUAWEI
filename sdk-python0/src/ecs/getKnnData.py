import time
import datetime


#get last time
def findLastTN(beforeSet,Type,date):
    lastTime='0'
    lastNum=-1
    t=0
    for ele in beforeSet:
        beforeDate=beforeSet[ele]
        if beforeDate>0 and ele>lastTime:
            lastTime=ele
            lastNum=beforeDate
    if lastTime!='0': 
        d=time.strptime(date,'%Y-%m-%d')
        l=time.strptime(lastTime,'%Y-%m-%d')
        t=int(str(datetime.datetime(d[0],d[1],d[2])-datetime.datetime(l[0],l[1],l[2])).\
           split()[0])
    return t,lastNum

#get each data
def analyze(date,dataSet,beforeSet,Type):
    #get message from date
    #week,(is jiejiari),last happen time,last happen num
    d=time.strptime(date,'%Y-%m-%d')
    week=d.tm_wday                #week
    lastTime,lastNum=findLastTN(beforeSet,Type,date)  #last happen analyze
    dataSet.append([week,lastTime,lastNum])


#deal data
def getKnnData(File,inputL):
    L=len(inputL)
    dayNum={}
    #init vm
    for i in range(L):
        dayNum[inputL[i]]={'X':[],'Y':{}}
    for line in File:
        #first get data
        data=line.split()
        Type=data[1]        #Type is the vm's type
        date=data[2].split(' ')[0]  #date is the date
        #count each day num
        for i in range(L):
            if date not in dayNum[inputL[i]]['Y']:
                #get feather
                analyze(date,dayNum[inputL[i]]['X'],dayNum[inputL[i]]['Y'],inputL[i])
                dayNum[inputL[i]]['Y'][date]=0
        if Type in dayNum:
            dayNum[Type]['Y'][date]+=1
    return dayNum

#get trainset
def getTrainSet(dayNum,ele,dateList):
    X=dayNum[ele]['X']
    dataY=dayNum[ele]['Y']
    Y=[]
    for e in dateList:
        Y.append(dataY[e])
    return X,Y

#get the happen last time and last num
def getLTN(date,dateList):
    lastNum=0
    lastTime=None
    for ele in dateList:
        Num=date['Y'][ele]
        if Num!=0:
            lastNum=Num
            lastTime=ele
    return lastTime,lastNum

#get p prediect
def getPPre(y,Py):
    L=len(y)
    p=[]
    for i in range(L):
        p.append(y[i]-Py[i])
    return p
