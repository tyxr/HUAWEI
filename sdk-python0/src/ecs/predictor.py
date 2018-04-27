from pinghua import *
from bag import *
import copy
from getKnnData import *
from knnPre import *
def predict_vm(ecs_lines, input_lines):
    # Do your work from here#
    result = []
    if ecs_lines is None:
        print 'ecs information is none'
        return result
    if input_lines is None:
        print 'input file information is none'
        return result
    #get inputSet
    inputSet=getInput(input_lines)
    inputL=getInputL(inputSet)
    VmList=getVmsList(inputSet)
    #get trainSet
      #get pinghua data
    data,dateList=getData(ecs_lines,inputL)
      #get knn data
    dataKnn=getKnnData(ecs_lines,inputL)
    #get predictSet
    pre=getPredictSet(data,dataKnn,inputL,inputSet[-2][0],inputSet[-1][0],dateList)
    sumPre=sum(pre.values())
    #deal result
    preResult,cpuW,memW,Ccpu,Cmem,V,dataName=dealResult(inputSet,pre,inputL,VmList)
    #bag problem
      #cpuW,memW,Ccpu,Cmem,V,dataName
    bagResult=[]
    #tanxin
    cpuW=cpuW[::-1]
    memW=memW[::-1]
    V=V[::-1]
    dataName=dataName[::-1]
    bagResult,sumPack=bag(Ccpu,Cmem,copy.deepcopy(cpuW),copy.deepcopy(memW),copy.deepcopy(dataName))
    result=getFinalResult(sumPre,preResult,sumPack,bagResult)
    return result

def getFinalResult(sumPre,preResult,sumPack,bagResult):
    result=[]
    result.append(sumPre)
    for line in preResult:
        result.append(str(line[0])+' '+str(line[1]))
    result.append('')
    result.append(sumPack)
    i=1
    for line in bagResult:
        s=str(int(i))
        for key,value in line.items():
            s=s+' '+str(key)+' '+str(int(value))
        i+=1
        result.append(s)
    return result

#deal result
def dealResult(inputSet,pre,inputL,List):
    Ccpu=int(inputSet[0][0])
    Cmem=int(inputSet[0][1])
    preSet=[]
    V=[]
    for ele in inputL:
        preSet.append([ele,pre[ele]])
    cpuW,memW,dataName=getBagSet(preSet,List)
    select=inputSet[-3][0]
    if select=='CPU':
        V=copy.deepcopy(cpuW)
    elif select=='MEM':
        V=memW
    else:
        print 'Error'
    return preSet,cpuW,memW,Ccpu,Cmem,V,dataName

#get inputSet
def getInput(inputLine):
    data=[]
    for line in inputLine:
        line=line.strip()
        if line=='':
            continue
        lineData=line.split()
        data.append(lineData)
    return data

#get trainSet
def getData(fileLine,inputL):
    L=len(inputL)
    dayNum={}
    dateList=[]
    #init 
    for i in range(L):
        dayNum[inputL[i]]={}
    for line in fileLine:
        line=line.strip()
        data=line.split()
        Type=data[1]
        date=data[2].split()[0]
        #init every day
        for i in range(L):
            if date not in dayNum[inputL[i]]:
                dayNum[inputL[i]][date]=0
        if Type in dayNum:
            dayNum[Type][date]+=1
        if date not in dateList:
            dateList.append(date)
    return dayNum,dateList

#get inputL
def getInputL(inputSet):
    inputL=[]
    L=len(inputSet)
    for i in range(2,L-3):
        inputL.append(inputSet[i][0])
    return inputL

#get predictSet
def getPredictSet(data,dataKnn,inputL,beginD,endD,dateList):
    pre={}
    w1=0.5
    w2=0.6
    for ele in inputL:
        #pinghua
        vL=getParameter(ele)
        s2,s1,s0,sD=pinghua(data[ele],vL,dateList)
        L=len(s1)
        p=predict(s0[L-1],s1[L-1],s2[L-1],sD[L-1],beginD,endD,vL)
        #Knn
        X,Y=getTrainSet(dataKnn,ele,dateList)
        lastTime,lastNum=getLTN(dataKnn[ele],dateList)
        knnPre=getKnnPredictSet(lastTime,lastNum,beginD,endD,X,Y)
        #predict
        Sum=sum(p)
        iSum=int(Sum)
        if Sum-iSum>0.5:
            iSum+=1
        ip=int(w1*iSum+w2*knnPre)
        pre[ele]=ip
     #   print ip
     #   print '*'*25
    return pre
        

