import copy
def bag(itemName,itemWeight,bagName,bagWeight,allCpu,allMem):
    #init bag
    resultSet={}
    resultLoc={}
    bagTypeLen=len(bagName)
    for i in range(bagTypeLen):
        resultSet[bagName[i]]=[]
        resultLoc[bagName[i]]=-1
    #copy the parameter
    copyItemName=copy.deepcopy(itemName)
    copyItemWeight=copy.deepcopy(itemWeight)
    copyBagName=copy.deepcopy(bagName)
    copyBagWeight=copy.deepcopy(bagWeight)
    #open the first bag
    allPara=getAllPara(itemWeight) 
    bagType,bagPara,bCpu,bMem=openBag(bagWeight,bagName,allPara)
    resultSet[bagType].append({})
    resultLoc[bagType]=0
    bagLoc=0
    end=0
    #cal item para
    itemPara=getItemPara(itemWeight,bagPara) #
    #get the item len
    itemL=len(itemName)
    while itemL>0 and end==0:
        flag=0
        i=0
        while i<itemL:
            itemLoc=itemPara[i][0]
            if bCpu>=copyItemWeight[itemLoc][0] and bMem>=copyItemWeight[itemLoc][1]:
                #put item in bag
                Name=copyItemName[itemLoc]
                if Name in resultSet[bagType][bagLoc]:
                    resultSet[bagType][bagLoc][Name]+=1
                else:
                    resultSet[bagType][bagLoc][Name]=1
                #updata bag weight
                bCpu-=copyItemWeight[itemLoc][0]
                bMem-=copyItemWeight[itemLoc][1]
                #update item
                del copyItemWeight[itemLoc]
                del copyItemName[itemLoc]
                #update flag
                flag=1
                #update item para
                if bMem==0 or bCpu==0:
                    flag=0
                    itemL-=1
                    break
                bagPara=float(bCpu)/float(bMem)
                itemPara=getItemPara(copyItemWeight,bagPara)
                #updata itemL
                itemL-=1
            else:
                i+=1
        #open new bag
        if flag==0:
            bagType=canPutOne(bagWeight,bagName,copyItemWeight,copyItemName,allCpu,allMem,resultSet)
            #the item can put in one bag
            if bagType!='no':
                resultSet[bagType].append({})
                resultLoc[bagType]+=1
                bagLoc=resultLoc[bagType]
                putAllInBag(resultSet,bagType,bagLoc,copyItemWeight,copyItemName)
                end=1
                break
            #the item can not put in one bag
            allPara=getAllPara(copyItemWeight)
            bagType,bagPara,bCpu,bMem=openBag(bagWeight,bagName,allPara)
            resultSet[bagType].append({})
            resultLoc[bagType]+=1
            bagLoc=resultLoc[bagType]
            itemPara=getItemPara(copyItemWeight,bagPara)
    return resultSet,resultLoc

#put all the item into bag
def putAllInBag(resultSet,bagType,bagLoc,copyItemWeight,copyItemName):
    #get len
    Len=len(copyItemName)
    bag=resultSet[bagType][bagLoc]
    for i in range(Len):
        if copyItemName[i] not in bag:
            bag[copyItemName[i]]=1
        else:
            bag[copyItemName[i]]+=1


#test if the item can put in one bag
def canPutOne(bagWeight,bagName,copyItemWeight,copyItemName,allCpu,allMem,resultSet):
    #get len
    bagL=len(bagName)
    itemL=len(copyItemName)
    #cal all item weigth
    allItemCpu=0
    allItemMem=0
    #can put bag
    canPutList=[]
    #init para
    bagType='no'
    for i in range(itemL):
        allItemCpu+=copyItemWeight[i][0]
        allItemMem+=copyItemWeight[i][1]
    #test
    for i in range(bagL):
        if bagWeight[i][0]>=allItemCpu and bagWeight[i][1]>=allItemMem:
            canPutList.append(i)
    if canPutList!=[]:
        bagType=chooseBestToPut(allCpu,allMem,resultSet,canPutList,bagName,bagWeight)
    return bagType

def chooseBestToPut(allCpu,allMem,resultSet,canPutList,bagName,bagWeight):
    L=len(canPutList)
    bL=len(bagName)
    bagCpu=0
    bagMem=0
    bagType='no'
    Max=-1
    if L==1:
        return bagName[canPutList[0]]
    else:
        #cal all bag use
        for j in range(bL):
            Name=bagName[j]
            bagNum=len(resultSet[Name])
            bagCpu+=bagNum*bagWeight[j][0]
            bagMem+=bagNum*bagWeight[j][1]
        for i in range(L):
            loc=canPutList[i]
            value=(allCpu/float(bagCpu+bagWeight[loc][0])+allMem/float(bagMem+bagWeight[loc][1]))*0.5
            if value>Max:
                Max=value
                bagType=bagName[loc]
        return bagType
            

    
#get item para
def getItemPara(itemWeight,bagPara):
    #cal parameter
    L=len(itemWeight)
    calPara={}
    for i in range(L):
        para=float(itemWeight[i][0])/float(itemWeight[i][1])
        calPara[i]=abs(bagPara-para)
    cLPara=sorted(calPara.items(),key=lambda item:item[1],reverse=False)
    return cLPara     

#get all item para
def getAllPara(itemWeight):
    cpu=0
    mem=0
    for ele in itemWeight:
        cpu+=ele[0]
        mem+=ele[1]
    return float(cpu)/float(mem)

#open new bag
def openBag(bagWeight,bagName,allPara):
    minDis=100000
    minLoc=-1
    minPara=-1
    bagLen=len(bagName)
    for i in range(bagLen):
        para=float(bagWeight[i][0])/float(bagWeight[i][1])
        Dis=abs(para-allPara)
        if Dis<minDis:
            minDis=Dis
            minLoc=i
            minPara=para
    bagType=bagName[minLoc]
    bagPara=minPara
    bCpu=bagWeight[minLoc][0]
    bMem=bagWeight[minLoc][1]
    return bagType,bagPara,bCpu,bMem


#ok
def getVmsList(data):
    List={}
    beginLoc=int(data[0][0])+2
    L=len(data)
    for i in range(beginLoc,L-2):
        Type=data[i][0]
        List[Type]={}
        List[Type]['cpu']=int(data[i][1])
        List[Type]['mem']=int(data[i][2])
    return List
        

#ok
def changeToData(name,List):
    return List[name]['cpu'],float(List[name]['mem']/1024)

#ok
def getBagSet(data,List):
    itemW=[]
    dataName=[]
    for line in data:
        num=line[1]
        if num!=0:
            for i in range(num):
                dataName.append(line[0])
                cpu,mem=changeToData(line[0],List)
                itemW.append([cpu,mem])
    return itemW,dataName
