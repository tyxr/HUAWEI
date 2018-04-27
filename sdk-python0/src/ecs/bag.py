def bag(c1,c2,w1,w2,dataName):
    #the bag set
    resultSet=[{}]
    #init
      #init len
    L=len(w1)
      #init loc:bag num
    bagloc=0
      #init c
    cCpu=c1
    cMem=c2
      #init para
    bagPara=float(cMem)/float(cCpu)
      #init itemPara
    itemPara=getItemPara(bagPara,w1,w2)
    while L>0:
        flag=0
        i=0
        while i<L:
            itemloc=itemPara[i][0]
            Name=dataName[itemloc]
            if w1[itemloc]<=cCpu and w2[itemloc]<=cMem:
                #put thing in bag
                if Name in resultSet[bagloc]:
                    resultSet[bagloc][Name]+=1
                else:
                    resultSet[bagloc][Name]=1
                #updata cpu and mem
                cCpu-=w1[itemloc]
                cMem-=w2[itemloc]
                #updata w
                del w1[itemloc]
                del w2[itemloc]
                del dataName[itemloc]
                #renew itemPara
                if cCpu>0 and cMem>0:
                    bagPara=float(cMem)/float(cCpu)
                    itemPara=getItemPara(bagPara,w1,w2)
                    L-=1
                    flag=1
                    i=0
                else:
                    flag=0
                    L-=1
                    break
            else:
                i+=1
        #use new bag
        if flag==0:
            resultSet.append({})
            bagloc+=1
            cCpu=c1
            cMem=c2
            bagPara=float(cMem)/float(cCpu)
            itemPara=getItemPara(bagPara,w1,w2)
    #avoid error
    bagLen=bagloc+1
    for i in range(bagLen):
        if resultSet[i]=={}:
            del resultSet[i]
            bagLen-=1
    return resultSet,int(bagLen)


def getItemPara(bagPara,w1,w2):
    #cal parameter
    L=len(w1)
    calPara={}
    for i in range(L):
        para=float(w2[i])/float(w1[i])
        calPara[i]=abs(bagPara-para)
    cLPara=sorted(calPara.items(),key=lambda item:item[1],reverse=False)
    return cLPara


def getVmsList(data):
    List={}
    L=len(data)
    for i in range(2,L-3):
        Type=data[i][0]
        List[Type]={}
        List[Type]['cpu']=int(data[i][1])
        List[Type]['mem']=int(data[i][2])
    return List
        


def changeToData(name,List):
    return List[name]['cpu'],(List[name]['mem']/1024)


def getBagSet(data,List):
    cpuW=[]
    memW=[]
    dataName=[]
    for line in data:
        num=line[1]
        if num!=0:
            for i in range(num):
                dataName.append(line[0])
                cpu,mem=changeToData(line[0],List)
                cpuW.append(cpu)
                memW.append(mem)
    return cpuW,memW,dataName
