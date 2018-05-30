from mat import *



xx = [871.5, 897.1, 904.3, 919.2, 935.0, 950.0, 965.0, 981.0,1028.0,1047.0,1061.0,1075.0,   \
              1086.0,1094.0,1110.0,1112.0,1125.0,1151.1,1159.4,1180.0,1195.6,1227.2,1243.6,  \
              1256.0,1128.0,1292.0,1296.0,1298.0,1302.1,1309.4,1317.0,1332.6,1235.2,1363.6,  \
              1385.1,1423.2,1456.4,1472.7,1488.0,1491.0,1501.0,1511.0,1520.0,1531.9,1538.6,1540.3]
x = []
y = []

X = []
E = []
xs = []




def formatData_0(data,pp):
    print 'formatData_0'
    global y
    tmpy = []
    
    for i in range(pp,len(data)):
        tmpy.append(data[i])
        tmp = []
        j_rev = range(i-pp,i)
        j_rev.reverse()
        for j in j_rev:
            tmp.append(data[j])
        x.append(tmp)
    y.append(tmpy)
    y = t(y)
    print 'formatData_0 end'



def LeastSquares(data,pp):
    print 'LeastSquares'
    tmpy = []
    formatData_0(data,pp)
    tx = t(x)
    for i in range(5):
        print tx[0][i]
    print 'tx end'
    temp = mulMat(tx,x)
    print 'temp_end'
    print len(temp)
    print len(temp[0])
    invx = inv(temp)
    print 'inv end'
    a = mulMat(mulMat(invx,tx),y)
    print 'mulmat end'
    a = t(a)
    print 'LeastSquares end'
    return a[0]

def getBiasSeries(data,a,p):
    calpn =data[:p]
    for i in range(p,len(data)):
        s = 0
        t = len(calpn)
        for j in range(p):
            s+= a[j]*calpn[t-j-1]
        calpn.append(s)
        t = len(calpn)
    var = []
    for i in range(p,len(calpn)):
        var.append(data[i]-calpn[i])
    return var
def formatData_1(data,bias,p,q,pp):
    global xs
    n = len(data)
    for i in range(pp,n-q):
        tmp = []
        for j in range(0,p):
            tmp.append(data[i+q-j])
        X.append(tmp)
    for i in range(len(bias)-q):
        tmp = []
        for j in range(q):
            tmp.append(bias[i+q-j])
        E.append(tmp)
    tmpx = []
    for i in range(pp+q,n):
        tmpx.append(data[i])
    xs.append(tmpx)
    xs = t(xs)
def get_parm_ab(data,bias,p,q,pp):
    formatData_1(data,bias,p,q,pp)
    c_XE = ConCols(t(X),t(E))
    r_XE = ConRows(X,E)
    invxe = inv(mulMat(c_XE,r_XE))
    tmp = mulMat(invxe,c_XE)
    ab = mulMat(tmp,xs)
    return t(ab)[0]
def predict(data,data_var,a,b,p,q,k):
    res = 0
    for i in range(len(data),k):
        s = 0
        t1 = len(data)
        for j in range(p):
            s+=a[j]*data[t1-j-1]
        t2 = len(data_var)
        for j in range(q):
            s-=b[j]*data_var[t2-j-1]
        data.append(s)
    return data[k-1]




'''def calPQ_N(data,data_var,a,b,p,q):
    n = len(data)
    calpn =data[:p]
    res = 0
    for i in range(p,n):
        s = 0
        t1 = n
        for j in range(p):
            s+=a[j]*data[t1-j-1]
        t2 = len(data_var)
        for j in range(q):
            s -= b[j]*data_var[t2-j-1]
        calpn.append(s)
    var = []
    for i in range(p,len(calpn)):
        var.append(data[i]-calpn[i])
    verpq = []
    for i in range(p,n):
        tmp = 0.0
        for j in range(q):
            tmp+=b[j]*var[t-j-1]
        var.append(tmp-var[t-p])
        varpq.append(tmp-var[t-p])
        cor = getAutoCor(varpq)'''
        
def main():
    print 'main'
    data = xx
    p = 7
    q = 8
    pp = 20
    ta = LeastSquares(data,pp)

    bias =getBiasSeries(data,ta,pp)

    
    ab = get_parm_ab(data,bias,p,q,pp)
    a = ab[:p]
    b = ab[p:]

    print(predict(data,bias,a,b,p,q,47))


if __name__=="__main__":
    main()
