import copy
def det(x):
    xx = copy.deepcopy(x)
    ddet = 1
    iteration = 0
    for i in range(len(xx)):
        if xx[i][i]==0:
            
            for j in range(i+1,len(xx)):
                if xx[j][i]!=0:
                    
                    xx[i],xx[j]=xx[j],xx[i]
                    
                    iteration+=1
        for k in range(i+1,len(xx)):
            if xx[i][i]==0:
                return 0
            yin = -1.0*xx[k][i]/xx[i][i]
            for u in range(len(xx[0])):
                xx[k][u] = xx[k][u]+xx[i][u]*yin
    for i in range(len(xx)):
        ddet = ddet*xx[i][i]
    if iteration%2==1:
        ddet = -ddet
    return ddet
def t(x):
    tx = [[] for i in x[0]]
    for i in x:
        for j in range(len(i)):
            tx[j].append(i[j])
    return tx
def mulMat(tx,x):#tx [m,n] x[n,p] txx[m,p]

    res = [[0] * len(x[0]) for i in range(len(tx))]
    for i in range(len(tx)):
        for j in range(len(x[0])):
            for k in range(len(x)):
                res[i][j] += tx[i][k] * x[k][j]
    return res
def det2(m):
    if len(m) <= 0:
        return None
    elif len(m) == 1:
        return m[0][0]
    else:
        s = 0
        for i in range(len(m)):
            n = [[row[a] for a in range(len(m)) if a != i] for row in m[1:]]  
            s += m[0][i] * det(n) * (-1) ** (i % 2)
        return s
def det3(m):  
    if len(m) <= 0:  
        return None  
    if len(m) == 1:  
        return m[0][0]  
    else:  
        s = 0  
        for i in range(len(m)):  
            n = [[row[a] for a in range(len(m)) if a != i] for row in m[1:]]  
            if i % 2 == 0:  
                s += m[0][i] * det(n)  
            else:  
                s -= m[0][i] * det(n)  
    return s
def det4(mat):
    
    
    n = len(mat[0])
    res = 1
    
    for col in range(n):
        row = col
        res *= mat[row][col]
        
        while mat[row][col] == 0 and row < n - 1:
            row += 1
        
        for i in range(row + 1, n):
            if mat[i][col] == 0:
                pass
            else:
                k = - mat[i][col] / mat[row][col]
                for j in range(col ,n):
                    mat[i][j] += mat[row][j] * k
    return res
def delMat(x,r,c):
    Ax = []
    for i in range(len(x)):
        temp = []
        for j in range(len(x[0])):
            if i!=r and j !=c :
                temp.append(x[i][j])
        if i!=r:
            Ax.append(temp)
    return Ax
def A(x):
    tmp = []
    res = []
    for i in range(len(x)):
        tp = [0 for _ in range(len(x[0]))]
        res.append(tp)
    for i in range(len(x)):
        for j in range(len(x[0])):
            tmp = x
            tmp = delMat(tmp,i,j)

            res[i][j]=(1 if (i+j)%2==0 else -1)*det(tmp)
    return t(res)


def inv(x):
    
    
    dets = float(det(x))

    res = A(x)
    
    
    for i in range(len(res)):
        for j in range(len(res[0])):
            res[i][j]/=dets
    return res
def ConRows(x,y):
    for i in range(len(y)):
        for j in range(len(y[0])):
            x[i].append(y[i][j])
    return x
def ConCols(x,y):
    for i in range(len(y)):
        row = []
        for j in range(len(y[0])):
            row.append(y[i][j])
        x.append(row)
    return x
def test_Mat():
    
    x =[[2,1,-1],[2,1,0],[1,-1,1]]

    Ax = inv(x)

    print Ax



if __name__ == '__main__':
    test_Mat()
