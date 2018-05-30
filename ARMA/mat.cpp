#include "head.h"

/**
 * a = inv(t(x) _*_ x) _*_ t(x) _*_ Y
 * t(x) : ��x��ת��
 * r(x) : ��x�����
 * -*- : ����˷�
 * inv(x): �������
 */

/**
 *����ת��
 */
vector<vector<Double> > t(vector<vector<Double> > x){
    //x��װ�þ���
    vector<vector<Double> > tx;
    //tx��ʼ������ֱ�ӷ����±�,����ԭ�����ת�õ���ʽ
    for(int i=0;i<x[0].size();i++){
        vector<Double> tmp(x.size(),0);
        tx.push_back(tmp);
    }

    for(int i=0;i<x.size();i++){
        for(int j=0;j<x[0].size();j++){
            tx[j][i] = x[i][j];
        }
    }
    return tx;
}

/**
 *����˷�
 */
vector<vector<Double> > mulMat(vector<vector<Double> > tx, vector<vector<Double> > x){
    vector<vector<Double> > res;
    //��ʼ���������ĸ�ʽrow(tx) X col(x)
    for(int i=0;i<tx.size();i++){
        vector<Double> tmp(x[0].size(),0);
        res.push_back(tmp);
    }

//    cout<<res.size()<<" "<<res[0].size()<<endl;
    for(int i=0;i<tx.size();i++){
        for(int j=0;j<x[0].size();j++){
            for(int k=0;k<x.size();k++){
               res[i][j] += tx[i][k] * x[k][j];
            }
        }
    }
    return res;
}

/**
 *���������ʽ�����б仯Ϊ�����Ǿ���
 */

Double det(vector<vector<Double> > x){
	//ֻ��һ��Ԫ��
	//if(x.size() == 1 && x[0].size() == 1) return x[0][0];
	
	Double det = 1;
	//��������ָ�������У��������б任������Ϊ�н�����
    int iter = 0;  //��¼�б任�Ĵ�����������
    for(int i=0;i<x.size();i++){
        if(x[i][i] == 0){
        	for(int j=i+1;j<x.size();j++){
	        	if(x[j][i] != 0){
	                swap(x[i],x[j]);//��������
    	            iter ++;
    	        }
    	    }
        }
        for(int k=i+1;k<x.size();k++){
            Double yin = -1 * x[k][i] / x[i][i] ;
            for(int u=0; u<x[0].size(); u++){
                x[k][u] = x[k][u] + x[i][u] * yin;
            }
        }
    }
	
	/**
   	cout<<"�����Ǿ���"<<endl;
   	for(int i=0;i<x.size();i++){
        for(int j=0;j<x[0].size();j++){
            cout<<x[i][j]<<" ";
        }
        cout<<endl;
    }**/
	for(int i=0;i<x.size();i++){//��Խ��ߵĻ� �� ����ʽ��ֵ
      det = det * x[i][i];
	}
	//�б任ż���η��Ų���
	if(iter%2 == 1)  det= -det;

	return det;
}

/**
 *ɾ������ĵ�r�У���c��
 */
vector<vector<Double> > delMat(vector<vector<Double> > x,int r,int c){
	vector<vector<Double> > Ax;
	for(int i=0;i<x.size();i++){
		vector<Double> tmp;
		for(int j=0;j<x[0].size();j++){
			if(i != r && j != c) tmp.push_back(x[i][j]);			
		}
		if(i != r) Ax.push_back(tmp);
	}	
	return Ax;
}


/**
 *�����İ������
 */
vector<vector<Double> > A(vector<vector<Double> > x){
	vector<vector<Double> > tmp(x),res;
	
	//tx��ʼ������ֱ�ӷ����±�,����ԭ�����ת�õ���ʽ
    for(int i=0;i<x.size();i++){
        vector<Double> tp(x[0].size(),0);
        res.push_back(tp);
    }
    
	for(int i=0;i<x.size();i++){
		for(int j=0;j<x[0].size();j++){
			tmp = x;
			tmp = delMat(tmp,i,j);
			res[i][j] = ((i+j)%2==0?1:-1) * det(tmp);
			
		}
	}
	return t(res);
}


/**
 *�������
 */
vector<vector<Double> > inv(vector<vector<Double> > x){
	vector<vector<Double> > res = A(x);
	Double dets = det(x);
	for(int i=0;i<res.size();i++){
		for(int j=0;j<res[0].size();j++){
			res[i][j] /= dets;
		}
	}
	return res;
}


/**
 *�ϲ���������ͬ�ľ���
 */ 
vector<vector<Double> > ConRows(vector<vector<Double> > x, vector<vector<Double> > y){
	//����ͬ�������
	for(int i=0;i<y.size();i++){
		for(int j=0;j<y[0].size();j++){
			x[i].push_back(y[i][j]);
		}
	}
	return x;
}

/**
 *�ϲ���������ͬ�ľ���
 */ 
vector<vector<Double> > ConCols(vector<vector<Double> > x, vector<vector<Double> > y){
	//����ͬ�������
	for(int i=0;i<y.size();i++){
		vector<Double> row;
		for(int j=0;j<y[0].size();j++){
			row.push_back(y[i][j]);
		}
		x.push_back(row);
	}
	return x;
}






/**
 *���Ծ�������ɹ�
 */
void test_Mat(){
    vector<vector<Double> > data,tdata,res,Ax;
    //data = getdata();
	Double x[] = {2,1,-1,2,1,0,1,-1,1};
	
    for(int i=0;i<3;i++){
        vector<Double> tmp;
        data.push_back(tmp);
        for(int j=0;j<3;j++){
            data[i].push_back(x[i*3+j]);
        }
    }
	
	/**
    tdata = t(data);

    for(int i=0;i<tdata.size();i++){
        for(int j=0;j<tdata[0].size();j++){
            cout<<tdata[i][j]<<" ";
        }
        cout<<endl;
    }

    res = mulMat(tdata,data);
    for(int i=0;i<res.size();i++){
        for(int j=0;j<res[0].size();j++){
            cout<<res[i][j]<<" ";
        }
        cout<<endl;
    }
	
	cout<<det(data)<<endl;
	*/
	
	//Ax = inv(data);
	cout<<det(data)<<endl;
	cout<<"�����:"<<endl;
	for(int i=0;i<Ax.size();i++){
        for(int j=0;j<Ax[0].size();j++){
            cout<<Ax[i][j]<<" ";
        }
        cout<<endl;
    }
}

/**
int main()
{

    vector<Double> data;


    test_Mat();
    return 0;
}

**/



