#include "head.h"
#include "source.cpp"
#include "mat.cpp"

/**
 * ���Զ�β��Եõ���õ�p,q 
 */
int Calculate_pq(vector<Double> data)
{
	cout<<"Calculate p q"<<endl;
    freopen("r_arma_plus.xls", "w", stdout);

    Double mean; //�������ݵľ�ֵ
    vector<Double> AutoCor;//�����ϵ��AutoCorrelation
    vector<Double> BiasCor;//ƫ���ϵ��

    AutoCor = getAutoCor(data); //�õ��������ϵ��
	BiasCor = getBiasCor(AutoCor); // �õ�ƫ���ϵ��
	//�����excel�ļ��У�ֱ����ͼ����ʾ�����ڵõ�p,q
	for(int k=0;k<AutoCor.size();k++){
        cout<<AutoCor[k]<<"\t";
    }
	cout<<endl;
    for(int k=0;k<BiasCor.size();k++){
        cout<<BiasCor[k]<<"\t";
    }

    return 0;
}



/**
 *����p�õ����ݾ�����ʽ
 */
vector<vector<Double> > y; //y = [data[p+1, ..., n]]
vector<vector<Double> > x;
/**
  [ data[p, ..., 1] ]
  [ data[p+1, ..., 2] ]  
  [ data[p+2, ..., 3] ]
  .
  .
  .
  [ data[n-1, ..., n-p] ]
 */


void formatData(vector<Double> data,int p){
	
	cout<<"formatData"<<endl;
	vector<Double> tmpy;
	
	for(int i=p;i<data.size();i++){
		tmpy.push_back(data[i]);
		vector<Double> tmp;
		for(int j=i-1;j>=i-p;j--){
			tmp.push_back(data[j]);
		}
		x.push_back(tmp);
	}
	y.push_back(tmpy);
	y = t(y);
}

/**
 *��С���˷������
 * a = inv(t(x) _*_ x) _*_ t(x) _*_ Y
 * e = sum(a) / (n-p)
 */
vector<Double> LeastSquares(vector<Double> data,int p){
	cout<<"leastSquars"<<endl;
	formatData(data,p);	
	
	vector<vector<Double> > a, tx,invx,tmp;
	tx = t(x);
	cout<<"Len tx:"<<tx.size()<<endl;
	for(int i=0;i<5;i++){
		cout<<"tx[0]["<<i<<"] ="<<tx[0][i]<<endl;
	}
	tmp = mulMat(tx, x);
	for(int i=0;i<5;i++){
		cout<<"tx[0]["<<i<<"] ="<<tmp[0][i]<<endl;
	}
	cout<<"tmp.size"<<tmp.size()<<endl;
	cout<<"tmp[0]size"<<tmp[0].size()<<endl;
	cout<<"det tmp"<<int(det(tmp))<<endl;
	invx = inv(tmp);
	a = mulMat(mulMat(invx,tx), y);
	a = t(a);
	return a[0];
}


/*
 *3�����������Ƶ� AR(PP)ģ�͵��Ƽ���в��� { ��[t] }
 */
vector<Double> getBiasSeries(vector<Double> data,vector<Double> a,int p){
	vector<Double> calPN(data.begin(),data.begin()+p);
	cout<<"getBiaseries"<<endl;
	for(int i=p;i<data.size();i++){
		Double s = 0;
		
		int t = calPN.size();
		for(int j=0;j<p;j++){
			
			s += a[j] * calPN[t-j-1];
		}
		calPN.push_back(s);
	}
	
//	cout<<calPN.size()<<endl;
	vector<Double> var;
	//����в�
	for(int i=p;i<calPN.size();i++){
		var.push_back(data[i] - calPN[i]);
	}
	
	cout<<"�в�"<<endl;
	for(int k=0;k<var.size();k++){
        cout<<var[k]<<"\t";
    }
    cout<<endl;
    return var;

}



/**
 *�Իع�ƽ����������
 *1����������ԭʼ���� x1 , ... , xn ���߽��Իع黬��ģ�� AR(pp) �����
 *	�ýϸ߽׵� AR(PP) ģ�� ȡ PP >> p,q���ܽϺõ����ԭʼ��������
 *2��ȷ��pp�󣬹���AR(pp)�в� ����1����2�� ... , ��nʹ������x[t] = ��[1]*x[t-1] + ��[2]x[t-2] + ... + ��[pp] x[t-pp] + ��[t]
 *3�����������Ƶ� AR(PP)ģ�͵��Ƽ���в��� { ��[t] }
 *4�����Ӳв��� { ��[t] } Ϊ�������� �������Իع�ģ��
 *5�����þ���ɵ� �� �� �� ����С���˹���
 *
 * ����ѡ�� p = 7�� q = 8, pp = 35��ԭ��Ӧ��ѡ��q=20�����ǵ���������ǿ��ܻᳬ��C++���Եľ��ȣ��������������ѡ��q = 8���м���
 */


/**
X = [ x[pp+q], x[pp+q-1], ..., x[pp+q-p+1] ]
	[ x[pp+q+1], x[pp+q], ..., x[pp+q-p+2] ]
	.
	.
	.
	.
	.
	[ x[n-1], x[n-2], ..., x[n-p] ]


E = [ e[p+q], e[p+q-1], ..., e[p+1] ]
	[ e[p+q+1], e[p+q], ..., e[p+2] ]
	.
	.
	.
	.
	.
	[ e[n-1],e[n-2], ..., e[n-q] ]


**/
vector<vector<Double> > X,E,xs;
void formatData(vector<Double> data,vector<Double> bias,int p, int q,int pp){
	cout<<"formatData2"<<endl;
	int n = data.size();
	for(int i=pp;i<n-q;i++){
		vector<Double> tmp;
		for(int j=0;j<p;j++){
			tmp.push_back(data[i+q-j]);
		}
		X.push_back(tmp);
	}
	
	for(int i=0;i<bias.size()-q;i++){
		vector<Double> tmp;
		for(int j=0;j<q;j++){
			tmp.push_back(bias[i+q-j]);
		}
		E.push_back(tmp);
	}
	
	vector<Double> tmpx;
	for(int i=pp+q;i<n;i++){
		tmpx.push_back(data[i]);
	}
	
	xs.push_back(tmpx);
	xs = t(xs);
}

/**
 *�õ�����,a,b
 */

vector<Double> getParm_ab(vector<Double> data,vector<Double> bias,int p, int q,int pp){
	cout<<"getParm_ab"<<endl;
	formatData(data,bias,p,q,pp);	
	
	vector<vector<Double> > ab, tx,invxe,tmp, r_XE, c_XE;
	c_XE = ConCols(t(X),t(E));
	r_XE = ConRows(X,E);
	
	invxe = inv(mulMat(c_XE,r_XE));
	tmp = mulMat(invxe,c_XE);
	
	ab = mulMat(tmp,xs);

	return t(ab)[0];	
}


/**
 *x[t] = sum[j: 0...p]{a[j]*data[t-j]} + e
 *
 *
 *ARMA ( p , q ) ģ�͵���ϼ���
 *���� ARMA ( p , q ) ģ������� AR ( p ) ģ�� MA (q ) ģ�͵ķ���������ͬ ���Ǽ�������ϲв������Ƿ�Ϊ�������� ��ͬ���� ���Ի����ϲ�
 *����������ʱʹ�ø��Բ�ͬ�����ģ�Ͷ���
 *��ʵ����
 *���� x[t] �Ƿ�Ϊ
 *ARMA ( p , q ) ʱ ֻ�����в��� { �� t } �Ƿ�������м��� ���в��еĹ���
 *ֵ { �� - t } ��������ֵ x1 , ... , xn ����ó�
 *�ж� { �� - t } �Ƿ��������
 *����������б�������еķ���
 *����Ϊ { x t } Ϊ ARMA ( p , q ) ����
 *����,������ { x t } ���� ARMA ( p , q ) ����
 *
 *1�����ģ�ͽ���H0
 *2������õ��Ĳ������õ�x t = ��[i,..,p]��[i]*x[t-i] + ��[t] - ��[k,...,q]��[k]*��[t-k]
 *3������ʽ����в�
 *4���� { ��[t] } ����Э�����
 *����ԭ��
 *�� { ��(��), k = 0,1,2, ... , n } �� Լ �� 68.3% �� �� �� �� �� �� ��� = �� 1 / n ��
 *Լ �� 95.4% �� �� �� �� �� �� �� �� = �� 2 / n ��
 *( ��[p+1] , ��[p+2] ,..., ��[n]) Ϊ������������ֵ ��ʱ����H0 ����ܾ�H0
 */

int calPQ_N(vector<Double> data,vector<Double> data_var,vector<Double> a,vector<Double> b,int p,int q){
	cout<<"calPQ_N"<<endl;
	int n = data.size();	
	//�õ��в�ڶ���x[t] - �Ʀ�[k]*x[t-k] t = [p, ... ,N]
	vector<Double> calPN(data.begin(),data.begin()+p);
		
	Double res = 0;
	for(int i=p;i<n;i++){
		Double s = 0;
		int t1 = data.size();
		for(int j=0;j<p;j++){
			s += a[j] * data[t1-j-1];
		}
		int t2 = data_var.size();
		for(int j=0;j<q;j++){
			s -= b[j] * data_var[t2-j-1];
		}		
		calPN.push_back(s);
	}	
	
	//cout<<calPN.size()<<endl;
	vector<Double> var;
	//����в�
	for(int i=p;i<calPN.size();i++){
		var.push_back(data[i] - calPN[i]);
	}
	
	vector<Double> varpq;//�õ���[p,.., n]�Ĳв����Э����ж�H0�Ƿ����
	for(int t=p;t<n;t++){
		Double tmp = 0;
		for(int j=0;j<q;j++){
			tmp += b[j]*var[t-j-1];
		}
		var.push_back(tmp - var[t-p]);
		varpq.push_back(tmp - var[t-p]);
	}
	
	vector<Double> Cor = getAutoCor(varpq);
	
	freopen("rvar_ARMA_p7_q_8.xls", "w", stdout);
	cout<<"�����ϵ��:"<<endl;
	for(int k=0;k<Cor.size();k++){
        cout<<Cor[k]<<"\t";
    }
    cout<<endl;
	
	//�����Ƿ��� 68.3% �� �� �� �� �� �� ��� = �� 1 / n ��
    //Լ �� 95.4% �� �� �� �� �� �� �� �� = �� 2 / n ��
    int k1 = 0,k2 = 0;
    Double p1 = 1.0 / Cor.size(),p2 = 2.0 / Cor.size();
    for(int k=0;k<Cor.size();k++){
    	if(Cor[k] >= -p1 && Cor[k] <= p1) k1++;	
    	if(Cor[k] >= -p2 && Cor[k] <= p2) k2++;	
    }
    cout<<"�� = �� 1 / n��"<<k1*1.0 / Cor.size()<<endl;
    cout<<"�� = �� 2 / n��"<<k2*1.0 / Cor.size()<<endl;
	return 0;	
	
}


/**
 *����ģ��Ԥ���Ժ�����ݣ�k��ʾ��k�����ݣ�����k����n,ע��ʱ�����У�������Ԥ��õ�n�����ܵõ�n+1��
 *���������k>n����Ԥ��[n,k]������λ�ã�����Ӵ�ԭ�����ϡ�
 *vector<Double> &data : ԭʼ����
 *vector<Double> &data_var: �������ݵõ��Ĳв�
 *a��b��e : ģ�Ͳ���
 *k : ��ĵ�k������
 */

Double predict(vector<Double> data,vector<Double> data_var,vector<Double> a,vector<Double> b,int p,int q,int k){
	cout<<"Double predict"<<endl;
	Double res = 0;
	for(int i=data.size();i<k;i++){
		Double s = 0;
		int t1 = data.size();
		for(int j=0;j<p;j++){
			s += a[j] * data[t1-j-1];
		}
		int t2 = data_var.size();
		for(int j=0;j<q;j++){
			s -= b[j] * data_var[t2-j-1];
		}		
		data.push_back(s);
	}
	return data[k-1];
}




Double xx[] = {871.5, 897.1, 904.3, 919.2, 935.0, 950.0, 965.0, 981.0,1028.0,1047.0,1061.0,1075.0,
              1086.0,1094.0,1110.0,1112.0,1125.0,1151.1,1159.4,1180.0,1195.6,1227.2,1243.6,
              1256.0,1128.0,1292.0,1296.0,1298.0,1302.1,1309.4,1317.0,1332.6,1235.2,1363.6,
              1385.1,1423.2,1456.4,1472.7,1488.0,1491.0,1501.0,1511.0,1520.0,1531.9,1538.6,1540.3
};

Double xd[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,18,19,20};

int main()
{
	cout<<"main()"<<endl;
	vector<Double> data;
	int p=7,q=8,pp=20;//һ��ע��p��q��ȡֵ��ͨ�����ݼ���󣬹��Ƴ����ġ�
	//��������
	for(int i=0;i<46;i++){
		data.push_back(xx[i]);
	}
	//����p,q,ͨ��ͼ����ʾ��ѡ��p = 7�� q = 20, pp = 38
    //	Calculate_pq(data);
	
	vector<Double> ta = LeastSquares(data,pp);
	cout<<"����ARģ�͵õ��Ĳ���ta����:  "<<ta.size()<<endl;
	for(int i=0;i<ta.size();i++){
		cout<<"ta["<<i<<"] = "<<ta[i]<<endl;
	}
	
	//�в�	
	vector<Double> bias = getBiasSeries(data,ta,pp);
	/**
	for(int i=0;i<bias.size();i++){
		cout<<"var["<<i<<"] = "<<bias[i]<<endl;
	}
	**/
	vector<Double> ab = getParm_ab(data,bias,p,q,pp);

	vector<Double> a(ab.begin(),ab.begin()+p);	
	vector<Double> b(ab.begin()+p,ab.begin()+p+q);
	cout<<"����a����:  "<<a.size()<<endl;
	for(int i=0;i<a.size();i++){
		cout<<"a["<<i<<"] = "<<a[i]<<endl;
	}
	cout<<"����b����:  "<<b.size()<<endl;
	for(int i=0;i<b.size();i++){
		cout<<"b["<<i<<"] = "<<b[i]<<endl;
	}
	
	//calPQ_N(data,bias,a,b,p,q);
	
	cout<< predict(data,bias,a,b,p,q,47)<<endl;; 
	
	
	return 0;
}



















