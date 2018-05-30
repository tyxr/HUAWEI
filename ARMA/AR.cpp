/**
 *������֪����ֵ x 1 , x 2 , L , x n �� AR ( p ) ģ���������� ��Ϊ�Իع�ģ������Իع�
 *ģ��������ݰ���
 *1 AR ( p ) ģ�ͽ��� p �Ĺ���
 *2 AR ( p ) ģ���в��� �� 1 , �� 2 , L , �� p �� �� 2 �Ĺ���
 *3 ��ģ������ϼ���
 */
#include "head.h"
#include "mat.cpp"
#include "source.cpp"


int Calculate_p(vector<Double> data)
{

    freopen("r_AR.xls", "w", stdout);

    Double mean; //�������ݵľ�ֵ
    vector<Double> AutoCor;//�����ϵ��AutoCorrelation
    vector<Double> BiasCor;//ƫ���ϵ��

    AutoCor = getAutoCor(data); //�õ��������ϵ��
    BiasCor = getBiasCor(AutoCor); // �õ�ƫ���ϵ��

    for(int k=0;k<BiasCor.size();k++){
        cout<<BiasCor[k]<<"\t";
    }

    return 0;
}



/**
 *ͳ�ƺ�������Կ�����k=16֮��|BiasCor[k]| < 5,���ѡ��p = 16֮�����������
 *֮��ʼ�����������������
 *ʹ����С���˷��������a �� e��
 * a = inv(t(x) _*_ x) _*_ t(x) _*_ Y
 * e = sum(a) / (n-p)
 * t(x) : ��x��ת��
 * r(x) : ��x�����
 * -*- : ����˷�
 * inv(x): �������
 */


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
	
	/**
	cout<<"X:"<<endl;
	for(int i=0;i<x.size();i++){
		for(int j=0;j<x[0].size();j++){
			cout<<x[i][j]<<" ";
		}
		cout<<endl;
	}
	
	cout<<"y:"<<endl;
	for(int i=0;i<y.size();i++){
		for(int j=0;j<y[0].size();j++){
			cout<<y[i][j]<<" ";
		}
		cout<<endl;
	}**/
	
}

/**
 *��С���˷������
 * a = inv(t(x) _*_ x) _*_ t(x) _*_ Y
 * e = sum(a) / (n-p)
 */
vector<Double> LeastSquares(vector<Double> data,int p){
	formatData(data,p);	
	
	vector<vector<Double> > a, tx,invx,tmp;
	tx = t(x);
	invx = inv(mulMat(tx, x));
	//cout<<invx.size()<<endl;
	
	/**
	cout<<"invx:"<<endl;
	cout<<invx.size()<<" "<<invx[0].size()<<endl;
	for(int i=0;i<invx.size();i++){
		for(int j=0;j<invx[0].size();j++){
			cout<<invx[i][j]<<" ";
		}
		cout<<endl;
	}
	
	cout<<"tx:"<<endl;
	cout<<tx.size()<<" "<<tx[0].size()<<endl;
	for(int i=0;i<tx.size();i++){
		for(int j=0;j<tx[0].size();j++){
			cout<<tx[i][j]<<" ";
		}
		cout<<endl;
	}
	**/
	a = mulMat(mulMat(invx,tx), y);
	a = t(a);
	return a[0];
}

/**
 *�õ�e
 */
Double getBias(vector<Double> data,vector<Double> a,int n,int p){
	Double sum = 0;
	vector<Double> calPN(data.begin(),data.begin()+p);

	
	for(int i=p;i<data.size();i++){
		Double s = 0;
		int t = calPN.size();
		for(int j=0;j<p;j++){
			s += a[j] * calPN[t-j-1];
		}
		calPN.push_back(s);
	}
	
	//cout<<calPN.size()<<endl;
	//����в�
	for(int i=p;i<calPN.size();i++){
		sum += (data[i] - calPN[i]);
	}
	
	return sum / (n-p);
}


/**
 *����ģ��
 *1���������H
 *2�����ݵõ��Ĳ�������� data[p+1, ..., n]
 *3������в
 *4����в����Э����ϵ�����ж�H�Ƿ����
 *
 *x[t] = sum[j: 0...p]{a[j]*data[t-j]} + e
 *�� { ��(��), k = 0,1,2, ... , n } �� Լ �� 68.3% �� �� �� �� �� �� ��� = �� 1 / n ��
 *Լ �� 95.4% �� �� �� �� �� �� �� �� = �� 2 / n ��
 *( ��[p+1] , ��[p+2] ,..., ��[n]) Ϊ������������ֵ ��ʱ����H0 ����ܾ�H0
 
 */

int calP_N(vector<Double> data,vector<Double> a,int p){

	int n = data.size();
	vector<Double> calPN(data.begin(),data.begin()+p);
	
	for(int i=p;i<data.size();i++){
		Double s = 0;
		int t = calPN.size();
		for(int j=0;j<p;j++){
			s += a[j] * calPN[t-j-1];
		}
		calPN.push_back(s);
	}
	
	//cout<<calPN.size()<<endl;
	vector<Double> var;
	//����в�
	for(int i=p;i<calPN.size();i++){
		var.push_back(data[i] - calPN[i]);
	}
	
    vector<Double> Avar;//�����ϵ��AutoCorrelation
    vector<Double> Bvar;//ƫ���ϵ��

    Avar = getAutoCor(var); //�õ��������ϵ��
    Bvar = getBiasCor(Avar); // �õ�ƫ���ϵ��
	/**
	for(int k=0;k<Avar.size();k++){
        cout<<Avar[k]<<"\t";
  
	*/
	freopen("rvar16.xls", "w", stdout);

	cout<<"�����ϵ��:"<<endl;
	for(int k=0;k<Avar.size();k++){
        cout<<Avar[k]<<"\t";
    }
    cout<<endl;
    
    //�����Ƿ��� 68.3% �� �� �� �� �� �� ��� = �� 1 / n ��
    //Լ �� 95.4% �� �� �� �� �� �� �� �� = �� 2 / n ��
    int k1 = 0,k2 = 0;
    Double p1 = 1.0 / Avar.size(),p2 = 2.0 / Avar.size();
    for(int k=0;k<Avar.size();k++){
    	if(Avar[k] >= -p1 && Avar[k] <= p1) k1++;	
    	if(Avar[k] >= -p2 && Avar[k] <= p2) k2++;	
    }	
    cout<<"�� = �� 1 / n��"<<k1*1.0 / Avar.size()<<endl;
    cout<<"�� = �� 2 / n��"<<k2*1.0 / Avar.size()<<endl;
    
    /**
     *���֮���ֶ��ڸ��������У�ֻ��
     *�� 16.7% �� �� �� �� �� �� ��� = �� 1 / n ��
     *�� 22.2% �� �� �� �� �� �� �� �� = �� 2 / n ��
     *��H0�ļ��費�ܱ����ܣ������ڻỻ�����ݽ��в��ԡ�
     */
    
    return 0;
}

/**
 *����ģ��Ԥ���Ժ�����ݣ�k��ʾ��k�����ݣ�����k����n,ע��ʱ�����У�������Ԥ��õ�n�����ܵõ�n+1��
 *���������k>n����Ԥ��[n,k]������λ�ã�����Ӵ�ԭ�����ϡ�
 */

Double predict(vector<Double> &data,vector<Double> a,int k,int p){
	Double res;
	for(int i=data.size();i<k;i++){
		Double s = 0;
		int t = data.size();
		for(int j=0;j<p;j++){
			s += a[j] * data[t-j-1];
		}
		data.push_back(s);
	}
	return data[k-1];
}


//����1987-2014�˿�: 35
Double xx[] = {871.5,897.1,904.3,919.2,935.0,950.0,965.0,981.0,1028.0,1047.0,1061.0,1075.0,
              1086.0,1094.0,1102.0,1112.0,1125.0,1251.1,1259.4,1240.0,1245.6,1257.2,1363.6,
            1385.1,1423.2,1456.4,1492.7,1538.0,1601.0,1676.0,1771.0,1860.0,1961.9,2018.6,2069.3
};

int main()
{
	//��������
	int p = 16;//ԭ����p����ѡ��[16-N],��������p����18ʱ������ᳬ�����ȣ��ֿ��ǵ�ģ�;����ܼ򵥣������������٣�ѡ��p=[15,16,17]
				//����ֱ���excel��¼��ͼ������
    vector<Double> data,a;
	for(int i=0;i<35;i++){
		data.push_back(xx[i]);
	}	
	
	//����p
	//Calculate_p(data);
	
	a = LeastSquares(data,p);	
	
	cout<<"����a����:  "<<a.size()<<endl;
	for(int i=0;i<a.size();i++){
		cout<<"a["<<i<<"] = "<<a[i]<<endl;
	}

	cout<<endl;
	
	//cout<<e<<endl;
	
	//�����㷨������,��������p=[16,17]ͨ����ͼ���۲��㷨������
	calP_N(data,a,p);
	
	//Ԥ���37�����ݵ�ֵ
	Double x = predict(data,a,37,p);	
	cout<<x<<endl;
    return 0;
}









