#include "head.h"
/**
 *���ļ����������㷨���õ��ĺ���
 */




/**
 *��Э���� AutoCov[k] = E((x[i] - u)(x[i-k] - u))
 *�����ϵ�� AutoCov[k] = AutoCov[k] / AutoCov[0]
 */
 
vector<Double> getAutoCov(vector<Double> data){
	//���������ϵ������
    int n = data.size();

    Double mean = 0; //���ݵľ�ֵ
    for(int i=0;i<n;i++){
        mean += data[i];
    }
    mean /= n;
	//cout<<"mean::"<<mean<<endl;
    //��ÿ�����ݶ���ȥ��ֵ�õ��µ�����
    vector<Double> prodata;
    
    for(int i=0;i<n;i++){
        prodata.push_back(data[i] - mean);
		//cout<<"prodata[i] "<<prodata[i]<<endl;
    }

    vector<Double> AutoCov(n,0);//��Э����AutoCovariance
    for(int k=0;k<n;k++){
        for(int i=0;i<n-k;i++){
            AutoCov[k] += prodata[i] * prodata[i+k];
        }
        AutoCov[k] /= n - k;
    }

	return AutoCov;
}



vector<Double> getAutoCor(vector<Double> data){
    
	vector<Double> AutoCor,AutoCov;//�����ϵ��AutoCorrelation��ע���±��0��ʼ
    AutoCov = getAutoCov(data);
    
    /**
    cout<<"AutoCor:��Э���begin"<<endl;
    for(int k=0;k<AutoCov.size();k++){
    	cout<<AutoCov[k]<<" ";
    }
    cout<<"AutoCor:��Э���end"<<endl;
    **/
    
    for(int k=0;k<data.size()-1;k++){
        AutoCor.push_back(AutoCov[k+1] / AutoCov[0]);
    }
	
	/**
	cout<<"AutoCor:��Э����ϵ����begin"<<endl;
    for(int k=0;k<AutoCor.size();k++){
    	cout<<AutoCor[k]<<" ";
    }
    cout<<"AutoCor:��Э����ϵ����end"<<endl;
    **/
    return AutoCor;
}


/**
 *�õ�ƫ���ϵ��BiasCor[k,k]
 *BiasCor[0,0] = AutoCor[0]
 *BiasCor[k,k] = (AutoCor[k-1] - sum[j:0...k-1]{AutoCor[k-j]*BiasCor[j,k-1]}) / (1 - sum[j:0...k-1]AutoCor[j]*BiasCor[j,k-1])
 *BiasCor[j,k] = BiasCor[j,k-1] - BiasCor[k,k]*BiasCor[k-j,k-1] j = 0...k
 *
 */
vector<Double> getBiasCor(vector<Double> AutoCor){
    //����BiasCor[i,j],Ϊ��ֱ�ӷ����±꣬���ȳ�ʼ��
    vector< vector<Double> > BiasCor;
    for(int i=0;i<AutoCor.size();i++){
        vector<Double> tmp(AutoCor.size(),0);
        BiasCor.push_back(tmp);
    }

    BiasCor[0][0] = AutoCor[0];

    for(int k=1;k<AutoCor.size();k++){
        BiasCor[k][k] = AutoCor[k];
        Double t1,t2;
        for(int j=0;j<=k-1;j++){
            t1 = AutoCor[k-j] * BiasCor[j][k-1];
            t2 = AutoCor[j] * BiasCor[j][k-1];

            BiasCor[j][k] = BiasCor[j][k-1] - BiasCor[k][k] * BiasCor[k-j][k-1];

        }
        BiasCor[k][k] = (BiasCor[k][k] - t1) / t2;
        for(int j=0;j<=k-1;j++){
            BiasCor[k][j] = BiasCor[j][k] = BiasCor[j][k-1] - BiasCor[k][k] * BiasCor[k-j][k-1];
        }
    }
    vector<Double> res;
    for(int k=0;k<AutoCor.size();k++){
        res.push_back(BiasCor[k][k]);
    }

    return res;

}
