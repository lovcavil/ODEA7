#pragma once
#include <vector>
typedef std::vector< double > vecd;
/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;



inline double pow2(double x)
{
	return x * x;
}
class bearCalc {
public:
	vecd Fbxi_;
	vecd Fbyi_;
	vecd Fbzi_;
	vecd Mxi_;
	vecd Myi_;
	vecd theta; //= zeros(Nb, 1);% t的函数 长Nb的向量
	vecd deltasum;// = zeros(Nb, 1); % 5自由度和t的函数 长Nb的向量
	vecd Q;// = zeros(Nb, 1); % 5自由度和t的函数 长Nb的向量
	//vecd fai_r_theta(Nb + 1);
	vecd deltazj;
	vecd deltayj;
	vecd A0pr;
	vecd alp0pr;
	vecd rmovement;
public:
	bearCalc(int Nb) {
		Fbxi_ = vecd(Nb + 1);
		Fbyi_ = vecd(Nb + 1);
		Fbzi_ = vecd(Nb + 1);
		Mxi_ = vecd(Nb + 1);
		Myi_ = vecd(Nb + 1);
		theta = vecd(Nb + 1); //= zeros(Nb, 1);% t的函数 长Nb的向量
		deltasum = vecd(Nb + 1);// = zeros(Nb, 1); % 5自由度和t的函数 长Nb的向量
		Q = vecd(Nb + 1);// = zeros(Nb, 1); % 5自由度和t的函数 长Nb的向量
		//vecd fai_r_theta(Nb + 1);
		deltazj = vecd(Nb + 1);
		deltayj = vecd(Nb + 1);
		A0pr = vecd(Nb + 1);
		alp0pr = vecd(Nb + 1);
		rmovement = vecd(Nb + 1);
	}

	vecd bearCalc_20190807_withDiag(double t, double xbi, double ybi, double zbi, double faixi, double faiyi, double Rj, double A0, double alp0, double Nb, double wcage, double Kc, double Pd, double dm, double Db, int install_disp, int bearing_num);
};

class bear {
public:
	double t;
	double xbi;
	double ybi;
	double zbi;
	double faixi;
	double faiyi;
	double Rj;
	double A0;
	double alp0;
	double Nb=8;
	double wcage;
	double Kc;
	double Pd;
	double dm;
	double Db; 
	int install_disp;
	int bearing_num;
	double Fbx;
	double Fby;
	double Fbz;
	double Mx;
	double My;
	static bearCalc m_bearCalc;
public:
	bear(double t_, double xbi_, double ybi_, double zbi_, double faixi_, double faiyi_,
		double Rj_, double A0_, double alp0_, double Nb_, double wcage_, double Kc_, double Pd_, double dm_, double Db_,
		int install_disp_, int bearing_num_)
		: t(t_), xbi(xbi_), ybi(ybi_), zbi(zbi_), faixi(faixi_), faiyi(faiyi_),
		Rj(Rj_), A0(A0_), alp0(alp0_), Nb(Nb_), wcage(wcage_), Kc(Kc_), Pd(Pd_), dm(dm_), Db(Db_),
		install_disp(install_disp_), bearing_num(bearing_num_) {};
	void calc();
};


