/*
 Copyright 2010-2012 Karsten Ahnert
 Copyright 2011-2013 Mario Mulansky
 Copyright 2013 Pascal Germroth
 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

//#define _MATH_DEFINES_DEFINED
#include<math.h>
#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <armadillo>
#include <corecrt_math_defines.h>
#include <thread>
#include "rapidjson/document.h"     // rapidjson's DOM-style API
#include "rapidjson/prettywriter.h" // for stringify JSON
#include <cstdio>
#include <fstream>
#include <chrono>   
#include "Source.h"
using namespace rapidjson;
using namespace std::chrono;

typedef std::vector< double > vecd;

/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

//[ rhs_class
/* The rhs of x' = f(x) defined as a class */
class harm_osc {

    double m_gam;

public:
    harm_osc( double gam ) : m_gam(gam) { }

    void operator() ( const state_type &x , state_type &dxdt , const double /* t */ )
    {
        dxdt[0] = x[1];
        dxdt[1] = -x[0] - m_gam*x[1];
    }
};
//]

std::string zerosSTR(int fileNum)
{
	std::stringstream ret;
	ret <<std::setfill('0')<<std::setw(5)<<fileNum;
	return ret.str().c_str();
	//std::to_string(fileNum)
}

void WriteFile(std::vector< state_type > states, std::vector< double > times,long long step,long long fileNum)
{
	std::ofstream outfile("./a/allMyNumbers" +zerosSTR(fileNum)+ ".csv", std::ios::out | std::ios::binary);
	std::vector<double>::iterator it2 = times.begin();
	outfile << "fs,s,time,zero,thetap,thetappr,xb1,xb1pr,yb1,yb1pr,zb1,zb1pr,faix1,faix1pr,faiy1,faiy1pr,x1,x1pr,y1,y1pr,z1,z1pr,theta1x,theta1xpr,theta1y,theta1ypr,theta1z,theta1zpr,xb2,xb2pr,yb2,yb2pr,zb2,zb2pr,faix2,faix2pr,faiy2,faiy2pr,xb3,xb3pr,yb3,yb3pr,zb3,zb3pr,faix3,faix3pr,faiy3,faiy3pr,x2,x2pr,y2,y2pr,z2,z2pr,theta2x,theta2xpr,theta2y,theta2ypr,theta2z,theta2zpr,xb4,xb4pr,yb4,yb4pr,zb4,zb4pr,faix4,faix4pr,faiy4,faiy4pr,thetag,thetagpr"<<std::endl;
	for (std::vector<state_type>::iterator it1 = states.begin();
		it1 != states.end() && it2 != times.end();
		++it1, ++it2)
	{
		outfile << step << "," << fileNum << "," << *it2 << ",";
		std::ostream_iterator<double> o_iter(outfile, ",");
		std::copy(it1->begin(), it1->end(), o_iter);
		outfile << std::endl;
	}
	states.clear();
	times.clear();
	return;
}

vecd bearCalc(double t, double  xbi, double  ybi, double zbi, double faixi, double faiyi,
	double Rj, double A0, double alp0, double Nb, double wcage, double Kc, double Pd,
	int install_disp);

ODE_MODEL::ODE_MODEL()
{

}

void ODE_MODEL::readini()
{
	std::ifstream in("cfg.ini");
	std::istreambuf_iterator<char> beg(in), end;
	std::string str(beg, end);
	Document document;  // Default template parameter uses UTF8 and MemoryPoolAllocator.

	char p[16000];
	strcpy_s(p, str.c_str());
	char buffer[sizeof(p)];
	memcpy(buffer, p, sizeof(p));
	if (document.ParseInsitu(buffer).HasParseError())
		;

	b = document["b"].GetDouble();
	bc = document["bc"].GetDouble();
	cb1 = document["cb1"].GetDouble();
	cb2 = document["cb2"].GetDouble();
	cb3 = document["cb3"].GetDouble();
	cb4 = document["cb4"].GetDouble();
	cm = document["cm"].GetDouble();
	cn1 = document["cn1"].GetDouble();
	cn2 = document["cn2"].GetDouble();
	csx1 = document["csx1"].GetDouble();
	csx2 = document["csx2"].GetDouble();
	csx3 = document["csx3"].GetDouble();
	csx4 = document["csx4"].GetDouble();
	csy1 = document["csy1"].GetDouble();
	csy2 = document["csy2"].GetDouble();
	csy3 = document["csy3"].GetDouble();
	csy4 = document["csy4"].GetDouble();
	csz1 = document["csz1"].GetDouble();
	csz2 = document["csz2"].GetDouble();
	csz3 = document["csz3"].GetDouble();
	csz4 = document["csz4"].GetDouble();
	ct1 = document["ct1"].GetDouble();
	ct2 = document["ct2"].GetDouble();
	Db = document["Db"].GetDouble();
	em = document["em"].GetDouble();
	er = document["er"].GetDouble();
	etas = document["etas"].GetDouble();
	etat = document["etat"].GetDouble();
	fi = document["fi"].GetDouble();
	fo = document["fo"].GetDouble();
	g = document["g"].GetDouble();
	Ix1 = document["Ix1"].GetDouble();
	Ix2 = document["Ix2"].GetDouble();
	Ix3 = document["Ix3"].GetDouble();
	Ix4 = document["Ix4"].GetDouble();
	Iy1 = document["Iy1"].GetDouble();
	Iy2 = document["Iy2"].GetDouble();
	Iy3 = document["Iy3"].GetDouble();
	Iy4 = document["Iy4"].GetDouble();
	J1 = document["J1"].GetDouble();
	J2 = document["J2"].GetDouble();
	Jg = document["Jg"].GetDouble();
	Jp = document["Jp"].GetDouble();
	km = document["km"].GetDouble();
	kn1 = document["kn1"].GetDouble();
	kn2 = document["kn2"].GetDouble();
	ksx1 = document["ksx1"].GetDouble();
	ksx2 = document["ksx2"].GetDouble();
	ksx3 = document["ksx3"].GetDouble();
	ksx4 = document["ksx4"].GetDouble();
	ksy1 = document["ksy1"].GetDouble();
	ksy2 = document["ksy2"].GetDouble();
	ksy3 = document["ksy3"].GetDouble();
	ksy4 = document["ksy4"].GetDouble();
	ksz1 = document["ksz1"].GetDouble();
	ksz2 = document["ksz2"].GetDouble();
	ksz3 = document["ksz3"].GetDouble();
	ksz4 = document["ksz4"].GetDouble();
	kt1 = document["kt1"].GetDouble();
	kt2 = document["kt2"].GetDouble();
	m1 = document["m1"].GetDouble();
	m2 = document["m2"].GetDouble();
	mb1 = document["mb1"].GetDouble();
	mb2 = document["mb2"].GetDouble();
	mb3 = document["mb3"].GetDouble();
	mb4 = document["mb4"].GetDouble();
	Pd = document["Pd"].GetDouble();
	r = document["r"].GetDouble();
	R = document["R"].GetDouble();
	R1 = document["R1"].GetDouble();
	R1 = document["R1"].GetDouble();
	R2 = document["R2"].GetDouble();
	R2 = document["R2"].GetDouble();
	ru1 = document["ru1"].GetDouble();
	ru2 = document["ru2"].GetDouble();
	Tdm = document["Tdm"].GetDouble();
	Tdr = document["Tdr"].GetDouble();
	Tgm = document["Tgm"].GetDouble();
	Tgr = document["Tgr"].GetDouble();
	xig = document["xig"].GetDouble();
	Z1 = document["Z1"].GetDouble();
	Z2 = document["Z2"].GetDouble();
	an = document["an"].GetDouble();
	al = document["al"].GetDouble();
	beta = document["beta"].GetDouble();
	betab = document["betab"].GetDouble();

	//calc some input
	beta = beta / 180. * M_PI;
	betab = betab / 180. * M_PI;
	an = an / 180. * M_PI;
	at = atan(tan(an) / cos(beta));
	al = al / 180. * M_PI;

	// read option
	flag_zero_thetap = document["flag_zero_thetap"].GetBool();
	flag_zero_thetag = document["flag_zero_thetag"].GetBool();
	flag_zero_fai = document["flag_zero_fai"].GetBool();
	flag_limit_fai = document["flag_limit_fai"].GetBool();
	flag_gear_off = document["flag_gear_off"].GetBool();
	flag_bearing_off = document["flag_bearing_off"].GetBool();
	flag_zero_cage_movement = document["flag_zero_cage_movement"].GetBool();
}

double ODE_MODEL::limit_fai(double fai)
{
	if (fai > (M_PI_2 / 18.))
	{
		fai = (M_PI_2 / 18.);
	}
	else if (fai < -(M_PI_2 / 18.))
	{
		fai = -(M_PI_2 / 18.);
	}
	return fai;
}

void ODE_MODEL::operator()(const state_type& X, state_type& dxdt, const double t)
{
	double thetap = X[1];
	double thetappr = X[2];
	double xb1 = X[3];
	double xb1pr = X[4];
	double yb1 = X[5];
	double yb1pr = X[6];
	double zb1 = X[7];
	double zb1pr = X[8];
	double faix1 = X[9];
	double faix1pr = X[10];
	double faiy1 = X[11];
	double faiy1pr = X[12];
	double x1 = X[13];
	double x1pr = X[14];
	double y1 = X[15];
	double y1pr = X[16];
	double z1 = X[17];
	double z1pr = X[18];
	double theta1x = X[19];
	double theta1xpr = X[20];
	double theta1y = X[21];
	double theta1ypr = X[22];
	double theta1z = X[23];
	double theta1zpr = X[24];
	double xb2 = X[25];
	double xb2pr = X[26];
	double yb2 = X[27];
	double yb2pr = X[28];
	double zb2 = X[29];
	double zb2pr = X[30];
	double faix2 = X[31];
	double faix2pr = X[32];
	double faiy2 = X[33];
	double faiy2pr = X[34];
	double xb3 = X[35];
	double xb3pr = X[36];
	double yb3 = X[37];
	double yb3pr = X[38];
	double zb3 = X[39];
	double zb3pr = X[40];
	double faix3 = X[41];
	double faix3pr = X[42];
	double faiy3 = X[43];
	double faiy3pr = X[44];
	double x2 = X[45];
	double x2pr = X[46];
	double y2 = X[47];
	double y2pr = X[48];
	double z2 = X[49];
	double z2pr = X[50];
	double theta2x = X[51];
	double theta2xpr = X[52];
	double theta2y = X[53];
	double theta2ypr = X[54];
	double theta2z = X[55];
	double theta2zpr = X[56];
	double xb4 = X[57];
	double xb4pr = X[58];
	double yb4 = X[59];
	double yb4pr = X[60];
	double zb4 = X[61];
	double zb4pr = X[62];
	double faix4 = X[63];
	double faix4pr = X[64];
	double faiy4 = X[65];
	double faiy4pr = X[66];
	double thetag = X[67];
	double thetagpr = X[68];
	if (flag_zero_fai) {
		faix1 = 0;
		faix1pr = 0;
		faiy1 = 0;
		faiy1pr = 0;

		faix2 = 0;
		faix2pr = 0;
		faiy2 = 0;
		faiy2pr = 0;

		faix3 = 0;
		faix3pr = 0;
		faiy3 = 0;
		faiy3pr = 0;

		faix4 = 0;
		faix4pr = 0;
		faiy4 = 0;
		faiy4pr = 0;
	}
	if (flag_zero_thetag) {
		thetag = 0;
	}
	if (flag_zero_thetap) {
		thetap = 0;
	}
	if (flag_limit_fai) {
		faix1 = limit_fai(faix1);
		faix2 = limit_fai(faix2);
		faix3 = limit_fai(faix4);
		faix4 = limit_fai(faix4);
		faiy1 = limit_fai(faiy1);
		faiy2 = limit_fai(faiy2);
		faiy3 = limit_fai(faiy3);
		faiy4 = limit_fai(faiy4);
	}
	if (flag_bearing_off) {
		xb1 = xb2 = xb3 = xb4 = 0;
		yb1 = yb2 = yb3 = yb4 = 0;
		zb1 = zb2 = zb3 = zb4 = 0;
		faix1 = faix2 = faix3 = faix4 = 0;
		faiy1 = faiy2 = faiy3 = faiy4 = 0;
		xb1pr = xb2pr = xb3pr = xb4pr = 0;
		yb1pr = yb2pr = yb3pr = yb4pr = 0;
		zb1pr = zb2pr = zb3pr = zb4pr = 0;
		faix1pr = faix2pr = faix3pr = faix4pr = 0;
		faiy1pr = faiy2pr = faiy3pr = faiy4pr = 0;
	}

	dm = (r + R) / 1000;// % 轴承的节圆直径 r + R unit : m
	A0 = (fo + fi - 1) * Db;
	alp0 = acos(1 - Pd / 2 / A0);
	Rj = 0.5 * dm + (fi * Db - 0.5 * Db) * cos(alp0); //% now suppose as unit : m
	double wi = 500. * 2. * M_PI / 60.;//% 轴承内圈角速度
	//	%% to check 含复杂滚动轴承建模的航空发动机整机振动耦合动力学模型_陈果 R r TO BE ro ri
	double 	wcage = wi * r / (r + R);
	if (flag_zero_cage_movement) {
		wcage = 0;
	}


	//%% n
	double omega1 = 2. * M_PI * 500. / 60.;
	double omega2 = 2. * M_PI * 500. / 60.;
	double faip = omega1 * t + thetap;
	double fai1 = omega1 * t + theta1z;
	double fai2 = omega2 * t + theta2z;
	double faig = omega2 * t + thetag;

	double cosfai1pr = -sin(fai1) * (omega1 + theta1zpr);
	double cosfai2pr = -sin(fai2) * (omega2 + theta2zpr);
	double sinfai1pr = cos(fai1) * (omega1 + theta1zpr);
	double sinfai2pr = cos(fai2) * (omega2 + theta2zpr);

	//?/%% torque
	double Tp = Tdm + Tdr * sin(omega1 * t + faip);
	double Tg = Tgm + Tgr * sin(omega2 * t + faig);
	//% Tp = 400;
	//% Tg = 400;

	//%% bearing

	double gama = Db * cos(alp0) / dm;
	double zrui = 1. / Db * (4. - 1. / fi + 2. * gama / (1. - gama)) / 1000.; //% to unit : mm - 1
	double zruo = 1. / Db * (4. - 1. / fo - 2. * gama / (1. + gama)) / 1000.;
	double Frui = (1. / fi + 2. * gama / (1. - gama)) / (4. - 1. / fi + 2. * gama / (1. - gama));
	double Fruo = (1. / fo - 2. * gama / (1. + gama)) / (4. - 1. / fo - 2. * gama / (1. + gama));
	double deltai_asterisk = -327.6145 + 1883.338 * Frui - 3798.1121 * pow(Frui, 2) + 3269.6154 * pow(Frui, 3) - 1026.96 * pow(Frui, 4);
	double deltao_asterisk = -327.6145 + 1883.338 * Fruo - 3798.1121 * pow(Fruo, 2) + 3269.6154 * pow(Fruo, 3) - 1026.96 * pow(Fruo, 4);

	double ki = 2.15e5 * pow(zrui, -0.5) * pow(deltai_asterisk, -1.5);
	double ko = 2.15e5 * pow(zruo, -0.5) * pow(deltao_asterisk, -1.5);
	Kc = pow((1. / (pow(ki, -1. / 1.5) + pow(ko, -1. / 1.5))), 1.5);

	//%% bear1
	vecd BearCalc;
	bear bear1 = bear(t, xb1, yb1, zb1, faix1, faiy1,// ... % var
		Rj, A0, alp0, Nb, wcage, Kc, Pd, dm, Db, -1, 1);
	bear bear2 = bear(t, xb2, yb2, zb2, faix2, faiy2,// ... % var
		Rj, A0, alp0, Nb, wcage, Kc, Pd, dm, Db, 1, 2);
	bear bear3 = bear(t, xb3, yb3, zb3, faix3, faiy3,// ... % var
		Rj, A0, alp0, Nb, wcage, Kc, Pd, dm, Db, -1, 3);
	bear bear4 = bear(t, xb4, yb4, zb4, faix4, faiy4,// ... % var
		Rj, A0, alp0, Nb, wcage, Kc, Pd, dm, Db, 1, 4);

	double Fbx1 = 0;
	double Fby1 = 0;
	double Fbz1 = 0;
	double Mx1 = 0;
	double My1 = 0;
	double Fbx2 = 0;
	double Fby2 = 0;
	double Fbz2 = 0;
	double Mx2 = 0;
	double My2 = 0;
	double Fbx3 = 0;
	double Fby3 = 0;
	double Fbz3 = 0;
	double Mx3 = 0;
	double My3 = 0;
	double Fbx4 = 0;
	double Fby4 = 0;
	double Fbz4 = 0;
	double Mx4 = 0;
	double My4 = 0;
	if (flag_bearing_off) {
		;
	}
	else {
		bear1.calc();
		Fbx1 = bear1.Fbx;
		Fby1 = bear1.Fby;
		Fbz1 = bear1.Fbz;
		Mx1 = bear1.Mx;
		My1 = bear1.My;

		bear2.calc();
		Fbx2 = bear2.Fbx;
		Fby2 = bear2.Fby;
		Fbz2 = bear2.Fbz;
		Mx2 = bear2.Mx;
		My2 = bear2.My;

		bear3.calc();
		Fbx3 = bear3.Fbx;
		Fby3 = bear3.Fby;
		Fbz3 = bear3.Fbz;
		Mx3 = bear3.Mx;
		My3 = bear3.My;

		bear4.calc();
		Fbx4 = bear4.Fbx;
		Fby4 = bear4.Fby;
		Fbz4 = bear4.Fbz;
		Mx4 = bear4.Mx;
		My4 = bear4.My;
	}



	double thetappr2 = (Tp - (ct1 * (thetappr - theta1zpr) + kt1 * (thetap - theta1z))) / Jp;

	double xb1pr2 = (Fbx1 - (csx1 * (xb1pr - x1pr) + cb1 * xb1pr + ksx1 * (xb1 - x1))) / mb1;
	double yb1pr2 = (-mb1 * g + Fby1 - (csy1 * (yb1pr - y1pr) + cb1 * yb1pr + ksy1 * (yb1 - y1))) / mb1;
	double zb1pr2 = (Fbz1 - (csz1 * (zb1pr - z1pr) + cb1 * zb1pr + ksz1 * (zb1 - z1))) / mb1;
	double faix1pr2 = (Mx1 - Iy1 * omega1 * faiy1pr) / Ix1;
	double faiy1pr2 = (My1 - Iy1 * omega1 * faix1pr) / Ix1;

	//%% bear2

	double xb2pr2 = (Fbx2 - (csx2 * (xb2pr - x1pr) + cb2 * xb2pr + ksx2 * (xb2 - x1))) / mb2;
	double yb2pr2 = (-mb2 * g + Fby2 - (csy2 * (yb2pr - y1pr) + cb2 * yb2pr + ksy2 * (yb2 - y1))) / mb2;
	double zb2pr2 = (Fbz2 - (csz2 * (zb2pr - z1pr) + cb2 * zb2pr + ksz2 * (zb2 - z1))) / mb2;
	double faix2pr2 = (Mx2 - Iy2 * omega1 * faiy2pr) / Ix2;
	double faiy2pr2 = (My2 - Iy2 * omega1 * faix2pr) / Ix2;


	double xb3pr2 = (Fbx3 - (csx3 * (xb3pr - x2pr) + cb3 * xb3pr + ksx3 * (xb3 - x2))) / mb3;
	double yb3pr2 = (-mb3 * g + Fby3 - (csy3 * (yb3pr - y2pr) + cb3 * yb3pr + ksy3 * (yb3 - y2))) / mb3;
	double zb3pr2 = (Fbz3 - (csz3 * (zb3pr - z2pr) + cb3 * zb3pr + ksz3 * (zb3 - z2))) / mb3;
	double faix3pr2 = (Mx3 - Iy3 * omega2 * faiy3pr) / Ix3;
	double faiy3pr2 = (My3 - Iy3 * omega2 * faix3pr) / Ix3;

	////%% bear4

	double xb4pr2 = (Fbx4 - (csx4 * (xb4pr - x2pr) + cb4 * xb4pr + ksx4 * (xb4 - x2))) / mb4;
	double yb4pr2 = (-mb4 * g + Fby4 - (csy4 * (yb4pr - y2pr) + cb4 * yb4pr + ksy4 * (yb4 - y2))) / mb4;
	double zb4pr2 = (Fbz4 - (csz4 * (zb4pr - z2pr) + cb4 * zb4pr + ksz4 * (zb4 - z2))) / mb4;
	double faix4pr2 = (Mx4 - Iy4 * omega2 * faiy4pr) / Ix4;
	double faiy4pr2 = (My4 - Iy4 * omega2 * faix4pr) / Ix4;

	double thetagpr2 = (-Tg - (ct2 * (thetagpr - theta2zpr) + kt2 * (thetag - theta2z))) / Jg;
	//%% gear

	double omegam = 2. * M_PI * 500. * Z1 / 60;
	double faim = 0.;
	double et = em + er * sin(omegam * t + faim);

	//%% 主动
	//	% deltat = (R1 * theta1z + R2 * theta2z) * cos(betab)...
	//	% -(R1 * theta1y + R2 * theta2y) * sin(betab) * cos(at)...
	//	% +(R1 * theta1x + R2 * theta2x) * sin(betab) * sin(at)...
	//	% +((x1 + ru1 * cos(fai1)) - (x2 + ru2 * cos(fai2))) * cos(al - at)...
	//	% +((x1 - ru1 * sin(fai1)) - (x2 + ru2 * sin(fai2))) * sin(al - at)...
	//	% +(z1 - z2) * tan(beta) - et;
	//%
	//	% deltatpr = (R1 * theta1zpr + R2 * theta2zpr) * cos(betab)...
	//	% -(R1 * theta1ypr + R2 * theta2ypr) * sin(betab) * cos(at)...
	//	% +(R1 * theta1xpr + R2 * theta2xpr) * sin(betab) * sin(at)...
	//	% +((x1pr + ru1 * cosfai1pr) - (x2pr + ru2 * cosfai2pr)) * cos(al - at)...
	//	% +((x1pr - ru1 * sinfai1pr) - (x2pr + ru2 * sinfai2pr)) * sin(al - at)...
	//	% +(z1pr - z2pr) * tan(beta) - et;

	double deltat = (R1 * theta1z - R2 * theta2z) * tan(beta)\
		- (R1 * theta1y - R2 * theta2y) * sin(al - at)\
		+ (R1 * theta1x - R2 * theta2x) * cos(al - at)\
		+ ((x1 + ru1 * cos(fai1)) - (x2 + ru2 * cos(fai2))) * cos(al - at)\
		+ ((y1 - ru1 * sin(fai1)) - (y2 + ru2 * sin(fai2))) * sin(al - at)\
		+ (z1 - z2) * tan(beta) - et;

	double deltatpr = (R1 * theta1zpr - R2 * theta2zpr) * tan(beta)\
		- (R1 * theta1ypr - R2 * theta2ypr) * sin(al - at)\
		+ (R1 * theta1xpr - R2 * theta2xpr) * cos(al - at)\
		+ ((x1pr + ru1 * cosfai1pr) - (x2pr + ru2 * cosfai2pr)) * cos(al - at)\
		+ ((y1pr - ru1 * sinfai1pr) - (y2pr + ru2 * sinfai2pr)) * sin(al - at)\
		+ (z1pr - z2pr) * tan(beta) - et;

	//if (deltat > 1 || deltat < -1) {
	//	deltat = deltat;
	//}

	double Fm = km * f(deltat, b) + cm * deltatpr;

	// test no gear
	if (flag_gear_off) {
		Fm = 0.;
	}

	//Fm = 1000.;
	double Fx = -Fm * cos(al - at);
	double Fy = -Fm * sin(al - at);
	double Fz = -Fm * tan(beta);

	arma::mat LEA = { {m1,						0,						m1 * ru1 * sin(fai1)},
				{0,							m1,						m1 * ru1 * cos(fai1)},
				{-m1 * ru1 * sin(fai1),		m1 * ru1 * sin(fai1),   J1 + m1 * pow(ru1, 2)	} };

	double dd1 = -Fx + m1 * ru1 * pow((omega1 + theta1zpr), 2) * cos(fai1)\
		- (csx1 * (x1pr - xb1pr) + csx2 * (x1pr - xb2pr) + ksx1 * (x1 - xb1) + ksx2 * (x1 - xb2));
	double dd2 = -m1 * g - Fy + m1 * ru1 * pow((omega1 + theta1zpr), 2) * cos(fai1)\
		- (csy1 * (y1pr - yb1pr) + csy2 * (y1pr - yb2pr) + ksy1 * (y1 - yb1) + ksy2 * (y1 - yb2));
	double dd3 = -Fm * R1 - ct1 * (theta1zpr - thetappr) - kt1 * (theta1z - thetap);

	arma::vec LEb;
	LEb << dd1 << dd2 << dd3;
	arma::vec LEx = arma::solve(LEA, LEb);

	double x1pr2 = LEx[0];
	double y1pr2 = LEx[1];
	double theta1zpr2 = LEx[2];

	double theta1xpr2 = (-Fm * R1 - cn1 * theta1xpr - kn1 * theta1x) / J1;
	double theta1ypr2 = (-Fm * R1 - cn1 * theta1ypr - kn1 * theta1y) / J1;
	double z1pr2 = (-Fz - (csz1 * (z1pr - zb1pr) + csz2 * (z1pr - zb2pr) + ksz1 * (z1 - zb1) + ksz2 * (z1 - zb2))) / m1;

	//%% 从动
	arma::mat LEA2 = { {m2,               0,                   m2 * ru2 * sin(fai2)},
	{0,                   m2,                  m2 * ru2 * cos(fai2)},
	{-m2 * ru2 * sin(fai2),   m2 * ru1 * sin(fai2),    J2 + m2 * pow(ru2, 2)} };

	dd1 = Fx + m2 * ru2 * pow((omega2 + theta2zpr), 2) * cos(fai2)\
		- (csx3 * (x2pr - xb3pr) + csx4 * (x2pr - xb4pr) + ksx3 * (x2 - xb3) + ksx4 * (x2 - xb4));
	dd2 = -m2 * g + Fy - m2 * ru2 * pow((omega2 + theta2zpr), 2) * cos(fai2)\
		- (csy3 * (y2pr - yb3pr) + csy4 * (y2pr - yb4pr) + ksy3 * (y2 - yb3) + ksy4 * (y2 - yb4));
	dd3 = Fm * R2 - ct2 * (theta2zpr - thetagpr) - kt2 * (theta2z - thetag);

	arma::vec LEb2;
	LEb2 << dd1 << dd2 << dd3;
	arma::vec LEx2 = arma::solve(LEA2, LEb2);

	double x2pr2 = LEx2[0];
	double y2pr2 = LEx2[1];
	double theta2zpr2 = LEx2[2];

	double theta2xpr2 = (Fm * R2 - cn2 * theta1xpr - kn2 * theta2x) / J2;
	double theta2ypr2 = (Fm * R2 - cn2 * theta1ypr - kn2 * theta2y) / J2;
	double z2pr2 = (Fz - (csz3 * (z2pr - zb3pr) + csz4 * (z2pr - zb4pr) + ksz3 * (z2 - zb3) + ksz4 * (z2 - zb4))) / m2;

	if (abs(z1 + z2) > 1) {
		dxdt[1] = X[2];
	}


	//%% map
	if (flag_bearing_off) {
		xb1pr2 = xb2pr2 = xb3pr2 = xb4pr2 = 0;
		yb1pr2 = yb2pr2 = yb3pr2 = yb4pr2 = 0;
		zb1pr2 = zb2pr2 = zb3pr2 = zb4pr2 = 0;
		faix1pr2 = faix2pr2 = faix3pr2 = faix4pr2 = 0;
		faiy1pr2 = faiy2pr2 = faiy3pr2 = faiy4pr2 = 0;
	}


	dxdt[1] = X[2];
	dxdt[2] = thetappr2;
	dxdt[3] = X[4];
	dxdt[4] = xb1pr2;
	dxdt[5] = X[6];
	dxdt[6] = yb1pr2;
	dxdt[7] = X[8];
	dxdt[8] = zb1pr2;
	dxdt[9] = X[10];
	dxdt[10] = faix1pr2;
	dxdt[11] = X[12];
	dxdt[12] = faiy1pr2;
	dxdt[13] = X[14];
	dxdt[14] = x1pr2;
	dxdt[15] = X[16];
	dxdt[16] = y1pr2;
	dxdt[17] = X[18];
	dxdt[18] = z1pr2;
	dxdt[19] = X[20];
	dxdt[20] = theta1xpr2;
	dxdt[21] = X[22];
	dxdt[22] = theta1ypr2;
	dxdt[23] = X[24];
	dxdt[24] = theta1zpr2;
	dxdt[25] = X[26];
	dxdt[26] = xb2pr2;
	dxdt[27] = X[28];
	dxdt[28] = yb2pr2;
	dxdt[29] = X[30];
	dxdt[30] = zb2pr2;
	dxdt[31] = X[32];
	dxdt[32] = faix2pr2;
	dxdt[33] = X[34];
	dxdt[34] = faiy2pr2;
	dxdt[35] = X[36];
	dxdt[36] = xb3pr2;
	dxdt[37] = X[38];
	dxdt[38] = yb3pr2;
	dxdt[39] = X[40];
	dxdt[40] = zb3pr2;
	dxdt[41] = X[42];
	dxdt[42] = faix3pr2;
	dxdt[43] = X[44];
	dxdt[44] = faiy3pr2;
	dxdt[45] = X[46];
	dxdt[46] = x2pr2;
	dxdt[47] = X[48];
	dxdt[48] = y2pr2;
	dxdt[49] = X[50];
	dxdt[50] = z2pr2;
	dxdt[51] = X[52];
	dxdt[52] = theta2xpr2;
	dxdt[53] = X[54];
	dxdt[54] = theta2ypr2;
	dxdt[55] = X[56];
	dxdt[56] = theta2zpr2;
	dxdt[57] = X[58];
	dxdt[58] = xb4pr2;
	dxdt[59] = X[60];
	dxdt[60] = yb4pr2;
	dxdt[61] = X[62];
	dxdt[62] = zb4pr2;
	dxdt[63] = X[64];
	dxdt[64] = faix4pr2;
	dxdt[65] = X[66];
	dxdt[66] = faiy4pr2;
	dxdt[67] = X[68];
	dxdt[68] = thetagpr2;

	if (flag_zero_fai) {
		dxdt[9] = 0;
		dxdt[10] = 0;
		dxdt[11] = 0;
		dxdt[12] = 0;
		dxdt[31] = 0;
		dxdt[32] = 0;
		dxdt[33] = 0;
		dxdt[34] = 0;
		dxdt[41] = 0;
		dxdt[42] = 0;
		dxdt[43] = 0;
		dxdt[44] = 0;
		dxdt[63] = 0;
		dxdt[64] = 0;
		dxdt[65] = 0;
		dxdt[66] = 0;
	}
	if (flag_zero_thetag) {
		dxdt[67] = 0;
		dxdt[68] = 0;
	}
	if (flag_zero_thetap) {
		dxdt[1] = 0;
		dxdt[2] = 0;
	}





	//for (int i = 1; i <= 67; i = i + 2) {
	//	if (X[i] > 10) {
	//		std::cout<<"?";
	//	}
	//}
	for (int i = 0; i <= 68; i++) {
		if (isnan(dxdt[i]) || isinf(dxdt[i]))
		{
			dxdt[0] = 0;
		}
	}
}


vecd bearCalc0(
	double t, double  xbi, double  ybi, double zbi, double faixi, double faiyi,
	double Rj, double A0, double alp0, double Nb, double wcage, double Kc, double Pd,
	int install_disp) {
	double Fbxi = 0;
	double Fbyi = 0;
	double Fbzi = 0;
	double Mxi = 0;
	double Myi = 0;
	vecd theta(Nb + 1); //= zeros(Nb, 1);% t的函数 长Nb的向量
	vecd deltasum(Nb + 1);// = zeros(Nb, 1); % 5自由度和t的函数 长Nb的向量
	vecd Q(Nb + 1);// = zeros(Nb, 1); % 5自由度和t的函数 长Nb的向量
	vecd fai_r_theta(Nb + 1);// = zeros(Nb, 1);
	vecd deltazj(Nb + 1);// = zeros(Nb, 1);
	vecd deltayj(Nb + 1);// = zeros(Nb, 1);
	vecd A0pr(Nb + 1);// = zeros(Nb, 1);
	vecd alp0pr(Nb + 1);// = zeros(Nb, 1);
	vecd Hp(Nb + 1);// = zeros(Nb, 1);
	double nA = sqrt(pow(xbi, 2) + pow(ybi, 2));
	if (nA == 0) {
		;
	}
	vecd x111(Nb + 1);
	vecd y111(Nb + 1);
	for (int i = 1; i <= Nb; i++) {
		theta[i] = wcage * t + 2. * M_PI / Nb * (i - 1);
		//deltazj[i] = zbi + Rj * (faixi * sin(theta[i]) - faiyi * cos(theta[i]));
		deltazj[i] = zbi + Rj * (faixi * sin(theta[i]) - faiyi * cos(theta[i]));
		//%% changed method
			//% deltayj = xbi * cos(theta[i]) + ybi * sin(theta[i]);% substituate
			//% A0pr = sqrt((A0 * sin(alp0) + deltazj). ^ 2 + (A0 * cos(alp0) + deltayj). ^ 2);
		//% deltasum[i] = A0pr - A0;
		

		//double nB = sqrt(pow(cos(theta[i]), 2) + pow(sin(theta[i]), 2));
		x111[i] = xbi * cos(theta[i]);
		y111[i] = ybi * sin(theta[i]);

		fai_r_theta[i] = acos((xbi * cos(theta[i]) + ybi * sin(theta[i])) / (nA * 1));
		if (isnan(fai_r_theta[i]) || isinf(fai_r_theta[i]))
		{
			//fai_r_theta[i]=0;
		}
		if ((xbi) == 0 && ybi == 0) {
			fai_r_theta[i] = theta[i];
		}
		if ((xbi * cos(theta[i]) + ybi * sin(theta[i])) > nA)
		{
			fai_r_theta[i] = 0;
		}else if ((xbi * cos(theta[i]) + ybi * sin(theta[i])) <- nA)
		{
			fai_r_theta[i] = M_PI;
		}



		deltayj[i] = nA * cos(fai_r_theta[i]) - 0.5 * Pd;

		//% if deltayj[i] < 0;
		//% deltayj[i] = 0;
		//% end
		A0pr[i] = sqrt(pow((A0 * sin(alp0) + deltazj[i]), 2) + pow((A0 * cos(alp0) + deltayj[i]), 2));
		deltasum[i] = A0pr[i] - A0;

		//% contact critia 修角度！！！！！！！！！！！！！！！！！！！！！！！！
		double tanalp0pr = (A0 * sin(alp0) + deltazj[i]) / (A0 * cos(alp0) + deltayj[i]);
		alp0pr[i] = atan(tanalp0pr) * (-1) * install_disp;


		Hp[i] = 0;

		if (install_disp == 1)// % z + Fz - (z - A = 0)
		{

			if ((A0 * sin(alp0) + deltazj[i]) <= 0)
			{
				alp0pr[i] = 0;
			}
		}
		else//% z - Fz + (z + A = 0)
		{

			if ((-A0 * sin(alp0) + deltazj[i]) >= 0)
			{
				alp0pr[i] = 0;
			}
		}
		if (A0 * cos(alp0) + deltayj[i] >= 0 && deltasum[i] >= 0)
		{
			Q[i] = Kc * pow(deltasum[i] * 1000., 1.5);
		}
		else {
			Q[i] = 0;
		}



		Fbxi = Fbxi + Q[i] * cos(alp0pr[i]) * cos(theta[i] + M_PI);
		Fbyi = Fbyi + Q[i] * cos(alp0pr[i]) * sin(theta[i] + M_PI);
		Fbzi = Fbzi + Q[i] * sin(alp0pr[i]) * install_disp;
		Mxi = Mxi + Rj * Q[i] * sin(alp0pr[i]) * sin(theta[i]);
		Myi = Myi -/*-*/ Rj * Q[i] * sin(alp0pr[i]) * cos(theta[i] );


	}
	vecd ret(5);
	ret[0] = Fbxi;
	ret[1] = Fbyi;
	ret[2] = Fbzi;

	//if (Mxi * faixi > 0) {
	//	Mxi = -Mxi;
	//}


	//if (Myi * faiyi > 0) {
	//	Myi = -Myi;
	//}


	ret[3] = Mxi;
	ret[4] = Myi;

	for (int i = 0; i <= 4; i++) {
		if (isnan(ret[i])||isinf(ret[i]))
		{
			ret[0] = Fbxi;
		}
	}



	return ret;

}


double f(double delta, double b) {
	double ret;
	if(delta > b )
	{
		ret = delta - b;
	}else if (delta < -b)
	{	
		ret = delta + b;
	}
	else {
		ret = 0;
	}
	return ret;
}
