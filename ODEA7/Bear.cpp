#include "Bear.h"
#include<math.h>
#include <corecrt_math_defines.h>

void bear::calc()
{
	vecd BearCalc=m_bearCalc.bearCalc_20190807_withDiag( t,  xbi,  ybi, zbi,  faixi,  faiyi,  Rj, A0, alp0, Nb, wcage, Kc, Pd, dm, Db, install_disp,  bearing_num);
	Fbx = BearCalc[0];
	Fby = BearCalc[1];
	Fbz = BearCalc[2];
	Mx = BearCalc[3];
	My = BearCalc[4];
}

vecd bearCalc::bearCalc_20190807_withDiag(double t, double xbi, double ybi, double zbi, double faixi, double faiyi, double Rj, double A0, double alp0, double Nb, double wcage, double Kc, double Pd, double dm, double Db, int install_disp, int bearing_num)
{
	double Fbxi = 0;
	double Fbyi = 0;
	double Fbzi = 0;
	double Mxi = 0;
	double Myi = 0;

	//double nA = sqrt(pow(xbi, 2) + pow(ybi, 2));
	if (abs(zbi) >= 0.1) {
		zbi = zbi;
	}
	double douter = dm + Db + Pd / 2;
	double rinner = (dm - Db - Pd / 2) / 2;
	for (int i = 1; i <= Nb; i++) {
		theta[i] = wcage * t + 2. * M_PI / Nb * (i - 1);
		double th = theta[i];
		double costh = cos(th);
		double sinth = sin(th);
		deltazj[i] = zbi + Rj * (-faixi * sinth + faiyi * costh);


		double root;
		double d = sqrt(pow2(rinner) - pow2(ybi * costh - xbi * sinth));
		//double root1 = pow((pow(rinner * cos(th), 2) - pow(ybi * cos(th), 2) + pow(rinner * sin(th), 2) - pow(xbi * sin(th), 2) + 2 * xbi * ybi * cos(th) * sin(th)), (0.5))
		//	+ xbi * cos(th) + ybi * sin(th);
		double root1 = d + xbi * costh + ybi * sinth;
		double root2 = -d + xbi * costh + ybi * sinth;
		if (root1 * root2 < 0) {
			root = root1 < 0 ? root2 : root1;
		}
		else {
			root = +xbi * costh + ybi * sinth;// check move too extremely?
		}
		rmovement[i] = root - rinner;
		deltayj[i] = rmovement[i] + A0 - Pd / 2;
		//deltayj[i] = rmovement[i]+A0 - Pd/2+Rj * (+faixi * sin(theta[i]) - faiyi * cos(theta[i]));
		if (deltayj[i] >= 0) {
			if (install_disp * deltazj[i] <= 0)//radical model
			{
				A0pr[i] = deltayj[i];
				deltasum[i] = A0pr[i] - A0;
				if (deltasum[i] > 0) {
					//A0pr[i]=sqrt(pow(( deltazj[i]), 2) + pow(( deltayj[i]), 2));
					Q[i] = Kc * pow(deltasum[i] * 1000., 1.5);
					alp0pr[i] = 0;
				}
				else {
					Q[i] = 0;
					alp0pr[i] = 0;
				}
			}
			else if (deltayj[i] > 0)			//angular model
			{
				A0pr[i] = sqrt(pow2((deltazj[i])) + pow2((deltayj[i])));
				deltasum[i] = A0pr[i] - A0;
				if (deltasum[i] > 0) {
					Q[i] = Kc * pow(deltasum[i] * 1000., 1.5);
					double tanalp0pr = deltazj[i] / deltayj[i];
					alp0pr[i] = atan(tanalp0pr);
				}
				else {
					Q[i] = 0;
				}

			}
			else {//push
				A0pr[i] = abs(deltazj[i]);
				deltasum[i] = A0pr[i] - A0;
				if (deltasum[i] > 0)
				{
					Q[i] = Kc * pow(deltasum[i] * 1000., 1.5);
					alp0pr[i] = -M_PI_2 * install_disp;
				}
				else {
					Q[i] = 0;
				}
			}
		}
		else //deltayj[i]<0**********
		{
			A0pr[i] = abs(deltazj[i]);
			deltasum[i] = A0pr[i] /*- A0*/;

			if (deltazj[i] * install_disp > 0)
			{
				Q[i] = Kc * pow(deltasum[i] * 1000., 1.5);
				alp0pr[i] = -M_PI_2 * install_disp;
			}
			else
			{

				Q[i] = 0;//alp0pr=0
			}

		}

		//if (install_disp == 1)// % z + Fz - (z - A = 0)
		//else//% z - Fz + (z + A = 0)
		Fbxi_[i] = -Q[i] * cos(alp0pr[i]) * costh;
		Fbyi_[i] = -Q[i] * cos(alp0pr[i]) * sinth;
		Fbzi_[i] = -Q[i] * sin(alp0pr[i]);
		Mxi_[i] = -(-Rj * Q[i] * abs(sin(alp0pr[i])) * sinth) * install_disp;
		Myi_[i] = - Rj * Q[i] * abs(sin(alp0pr[i])) * costh * install_disp;

		Fbxi = Fbxi + Fbxi_[i];

		Fbyi = Fbyi + Fbyi_[i];
		Fbzi = Fbzi + Fbzi_[i];
		Mxi = Mxi + Mxi_[i];
		Myi = Myi + Myi_[i];
	}
	vecd ret(5);

	//if (Mxi * faixi > 0) { Mxi = -Mxi; }
	//if (Myi * faiyi > 0) { Myi = -Myi; }

	ret[0] = Fbxi;
	ret[1] = Fbyi;
	ret[2] = Fbzi;
	ret[3] = Mxi;
	ret[4] = Myi;


	if (zbi * Fbzi > 0 || xbi * Fbxi > 0 || ybi * Fbyi > 0) {
		ret[4] = ret[4];
	}
	if (abs(faixi) > 0.7 || abs(faiyi) > 0.7) {
		ret[4] = ret[4];
	}
	if (abs(xbi) > 0.1 || abs(faiyi) > 0.1) {
		ret[4] = ret[4];
	}



	for (int i = 0; i <= 4; i++) {
		if (isnan(ret[i]) || isinf(ret[i]))
		{
			ret[0] = Fbxi;
		}
	}
	return ret;
}
bearCalc bear::m_bearCalc=bearCalc(20);