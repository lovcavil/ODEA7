#include <boost/numeric/odeint.hpp>
#include <corecrt_math_defines.h>
#include <chrono>
#include <fstream>
#include "Source.h"
typedef std::vector< double > vecd;
/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

using namespace std;
using namespace boost::numeric::odeint;
using namespace std::chrono;

void test_bearing(ODE_MODEL model);

int main(int argc, char** argv)
{

	state_type x(69);

	auto start = system_clock::now();
	double totaltime;
	if (argc == 1) {
		totaltime = 2.;
	}
	else {
		auto c = argv[1];
		totaltime = atof(c);
	}
	cout << "time=" << totaltime << "s" << std::endl;

	ODE_MODEL model = ODE_MODEL();
	model.readini();

	model(x, x,0);
	test_bearing(model);

	//[ state_initialization

	for (int i=0;i<=68;i++)
	{
		x[i] = 0;
	}
	//]

	size_t steps;

	//[ integrate_observ
	vector<state_type> x_vec;
	vector<double> times;

	runge_kutta_cash_karp54 < state_type > stepper;
	push_back_state_and_time obs = push_back_state_and_time(x_vec, times);
	steps=integrate(model, x, 0.0, totaltime, 0.001,obs);

	std::vector< state_type > new_states;
	std::vector< double > new_times;

	x_vec.swap(new_states);
	times.swap(new_times);
	WriteFile(new_states,new_times,99999999999,99999);

	// do something...
	auto end = system_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	cout << "»¨·ÑÁË"
		<< double(duration.count()) * microseconds::period::num / microseconds::period::den << "Ãë" << endl;

	std::ofstream outfile(std::to_string(double(duration.count()) * microseconds::period::num / microseconds::period::den) + ".txt", std::ios::out | std::ios::binary);



}
void vali(double xbi, double ybi, double zbi, double faixi, double faiyi, int disp, ODE_MODEL model);
void test_bearing(ODE_MODEL model) {
/*
	vecd bearCalc_20190807_withDiag(
		double t, double  xbi, double  ybi, double zbi, double faixi, double faiyi,
		double Rj, double A0, double alp0, double Nb, double wcage, double Kc, double Pd, double dm, double Db,
		int install_disp)*/
	
	double Rj = model.Rj;
	double A0 = model.A0;
	double alp0 = model.alp0;
	double Kc = model.Kc;
	double Pd = model.Pd;
	double dm = model.dm;
	double xbi = 0;
	double ybi = 0;
	double zbi = 0;
	double faixi = 0;
	double faiyi = 0;
	vali(xbi, ybi, zbi, faixi, faiyi,-1,model);
	vali(-0.1, 0, 0, 0, 0,-1,model);
	vali(-0.1, 0, 0, 0, 0,1,model);
	vali(0.1, 0, 0, 0, 0, -1, model);
	vali(0.1, 0, 0, 0, 0, 1, model);

	vali(0, -0.1, 0, 0, 0, -1, model);
	vali(0, -0.1, 0, 0, 0, 1, model);
	vali(0, 0.1, 0, 0, 0, -1, model);
	vali(0, 0.1, 0, 0, 0, 1, model);

	vali(0,0, -0.1, 0,  0, -1, model);
	vali(0,0, -0.1, 0,  0, 1, model);
	vali(0,0, 0.1, 0, 0, -1, model);
	vali(0,0, 0.1, 0, 0, 1, model);

	vali(0, 0, 0, -0.1, 0, -1, model);
	vali(0, 0, 0, -0.1, 0, 1, model);
	vali(0, 0, 0, 0.1, 0, -1, model);
	vali(0, 0, 0, 0.1, 0, 1, model);

	vali(0, 0, 0, 0, -0.1, -1, model);
	vali(0, 0, 0, 0, -0.1, 1, model);
	vali(0, 0, 0, 0, 0.1, -1, model);
	vali(0, 0, 0, 0, 0.1, 1, model);


	vecd BearCalc;
	//-0.01039 - 657.59 - 0.012791 - 840.439 - 8.9113 - 126963	1.94059 - 2624.07 - 2.08627	4899.45
	//BearCalc=bearCalc_20190807_withDiag(0, -0.01039, - 0.012791, - 8.9113,1.94059, - 2.08627, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, -1);
	//BearCalc = bearCalc_20190807_withDiag(0.00319008, 0, 0, 0, 0, 1,
	//	Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, 1, 4);
	//BearCalc = bearCalc_20190807_withDiag(0.00319008, -0.000000730201, 0.0000847991, 0.0000308555, 0.928009, 0.0131197,
	//	Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, 1, 4);
	
	//BearCalc=bearCalc_20190807_withDiag(0.00319008, -0.00000127722, 0.0000779963,-0.0000293413, 1.00661, -0.00666389,
	//	Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, -1,1);
	//BearCalc = bearCalc_20190807_withDiag(0.495128, -0.0000677209, 0.000403131, 0.000241044, -0.0068632, -0.000989304,
	//	Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, 1, 4);


		//BearCalc = bearCalc_20190807_withDiag(0.499086,0.229548,-1.24563, -8.63832, 0.000241044, 0.0419844,
		//Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, -1, 1);


	//BearCalc=bearCalc_20190807_withDiag(0, 0, 0, 0, 0.1, 0, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, -1,0);
	//BearCalc=bearCalc_20190807_withDiag(0, 0, 0, 0, -0.1, 0, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, -1,0);
	//BearCalc = bearCalc_20190807_withDiag(0, 0, 0, 0, 0, 0.1, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, -1,0);
	//BearCalc = bearCalc_20190807_withDiag(0, 0, 0, 0, 0, -0.1, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, -1,0);
	//BearCalc = bearCalc_20190807_withDiag(0, 0, 0, 0, 0.1, 0, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, 1,0);
	//BearCalc = bearCalc_20190807_withDiag(0, 0, 0, 0, -0.1, 0, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, 1,0);
	//BearCalc = bearCalc_20190807_withDiag(0, 0, 0, 0, 0, 0.1, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, 1,0);
	//BearCalc = bearCalc_20190807_withDiag(0, 0, 0, 0, 0, -0.1, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, 1,0);

	//BearCalc = bearCalc_20190807_withDiag(0, -0.001, 0, 0, 0, 0, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, -1);


	//BearCalc = bearCalc_20190807_withDiag(0,-0.001, -0.001, 0,  0, 0, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, -1);

	//BearCalc = bearCalc_20190807_withDiag(0,-0.001, -0.001, -0.001, 0, 0, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, -1);

	//BearCalc = bearCalc_20190807_withDiag(0,-0.0001, -0.0001,0.1, 0, 0, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, -1,1);
	//BearCalc = bearCalc_20190807_withDiag(0,-0.0001, -0.0001,0.1, 0, 0, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, 1,2);
	//BearCalc = bearCalc_20190807_withDiag(0,-0.0001, -0.0001,-0.1, 0, 0, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, -1,3);
	//BearCalc = bearCalc_20190807_withDiag(0,-0.0001, -0.0001,-0.1, 0, 0, Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, 1,4);

}

void vali(double xbi, double ybi, double zbi, double faixi, double faiyi,int disp,ODE_MODEL model) {
	double Rj = model.Rj;
	double A0 = model.A0;
	double alp0 = model.alp0;
	double Kc = model.Kc;
	double Pd = model.Pd;
	double dm = model.dm;
	vecd BearCalc = bearCalc(10).bearCalc_20190807_withDiag(0,xbi,  ybi,  zbi,  faixi,  faiyi,
		Rj, A0, alp0, 8, 0, Kc, Pd, dm, 0.01, disp, 1);
	double Fbx = BearCalc[0];
	double Fby = BearCalc[1];
	double Fbz = BearCalc[2];
	double Mx = BearCalc[3];
	double My = BearCalc[4];

	if (xbi*Fbx>0)
	{
		std::cout << "x";
	}
	if (ybi * Fby > 0)
	{
		std::cout << "y";
	}
	if (zbi * Fbz > 0)
	{
		std::cout << "z";
	}
	if (faixi * Mx > 0)
	{
		std::cout << "fx";
	}
	if (faiyi * My > 0)
	{
		std::cout << "fy";
	}


}