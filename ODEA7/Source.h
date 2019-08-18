#pragma once
#include <thread>
#include "Bear.h"
typedef std::vector< double > vecd;
/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;
void WriteFile(std::vector< state_type > states, std::vector< double > times, long long step, long long fileNum);



double f(double delta, double b);
//[ integrate_observer
struct push_back_state_and_time
{
	std::vector< state_type >& m_states;
	std::vector< double >& m_times;
	long long step = 0;
	long long file_step = 1;
	int fileNUM = 0;
	push_back_state_and_time(std::vector< state_type >& states, std::vector< double >& times)
		: m_states(states), m_times(times) { }

	void operator()(const state_type& x, double t)
	{
		m_states.push_back(x);
		m_times.push_back(t);

		if (step % 100000 == 0)
		{
			std::cout << m_times[file_step - 1] << '\t' << m_states[file_step - 1][0] << '\t' << m_states[file_step - 1][1] << '\n';
			std::vector< state_type > new_states;
			std::vector< double > new_times;

			m_states.swap(new_states);
			m_times.swap(new_times);
			std::thread second(WriteFile, new_states, new_times, step, fileNUM);
			second.detach();
			file_step = 0;
			m_states.clear();
			m_times.clear();
			fileNUM++;
		}

		step++;
		file_step++;

		//if (step == 25643) {
		//	step=step;
		//}
	}


};
//]

struct write_state
{
	void operator()(const state_type& x) const
	{
		std::cout << x[0] << "\t" << x[1] << "\n";
	}
};

class ODE_MODEL {
public:
	double g = 9.8;
	double kt1 = 8e8; double  kt2 = 1.5e8;
	double ct1 = 400; double  ct2 = 400;
	double kn1 = 8e8; double kn2 = 8e8; double cn1 = 400; double cn2 = 400;
	double cb1 = 500; double cb2 = 500; double cb3 = 500; double cb4 = 500;
	double etat = 0.1;
	double etas = 0.1;
	double ksx1 = 6e8; double ksy1 = 6e8; double ksz1 = 6e8; double ksx2 = 6e8; double ksy2 = 6e8; double ksz2 = 6e8;
	double ksx3 = 0.16e8; double ksy3 = 0.16e8; double ksz3 = 0.16e8; double ksx4 = 0.16e8; double ksy4 = 0.16e8; double ksz4 = 0.16e8;
	double csx1 = 500; double csy1 = 500; double csz1 = 500; double csx2 = 500; double csy2 = 500; double csz2 = 500;
	double csx3 = 400; double csy3 = 400; double csz3 = 400; double csx4 = 400; double csy4 = 400; double csz4 = 400;
	double Ix1 = 5e1; double Iy1 = 5e1;
	double Ix2 = 5e1; double Iy2 = 5e1;
	double Ix3 = 5e1; double Iy3 = 5e1;
	double Ix4 = 5e1; double Iy4 = 5e1;
	double mb1 = 115.2; double mb2 = 115.2; double mb3 = 7.32; double mb4 = 7.32;
	double m1 = 667; double m2 = 141;
	double Jp = 20; double Jg = 5; double J1 = 20; double J2 = 20;
	double r = 20;
	double R = 47;
	double fo = 0.515;
	double fi = 0.525;
	double Db = 10. / 1000.; //% ¹öÖéÖ±¾¶ unit : m
	double Pd = 0.015 / 1000; //% unit:m
	double Tdm = 400;
	double Tdr = 300;
	double Tgm = 400;
	double Tgr = 300;
	int	Nb = 8;
	//beta = 15 / 180 * pi;
	//betab = 13.8 / 180 * pi;
	double beta = 0.;
	double betab = 0.;
	double an = 20. / 180. * M_PI;
	double at = atan(tan(an) / cos(beta));
	double al = 120. / 180. * M_PI;
	double xig = 0.1;
	double em = 1e-5;
	double er = 1e-5;
	double b = 4e-5;
	double ru1 = 2e-5;
	double ru2 = 2e-5;

	double bc = 1.65e-2;
	double Z1 = 30;
	double Z2 = 30;

	double cm = 1.2e3;
	double km = 500e6;
	/*%% m = 2£¿*/
	double R1 = 0.03;
	double R2 = 0.03;

	double Rj;
	double A0;
	double alp0;
	double dm;
	double Kc;


	bool flag_zero_thetap = false;
	bool flag_zero_thetag = false;
	bool flag_zero_fai = false;
	bool flag_gear_off = false;
	bool flag_bearing_off = false;
	bool flag_zero_cage_movement = false;
	bool flag_limit_fai = false;

public:
	ODE_MODEL();
	void readini();
	double limit_fai(double fai);
	void operator() (const state_type& X, state_type& dxdt, const double  t);
};


