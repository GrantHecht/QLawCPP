// Coded by Donald Ellison
// March 5th 2014

#include "qlaw.h"

void gaussvar(double X[], double dX[], double L, double & stepsize, double & L_bin, int & throttle_lock, int & coast, options * probdata)
{
	// state vector elements
	double p = X[0];
	double e = X[1];
	double inc = X[2];
	double ape = X[3];
	double ran = X[4];
	double mass = X[5];

	// other orbital parameters
	double tru = L - ran - ape;
	double sma = p / (1.0 - e*e);
	double h = sqrt(probdata->mu*p);
	double r = p / (1.0 + e*cos(tru));
	double g = e*sin(ape + ran);
	double w = 1.0 + e*cos(ape + ran)*cos(L) + g*sin(L);
	double k = tan(inc / 2.0)*sin(ran);

	double current_thrust = probdata->thrust;
	double f = current_thrust / mass; //thrust acceleration
	double fr, ftheta, fh;


	// "Optimal" Lyapunov Control
	double u[3];
	lyapunov(u, sma, e, inc, ape, ran, tru, f, probdata);
	//lyapunov_rapomax(u, sma, e, inc, ape, ran, tru, f, probdata);

	// Lyapunov does not result in physical control values, so normalize to current thrust level
	double normu = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
	fr = u[0] / normu*f;
	ftheta = u[1] / normu*f;
	fh = u[2] / normu*f;

	//Effectivity Calculation (coasting) 
	//This is REALLY slow right now, I need to think of a faster way to implement this
	//In MATLAB this is trivial (dot operator)

	if (probdata->effectivity)
	{
		
		//Once we turn the thruster on we want it to fire for a minimum
		//of 10.0 degrees in true longitude to avoid chattering
		if (L_bin * 180.0 / probdata->APPLE_PI > 10.0)
		{
			L_bin = 0.0;
			throttle_lock = 0;
		}
		
		double eta_abs = 0.0, eta_rel = 0.0;

		//calculate how effective thrusting at this location on the osculating orbit
		//is at driving Q to zero
		effectivity(eta_abs, eta_rel, sma, e, inc, ape, ran, tru, f, L_bin, probdata);

		//if we are not very effective then turn the thruster off
		//the minimum 10.0 degree rule overrides this though
		if (probdata->rel_effectivity)
		{
			if (eta_rel < probdata->eta_rel_tol && throttle_lock == 0)
			{
				fr = 0.0;
				ftheta = 0.0;
				fh = 0.0;
				current_thrust = 0.0;
				coast = 1;
			}
			else
				coast = 0;
		}
		
		if (probdata->abs_effectivity)
		{
			if (eta_abs < probdata->eta_abs_tol && throttle_lock == 0)
			{
				fr = 0.0;
				ftheta = 0.0;
				fh = 0.0;
				current_thrust = 0.0;
				coast = 1;
			}
			else
				coast = 0;
		}

		//If we don't coast, lock the thruster until 10.0 degrees have passed
		if (coast == 0)
			throttle_lock = 1;
	}
	
	// Gauss Variational Equations
	// Integrate w.r.t. true longitude L
	double dL_dt = sqrt(probdata->mu*p)*(w / p)*(w / p) + 1.0 / w*sqrt(p / probdata->mu)*(tan(inc / 2.0)*cos(ran)*sin(L) - k*cos(L))*fh;
	double dt_dL = 1.0 / dL_dt;
	
	// Calculate the state vector gradient w.r.t. time
	dX[0] = 2.0*p / w*sqrt(p / probdata->mu)*ftheta;
	dX[1] = 1.0 / h*(p*sin(tru)*fr + ((p + r)*cos(tru) + r*e)*ftheta);
	dX[2] = r*cos(tru + ape) / h*fh;
	dX[3] = 1.0 / (e*h)*(-p*cos(tru)*fr + (p + r)*sin(tru)*ftheta) - (r*sin(tru + ape)*cos(inc)) / (h*sin(inc))*fh;
	dX[4] = r*sin(tru + ape) / (h*sin(inc))*fh;
	dX[5] = -current_thrust / (probdata->Isp*probdata->grav);

	// Convert integration variable from time to true longitude
	for (size_t j = 0; j < 6; ++j)
		dX[j] = dX[j]*dt_dL;
	
	// Add time derivative to the gradient vector so that it may be integrated along with the COE's
	dX[6] = dt_dL;

}