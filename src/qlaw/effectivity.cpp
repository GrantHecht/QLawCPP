#include "qlaw.h"

void effectivity(double & eta_abs, double & eta_rel, double & sma, double & e, double & inc, double & ape, double & ran, double & tru, double & f, double L_bin, options * probdata)
{

	double Qdot_n, Qdot_test, Qdot_nn, Qdot_nx;
	double trutest = 0.0;

	Qdot_n = minQdot(sma, e, inc, ape, ran, tru, f, probdata);
	Qdot_nn = Qdot_n;
	Qdot_nx = 0.0;

	while (trutest <= 2.0 * probdata->APPLE_PI)
	{
		Qdot_test = minQdot(sma, e, inc, ape, ran, trutest, f, probdata);
		if (Qdot_test < Qdot_nn)
			Qdot_nn = Qdot_test;
		else if (Qdot_test > Qdot_nx)
			Qdot_nx = Qdot_test;
		trutest = trutest + probdata->effectivity_sample; //sample the osculating orbit every so often
	}

	if (Qdot_n - Qdot_nn < 1.0e-10)
		eta_abs = 1.0;
	else
		eta_abs = Qdot_n / Qdot_nn;
	
	if ((Qdot_n - Qdot_nx) - (Qdot_nn - Qdot_nx) < 1.0e-10)
		eta_rel = 1.0;
	else
		eta_rel = (Qdot_n - Qdot_nx) / (Qdot_nn - Qdot_nx);
	
}