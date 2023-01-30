#include "qlaw.h"

bool check_convergence(double * X_current, options * probdata)
{
	
	bool smacon = false;
	bool econ   = false;
	bool inccon = false;
	bool apecon = false;
	bool rancon = false;

	// sma is not in state, p is
	double sma_current = X_current[0] / (1.0 - X_current[1] * X_current[1]);

	// check the convergence of each COE
	// if a COE's final value is free then set its convergence flag to true
	if (probdata->sma_switch == 0)
		smacon = true;
	else if (fabs(sma_current - probdata->sma_t) <= probdata->sma_tol)
		smacon = true;

	if (probdata->e_switch == 0)
		econ = true;
	else if (fabs(X_current[1] - probdata->e_t) <= probdata->e_tol)
		econ = true;

	if (probdata->inc_switch == 0)
		inccon = true;
	else if (fabs(X_current[2] - probdata->inc_t) <= probdata->inc_tol)
		inccon = true;

	if (probdata->ape_switch == 0)
		apecon = true;
	else if (X_current[3] < 0.0)
	{
		if(fabs(X_current[3] + 2.0*probdata->APPLE_PI - probdata->ape_t) <= probdata->ape_tol)
			apecon = true;
	}
	else if (X_current[3] >= 0.0)
	{
		if (fabs(X_current[3] - probdata->ape_t) <= probdata->ape_tol)
			apecon = true;
	}


	if (probdata->ran_switch == 0)
		rancon = true;
	else if (X_current[4] < 0.0)
	{
		if (fabs(X_current[4] + 2.0*probdata->APPLE_PI - probdata->ran_t) <= probdata->ran_tol)
			apecon = true;
	}
	else if (X_current[4] >= 0.0)
	{
		if (fabs(X_current[4] - probdata->ran_t) <= probdata->ran_tol)
			apecon = true;
	}

	// if all convergence checks pass then the target orbit has been reached to within tolerences
	if (smacon && econ && inccon && apecon && rancon)
		return true;
	else
		return false;
}