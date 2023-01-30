#include "qlaw.h"



bool read_options_file(std::string filename, options * probdata)
{
	std::ifstream settings;

	char file_buffer[256], choice[40];
	std::string stringvalue;
	double value;
	char peek;
	int lbound_index = 0, ubound_index = 0;
	int linenumber = 0;

	settings.open(filename);

	if (settings.is_open())
	{
		while (!settings.eof())
		{
			++linenumber;
			peek = settings.peek();
			if (peek == '#' || peek == '\n' || peek == '\r')
			{
				settings.getline(file_buffer, 256); //grab the entire line and don't care where it goes.
				continue;
			}

			// grab options file identifier
			settings >> choice;

			// grap the value next to the identifier
			settings >> value;

			// problem physical constants and scaling
			if (!strcmp(choice, "mu"))
			{
				probdata->mu = value;
				continue;
			}
			else if (!strcmp(choice, "DU"))
			{
				probdata->DU = value;
				continue;
			}
			else if (!strcmp(choice, "TU"))
			{
				probdata->TU = value;
				continue;
			}
			

			// Spacecraft parameters
			else if (!strcmp(choice, "thrust"))
			{
				probdata->thrust = value;
				continue;
			}
			else if (!strcmp(choice, "Isp"))
			{
				probdata->Isp = value;
				continue;
			}
			else if (!strcmp(choice, "initial_mass"))
			{
				probdata->initial_mass = value;
				continue;
			}

			// initial spacecraft COE's
			else if (!strcmp(choice, "sma0"))
			{
				probdata->sma0 = value;
				continue;
			}
			else if (!strcmp(choice, "e0"))
			{
				probdata->e0 = value;
				continue;
			}
			else if (!strcmp(choice, "inc0"))
			{
				probdata->inc0 = value;
				continue;
			}
			else if (!strcmp(choice, "ape0"))
			{
				probdata->ape0 = value;
				continue;
			}
			else if (!strcmp(choice, "ran0"))
			{
				probdata->ran0 = value;
				continue;
			}
			else if (!strcmp(choice, "tru0"))
			{
				probdata->tru0 = value;
				continue;
			}

			// COE switches
			else if (!strcmp(choice, "sma_switch"))
			{
				probdata->sma_switch = value;
				continue;
			}
			else if (!strcmp(choice, "e_switch"))
			{
				probdata->e_switch = value;
				continue;
			}
			else if (!strcmp(choice, "inc_switch"))
			{
				probdata->inc_switch = value;
				continue;
			}
			else if (!strcmp(choice, "ape_switch"))
			{
				probdata->ape_switch = value;
				continue;
			}
			else if (!strcmp(choice, "ran_switch"))
			{
				probdata->ran_switch = value;
				continue;
			}
			
			// target COE's
			else if (!strcmp(choice, "sma_t"))
			{
				probdata->sma_t = value;
				continue;
			}
			else if (!strcmp(choice, "e_t"))
			{
				probdata->e_t = value;
				continue;
			}
			else if (!strcmp(choice, "inc_t"))
			{
				probdata->inc_t = value;
				continue;
			}
			else if (!strcmp(choice, "ape_t"))
			{
				probdata->ape_t = value;
				continue;
			}
			else if (!strcmp(choice, "ran_t"))
			{
				probdata->ran_t = value;
				continue;
			}

			// convergence tolerances
			else if (!strcmp(choice, "sma_tol"))
			{
				probdata->sma_tol = value;
				continue;
			}
			else if (!strcmp(choice, "e_tol"))
			{
				probdata->e_tol = value;
				continue;
			}
			else if (!strcmp(choice, "inc_tol"))
			{
				probdata->inc_tol = value;
				continue;
			}
			else if (!strcmp(choice, "ape_tol"))
			{
				probdata->ape_tol = value;
				continue;
			}
			else if (!strcmp(choice, "ran_tol"))
			{
				probdata->ran_tol = value;
				continue;
			}

			// Q-law gains and parameters
			else if (!strcmp(choice, "m_petro"))
			{
				probdata->m_petro = value;
				continue;
			}
			else if (!strcmp(choice, "n_petro"))
			{
				probdata->n_petro = value;
				continue;
			}
			else if (!strcmp(choice, "r_petro"))
			{
				probdata->r_petro = value;
				continue;
			}
			else if (!strcmp(choice, "b_petro"))
			{
				probdata->b_petro = value;
				continue;
			}
			else if (!strcmp(choice, "k_petro"))
			{
				probdata->k_petro = value;
				continue;
			}
			else if (!strcmp(choice, "c_petro"))
			{
				probdata->c_petro = value;
				continue;
			}
			else if (!strcmp(choice, "Wp"))
			{
				probdata->Wp = value;
				continue;
			}
			else if (!strcmp(choice, "Wc"))
			{
				probdata->Wc = value;
				continue;
			}
			else if (!strcmp(choice, "rpermin"))
			{
				probdata->rpermin = value;
				continue;
			}
			else if (!strcmp(choice, "rapomax"))
			{
				probdata->rapomax = value;
				continue;
			}
			else if (!strcmp(choice, "Wsma"))
			{
				probdata->Wsma = value;
				continue;
			}
			else if (!strcmp(choice, "We"))
			{
				probdata->We = value;
				continue;
			}
			else if (!strcmp(choice, "Winc"))
			{
				probdata->Winc = value;
				continue;
			}
			else if (!strcmp(choice, "Wape"))
			{
				probdata->Wape = value;
				continue;
			}
			else if (!strcmp(choice, "Wran"))
			{
				probdata->Wran = value;
				continue;
			}

			//effectivity settings
			else if (!strcmp(choice, "effectivity"))
			{
				probdata->effectivity = value;
				continue;
			}
			else if (!strcmp(choice, "abs_effectivity"))
			{
				probdata->abs_effectivity = value;
				continue;
			}
			else if (!strcmp(choice, "rel_effectivity"))
			{
				probdata->rel_effectivity = value;
				continue;
			}
			else if (!strcmp(choice, "eta_rel_tol"))
			{
				probdata->eta_rel_tol = value;
				continue;
			}
			else if (!strcmp(choice, "eta_abs_tol"))
			{
				probdata->eta_abs_tol = value;
				continue;
			}
			else if (!strcmp(choice, "effectivity_sample"))
			{
				probdata->effectivity_sample = value;
				continue;
			}

			//integrator settings
			else if (!strcmp(choice, "Lf"))
			{
				probdata->Lf = value;
				continue;
			}
			else if (!strcmp(choice, "stepsize"))
			{
				probdata->stepsize = value;
				continue;
			}
			else if (!strcmp(choice, "t0"))
			{
				probdata->t0 = value;
				continue;
			}

			// if keyword is not recognized
			else
			{
				std::cout << "Error on line " << linenumber << std::endl;
				std::cout << "Could not recognize keyword: " << choice << std::endl;
				return false;
			}

			
		}
		settings.close();
	}

	else  //file didn't open
	{
		std::cout << "File error! Couldn't find or open the file: " << filename << std::endl;
		return false;
	}// end if file open

	return true;
}