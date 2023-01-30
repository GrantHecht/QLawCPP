/*
Q law originally developed by Anastassios E. Petropoulos (JPL)

c++ implementation by Donald Ellison (UIUC)
March 5th 2014

Transferred to Jacob Englander and Matthew Vavrina (GSFC) March 16th 2014
*/

#include "qlaw.h"

int main(int argc, char* argv[])
{
	std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n";
	std::cout << "%" << "\n";
	std::cout << "% Q law low thrust transfer tool" << "\n";
	std::cout << "%" << "\n";
	std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n \n";

	// data structure
	options probdata;

	// read options file
	std::string options_file;
	bool options_read = false;

	if (argc == 1)
		options_file = "options.txt";
	else if (argc == 2)
		options_file.assign(argv[1]);

	options_read = read_options_file(options_file, &probdata);

	if (options_read)
	{
		std::cout << "options file successfully read \n";
		std::cout << std::endl;
	}

	else
	{
		std::cout << "options file not read.....press ENTER to abort \n";
		getchar();
		return 0;
	}


	// setup output data file
	std::filebuf fb;
	char request;
	std::cout << "You are about to delete the contents of the output file? Are you sure you want to proceed? y/n: ";
	std::cin >> request;
	std::cin.get();

	
	if (request == 'y')
		fb.open("qlaw_output.txt", std::ios::out);
	else
	{
		std::cout << "Hit ENTER to exit...";
		getchar();
		return 0;
	}

	std::ostream outfile(&fb);
	std::cout << std::endl;
	
	// Set up output file headers
	outfile.width(21); outfile << std::left << "p";
	outfile.width(21); outfile << std::left << "e";
	outfile.width(21); outfile << std::left << "inc";
	outfile.width(21); outfile << std::left << "ape";
	outfile.width(21); outfile << std::left << "ran";
	outfile.width(21); outfile << std::left << "mass";
	outfile.width(21); outfile << std::left << "time";
	outfile.width(21); outfile << std::left << "true longitude";
	outfile << "\n \n";



	probdata.APPLE_PI = 3.14159265358979323;

	//%%%%%%%%%%%%%%%%%%%%%%%
	//SCALING OF USER OPTIONS
	//%%%%%%%%%%%%%%%%%%%%%%%
	probdata.grav = 9.80665*1.0e-3*probdata.TU*probdata.TU / probdata.DU; // DU/TU^2

	// Spacecraft parameters
	probdata.thrust = probdata.thrust / 1000.0*probdata.TU*probdata.TU / probdata.DU; //thrust in kg*DU/TU^2
	probdata.Isp = probdata.Isp / probdata.TU; //TU's

	// initial coe's
	double sma = probdata.sma0 / probdata.DU; //DU's
	double e = probdata.e0;
	double inc = probdata.inc0*(probdata.APPLE_PI / 180.0); //radians
	double ape = probdata.ape0*(probdata.APPLE_PI / 180.0); //radians
	double ran = probdata.ran0*(probdata.APPLE_PI / 180.0); //radians
	double tru = probdata.tru0*(probdata.APPLE_PI / 180.0); //radians
	
	double L0 = ape + ran + tru; //initial true longitude
	double p = sma*(1.0 - e*e);
	
	
    // target coe's
	probdata.sma_t = probdata.sma_t / probdata.DU; //DU's
	//probdata.e_t = 0.01;
	probdata.inc_t = probdata.inc_t*(probdata.APPLE_PI / 180); //radians
	probdata.ape_t = probdata.ape_t*(probdata.APPLE_PI / 180);  //radians
	probdata.ran_t = probdata.ran_t*(probdata.APPLE_PI / 180); //radians


	// Convergence tolerances
	probdata.sma_tol = probdata.sma_tol / probdata.DU; //DU's
	//probdata.e_tol = 0.001;
	probdata.inc_tol = probdata.inc_tol*probdata.APPLE_PI / 180.0; //radians
	probdata.ape_tol = probdata.ape_tol*probdata.APPLE_PI / 180.0; //radians
	probdata.ran_tol = probdata.ran_tol*probdata.APPLE_PI / 180.0; //radians
	
	

	// Q-law gains and parameters require tuning, suggest NSGAII
	probdata.rpermin = probdata.rpermin / probdata.DU;
	probdata.rapomax = probdata.rapomax / probdata.DU;
	
	// scale effectivity sample size
	probdata.effectivity_sample = probdata.effectivity_sample*probdata.APPLE_PI / 180.0;
	


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// PERFORM FEEDBACK CONTROL SIMULATION
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	// The final value of the true longitude is a bit arbitrary
	// Set it to a "large enough" value and let the convergence check stop the simulation
	double L;
	double L_bin = 0.0;
	probdata.nseg = (probdata.Lf - L0) / probdata.stepsize;
	probdata.ns = 7;

	// What wall-clock time do you want to start integrating from?
	probdata.t0 = probdata.t0 / probdata.TU;

	double dX[7];
	double error = 0.0;
	int throttle_lock = 0;
	int coast = 0;

	double X_next[7];
	double X_current[7] = { p, e, inc, ape, ran, probdata.initial_mass, probdata.t0 };
	double Xf[7];
	double sma_current;
	double L_current = L0;
	double percentComplete = 0.0;
	double accumulate = 1.0;
	bool converged;

	//std::vector< std::vector<double> > Xbucket; //place to store COE history
	//std::vector< std::vector<double> > Lbucket; //place to store independent variable history
	std::vector< std::vector< double> > writebin(ceil(probdata.nseg / 100),std::vector<double>(8, 0));
	//std::vector< std::vector< float > > item(2, std::vector<float>(2, 0));
	//writebin.resize(ceil(probdata.nseg / 100));

	// Write initial COE's to file
	outfile.precision(10);
	outfile.width(20); outfile << std::left << p << " ";
	outfile.width(20); outfile << std::left << e << " ";
	outfile.width(20); outfile << std::left << inc << " ";
	outfile.width(20); outfile << std::left << ape << " ";
	outfile.width(20); outfile << std::left << ran << " ";
	outfile.width(20); outfile << std::left << probdata.initial_mass << " ";
	outfile.width(20); outfile << std::left << probdata.t0 << " ";
	outfile.width(20); outfile << std::left << L0 << " ";
	outfile << "\n";

	//Adaptive step sizing doesn't seem to work very well
	//L = adaptive_step_int(X_current, Xf, f1, L0, probdataptr);

	// Fixed step integration routine
	std::cout << "Simulation Progress: \n";
	int writerow = 0;
	for (size_t j = 0; j < probdata.nseg; ++j)
	{
		// If we run out of mass, the integration will blow up
		if (X_current[5] <= 0.0)
		{
			std::cout << "Out of mass" << std::endl;
			break;
		}

		// Take an integration step
		try
		{
			gaussvar(X_current, dX, L_current, probdata.stepsize, L_bin, throttle_lock, coast, &probdata);
			rk8713M(X_current, X_next, dX, probdata.stepsize, L_current, probdata.ns, error, L_bin, throttle_lock, coast, &probdata);
		}
		catch (int e)
		{
			break;
		}

		// Keep inclination and eccentricity away from 0
		if (abs(X_next[2]) < 1.0e-4)
			X_next[2] > 0.0 ? 1.0e-4 : -1.0e-4;
		if (X_next[1] < 1.0e-4)
			X_next[1] = 1.0e-4;

		for (size_t m = 0; m < probdata.ns; ++m)
			X_current[m] = X_next[m];

		
		// advance independent integration variable
		L_current = L_current + probdata.stepsize;
		L_bin = L_bin + probdata.stepsize;

		// Store the state vector and current true longitude every so often for writing to file later
		if (j % 100 == 0)
		{
			writebin[writerow][0] = X_current[0];
			writebin[writerow][1] = X_current[1];
			writebin[writerow][2] = X_current[2];
			writebin[writerow][3] = X_current[3];
			writebin[writerow][4] = X_current[4];
			writebin[writerow][5] = X_current[5];
			writebin[writerow][6] = X_current[6];
			writebin[writerow][7] = L_current;
			++writerow;
		}

		// Check for convergence to target COE's
		// We don't know how long the transfer will take so this will stop the integration if we converge to the target orbit

		converged = check_convergence(X_current, &probdata);
		if (converged)
		{
			std::cout << std::endl;
			std::cout << "CONVERGED" << std::endl << std::endl;
			break;
		}
		
		

		// Calculate progress of integration
		// This will not reach 100% if the orbit converges before the integration is complete
		percentComplete = L_current / probdata.Lf*100.0;
		if (percentComplete >= accumulate)
		{
			std::cout << int(floor(percentComplete)) <<"% " <<"of trajectory computed" <<std::endl;
			accumulate = accumulate + 1.0;
		}

		
	} // end integration loop

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// WRITE THE SIMULATION HISTORY TO FILE FOR VISUALIZATION AND PLOTTING
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	outfile.precision(10);
	int m;
	for (size_t k = 0; k < writerow; ++k)
	{
		for (m = 0; m < probdata.ns; ++m)
		{
			if (writebin[k][m] == writebin[k][m])
			{
				outfile.width(20);
				outfile << std::left << writebin[k][m] << " ";
			}
		}
		outfile.width(20); outfile << std::left << writebin[k][m];
		outfile << '\n';
	}

	// Close the output file
	fb.close();

	// Final COE values
	double p_f = X_current[0];
	double e_f = X_current[1];
	double inc_f = X_current[2];
	double ape_f = X_current[3];
	double ran_f = X_current[4];
	double mass_f = X_current[5];
	double t_f = X_current[6];

	double tru_f = L_current - ran_f - ape_f;
	double sma_f = p_f / (1.0 - e_f*e_f);

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// REPORT SOME INFORMATION BACK TO THE USER
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	std::cout << "SIMULATION RESULTS: \n \n";
	std::cout << "final semimajor axis: "        << sma_f*probdata.DU << " km" << std::endl;
	std::cout << "final eccentricity: "               << e_f   << std::endl;
	std::cout << "final inclination: "          << inc_f*180.0 / probdata.APPLE_PI << " deg" << std::endl;
	std::cout << "final argument of periapse: " << ape_f*180.0 / probdata.APPLE_PI << " deg" << std::endl;
	std::cout << "final RAN: "                  << ran_f*180.0 / probdata.APPLE_PI << " deg" <<std::endl;
	std::cout << "Final mass: "                 << mass_f << " kg" << std::endl;
	std::cout << "Flight time: "                << t_f * probdata.TU / 86400.0 << " days" << std::endl << std::endl;
	
	if (probdata.sma_switch)
		std::cout << "semimajor axis targeting error: "        << fabs(probdata.sma_t - sma_f)*probdata.DU << " km" << std::endl;
	else
		std::cout << "semimajor axis targeting error: "        << "NA"                                     << std::endl;

	if (probdata.e_switch)
		std::cout << "eccentricity targeting error: " << fabs(probdata.e_t - e_f) << std::endl;
	else
		std::cout << "eccentricity targeting error: " << "NA"                     << std::endl;

	if (probdata.inc_switch)
		std::cout << "inclination targeting error: " << fabs(probdata.inc_t - inc_f)*180.0 / probdata.APPLE_PI << " deg" << std::endl;
	else
		std::cout << "inclination targeting error: " << "NA"                                                   << std::endl;

	if (probdata.ape_switch)
		std::cout << "argument of periapsis targeting error: " << fabs(probdata.ape_t - ape_f)*180.0 / probdata.APPLE_PI << " deg" << std::endl;
	else
		std::cout << "argument of periapsis targeting error: "                                                           << "NA" << std::endl;

	if (probdata.ran_switch)
		std::cout << "RAN targeting error: " << fabs(probdata.ran_t - ran_f)*180.0 / probdata.APPLE_PI << " deg" << std::endl;
	else
		std::cout << "RAN targeting error: "                                                           << "NA" << std::endl;

	getchar();

} //end main