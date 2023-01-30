// Coded by Donald Ellison
// March 5th 2014

#ifndef QLAWHEADER
#define QLAWHEADER

#include <vector>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <math.h>

// Data structure
struct options{
	double APPLE_PI;
	double DU;
	double TU;
	double mu;
	double grav;
	double thrust;
	double Isp;
	double m_petro;
	double n_petro;
	double r_petro;
	double b_petro;
	double k_petro;
	double c_petro;
	double Wp;
	double Wc;
	double rpermin;
	double rapomax;
	double Wsma;
	double We;
	double Winc;
	double Wape;
	double Wran;
	double nseg;
	double stepsize;
	int ns;
	double resumeError;
	double precisionError;
	double precisionTarget;
	int sma_switch;
	int e_switch;
	int inc_switch;
	int ape_switch;
	int ran_switch;
	double sma0;
	double e0;
	double inc0;
	double ape0;
	double ran0;
	double tru0;
	double initial_mass;
	double sma_t;
	double e_t;
	double inc_t;
	double ape_t;
	double ran_t;
	double sma_tol;
	double e_tol;
	double inc_tol;
	double ape_tol;
	double ran_tol;
	int effectivity;
	double abs_effectivity;
	double rel_effectivity;
	double eta_rel_tol;
	double eta_abs_tol;
	double effectivity_sample;
	double Lf;
	double t0;
};


// Function Prototypes
double adaptive_step_int(double x_left[], double *x_right, double f1[], double t0, void * probdataptr);
void rk8713M(double x_left[], double x_right[], double f1[], double & h, double t, int & ns, double & error, double & L_bin, int & throttle_lock, int & coast, options * probdata);
void gaussvar(double X[], double dX[], double L, double & stepsize, double & L_bin, int & throttle_lock, int & coast, options * probdata);
void lyapunov(double u[], double & sma, double & e, double & inc, double & ape, double & ran, double & tru, double & f, options * probdata);
void lyapunov_rapomax(double u[], double & sma, double & e, double & inc, double & ape, double & ran, double & tru, double & f, options * probdata);
void effectivity(double & eta_abs, double & eta_rel, double & sma, double & e, double & inc, double & ape, double & ran, double & tru, double & f, double L_bin, options * probdata);
double minQdot(double & sma, double & e, double & inc, double & ape, double & ran, double & tru, double & f, options * probdata);
bool read_options_file(std::string filename, options * probdata);
bool check_convergence(double X_current[], options * probdata);

#endif