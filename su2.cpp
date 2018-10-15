#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include <functional>
#include <map>
#include "link.h"
#include "lattpoint.h"
using namespace std;

struct LattLink {
	Link U[4];
};

const int L = 10;		// grid width
const int L4 = L*L*L*L;	// number of lattice points
const int N = L4*4;		// number of link variables
const int thermSteps = 20;// total number of monte-carlo steps to take
const int MCSteps = 10;

map<LattPoint, LattLink> latt;	// link variables, identified by a point and index in 4D

////////////////////////// random number generators
random_device rd;
mt19937 mtunif(rd());
mt19937 mtgaus(rd());
uniform_real_distribution<double> unifDistr(0,1);
normal_distribution<double> gausDistr(0,1);

auto randUnif = bind(unifDistr, mtunif);
auto randGaus = bind(gausDistr, mtgaus);

////////////////////////// function declarations=
void initialize();
void HeatBathStep(const double beta);
void oneMonteCarloStepPerLink(const double beta);
double wilsonLoop(const int S);


int main(int argc, char* argv[]) {
	// open data file
	ofstream dataFile("su2.data");
	dataFile << "(beta, W(1x1)_mean, W(1x1)_std, W(2x2)_mean, W(2x2)_std, W(3x3)_mean, W(3x3)_std, "
		<< "W(4x4)_mean, W(4x4)_std, W(5x5)_mean, W(5x5)_std)" << '\n' << '\n';

	
	// run MC algorithm and calculate correlation length for beta values between
	// 0.1 and 3.0, with a step size of 0.1
	for(int i = 1; i <= 30; i++) {
		double beta = 0.1 * i;
		cout << " Calculating correlation length with beta = " << beta << '\n' << flush;
				
		initialize();
	
		// allow the system to come to thermal equilibrium at the given temperature
		cout << " Performing " << thermSteps
			<< " steps to thermalize the system " << flush;
		for(int l = 0; l < thermSteps; l++) {
			oneMonteCarloStepPerLink(beta);
			if((l+1) % (thermSteps/5) == 0)
				cout << "." << flush;
		}
	
		// perform production steps
		cout << " Done\n Performing production steps " << flush;
		double loop[5][MCSteps] = {0};
		double loopAv[5] = {0};
		for(int l = 0; l < MCSteps; l++) {
			oneMonteCarloStepPerLink(beta);
			
			for(int s = 0; s < 5; s++) {
				loop[s][l] = wilsonLoop(s+1);
				loopAv[s] += loop[s][l];
			}
			//if((l+1) % (MCSteps/20) == 0)
			cout << "." << flush;
		}
		cout << " Done\n\n" << flush;
		for(int s = 0; s < 5; s++)
			loopAv[s] /= MCSteps - 1;
	
		// calculate statistics of correlation length
		double loopStd[5] = {0};
		for(int s = 0; s < 5; s++) {
			for(int l = 0; l < MCSteps; l++) {
				double err = loopAv[s] - loop[s][l];
				loopStd[s] += err*err;
			}

			loopStd[s] /= MCSteps;
			loopStd[s] = sqrt(loopStd[s]);
		}
		
		// write data to file
		dataFile << "(" << beta;
		for(int s = 0; s < 5; s++)
			dataFile << ", " << loopAv[s] << ", " << loopStd[s];
		dataFile << ")" << '\n';
	}
	
	dataFile.close();
}


void initialize() {
	// randomly generate a starting configuration
	for(int i = 0; i < L4; i++) {
		// define Lattpoint
		int xi[4];

		xi[0] =  i        % 10;
		xi[1] = (i/10   ) % 10;
		xi[2] = (i/100  ) % 10;
		xi[3] = (i/1000 ) % 10;

		LattPoint Xi(xi, L);

		// prepare link variable per index
		LattLink Ui;

		for(int mu = 0; mu < 4; mu++) {
			// pick a point on the 3-sphere by choosing 4 normally distributed variables 
			double x0 = randGaus();
			double x1 = randGaus();
			double x2 = randGaus();
			double x3 = randGaus();

			double norm = sqrt(x0*x0 + x1*x1 + x2*x2 + x3*x3);

			double a[4];
			a[0] = x0/norm;
			a[1] = x1/norm;
			a[2] = x2/norm;
			a[3] = x3/norm;

			Ui.U[mu].setA(a);
		}

		latt[Xi] = Ui;
	}
}

void HeatBathStep(const double beta) {
	// choose a random point and index
	int xi[4];
	xi[0] = int(randUnif()*L);
	xi[1] = int(randUnif()*L);
	xi[2] = int(randUnif()*L);
	xi[3] = int(randUnif()*L);
	int mu = int(randUnif()*4);

	// create random lattpoint
	LattPoint Xi(xi, L);

	// create index mu direction
	int xmu[4];

	for(int i = 0; i < 4; i++) {
		if(i == mu) xmu[i] = 1;
		else xmu[i] = 0;
	}

	LattPoint Xmu(xmu, L);

	// Utilde is the sum of 6 products of 3 link variables which interact with the random link
	double UtildeA[4] = {0};
	Link Utilde(UtildeA);

	// calculate Utilde
	for(int nu = 0; nu < 4; nu++) {
		if(mu == nu) continue;

		// create index nu directions
		int xnu[4];

		for(int i = 0; i < 4; i++) {
			if(i == nu) xnu[i] = 1;
			else xnu[i] = 0;
		}

		LattPoint Xnu(xnu, L);

		Utilde += latt[Xi+Xmu].U[nu] * latt[Xi+Xnu].U[mu].inv() * latt[Xi].U[nu].inv() +
			latt[Xi+Xmu-Xnu].U[nu].inv() * latt[Xi-Xnu].U[mu].inv() * latt[Xi-Xnu].U[nu];
	}

	// calculate k = sqrt(det(\tilde{U}))
	double k = sqrt(Utilde.det());
	double w = exp(-2 * beta * k);

	double a[4];
	
	// choose a0 with P(a0) ~ sqrt(1 - a0^2) * exp(beta * k * a0)
	do {
		double xtrial = randUnif()*(1.0 - w) + w;
		a[0] = 1 + log(xtrial)/(beta * k);
	} while(sqrt(1-a[0]*a[0]) < randUnif());

	// choose \vec{a} randomly
	double r = sqrt(1-a[0]*a[0]);

	double a1 = randGaus();
	double a2 = randGaus();
	double a3 = randGaus();

	double norm = sqrt(a1*a1 + a2*a2 + a3*a3);

	a[1] = a1 * r / norm;
	a[2] = a2 * r / norm;
	a[3] = a3 * r / norm;

	Link UUbar(a);

	Link newU = UUbar * Utilde.scale(1/k).inv();

	latt[Xi].U[mu] = newU;
}

void oneMonteCarloStepPerLink(const double beta){
	for(int i = 0; i < 2 * N; i++){
		HeatBathStep(beta);
	}
}

double wilsonLoop(const int S) {
	double loopSum = 0;

	// loop over all points in the lattice
	for(int i = 0; i < L4; i++) {
		// define Lattpoint
		int xi[4];

		xi[0] =  i        % 10;
		xi[1] = (i/10   ) % 10;
		xi[2] = (i/100  ) % 10;
		xi[3] = (i/1000 ) % 10;

		LattPoint Xi(xi, L);

		// 6 unique loops in the positive direction at any point
		for(int mu = 0; mu < 4; mu++) {
			for(int nu = 0; nu < mu; nu++) {
				// create index directions
				int xmu[4];
				int xnu[4];

				for(int i = 0; i < 4; i++) {
					if(i == mu) xmu[i] = 1;
					else xmu[i] = 0;

					if(i == nu) xnu[i] = 1;
					else xnu[i] = 0;
				}

				LattPoint Xmu(xmu, L);
				LattPoint Xnu(xnu, L);

				Link loop;

				// build our loop from the 4*S links on the perimeter
				for(int s = 0; s < S; s++)
					loop *= latt[Xi+Xmu.scale(s)].U[mu];
				for(int s = 0; s < S; s++)
					loop *= latt[Xi+Xmu.scale(S)+Xnu.scale(s)].U[nu];
				for(int s = S-1; s >= 0; s--)
					loop *= latt[Xi+Xmu.scale(s)+Xnu.scale(S)].U[mu].inv();
				for(int s = S-1; s >= 0; s--)
					loop *= latt[Xi+Xnu.scale(s)].U[nu].inv();
				loopSum += loop.tr()/2;
			}
		}
	}

	return loopSum / (L4 * 6);
}