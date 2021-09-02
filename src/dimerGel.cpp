/*

	Jack Treado
	08/2021
	dimerGel.cpp

	Generate gel of dimers via athermal extension

	compile: g++ -O3 src/dimerGel.cpp -o dgel.o
	./dgel.o 32 0.1 1e-4 1e-4 0 0.05 1e-8 1 pos.test

*/

// preprocessor directives
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>

# define NDIM 3
# define NNN 13

// namespace
using namespace std;



// GLOBAL CONSTANTS
const double PI 			= 4*atan(1);
const int w 				= 6;
const int wnum 				= 25;
const int pnum 				= 12;

const double phi0 			= 1.1;
const double phimin 		= 0.3;
const double timeStepMag 	= 0.01;
const double dr 			= 0.05;
const double Umin 			= 1e-16;
const double l0min 			= 1e-4;
const double kl 			= 0.5;

const double alpha0      	= 0.2;
const double finc        	= 1.1;
const double fdec        	= 0.5;
const double falpha      	= 0.99;

const double phiskip 		= 0.01;
const int NSKIP 			= 1e3;
const int NMIN        		= 100;
const int NNEGMAX     		= 2000;
const int NDELAY      		= 1000;
const int itmax       		= 1e7;  


// function prototypes
void printPos(ofstream &posout, vector<double> &pos, vector<double> &radii, vector<double> &L, vector<double> &S, vector<int> &z, int N);
int cmindex(int i, int j, int N);
void addContacts(int i, int j, vector<bool> &cij, vector<int> &z, int N);
void cellLinkedList(vector<int> &sc, vector<double> &lc, vector< vector<int> > &nn, int &NCELLS, vector<double> &L, double boxsize);

// MAIN
int main(int argc, char const *argv[])
{
	// variables used throughout main
	int i, j, d, ind;
	double sij, rij, dp, ftmp, atmp;
	vector<double> vij(NDIM,0.0);

	// parameters to be read in 
	int N, seed;
	double l2, dl0, dphi, dg, del, Ftol;

	// read in parameters from command line input
	string N_str 		= argv[1];
	string dl0_str 		= argv[2];
	string dphi_str 	= argv[3];
	string dg_str 		= argv[4]; 	// NOTE: in "units" of delta_phi, so will scale with dphi parameter
	string del_str 		= argv[5];
	string l2_str 		= argv[6];
	string Ftol_str 	= argv[7];
	string seed_str 	= argv[8];
	string posFile 		= argv[9];

	stringstream Nss(N_str);
	stringstream dl0ss(dl0_str);
	stringstream dphiss(dphi_str);
	stringstream dgss(dg_str);
	stringstream delss(del_str);
	stringstream l2ss(l2_str);
	stringstream Ftolss(Ftol_str);
	stringstream seedss(seed_str);

	Nss >> N;
	dl0ss >> dl0;
	dphiss >> dphi;
	dgss >> dg;
	delss >> del;
	l2ss >> l2;
	Ftolss >> Ftol;
	seedss >> seed;

	// number of monomers
	int NMTOT = 2*N;

	// open xyz file
	ofstream posout;
	posout.open(posFile.c_str());
	if (!posout.is_open()){
		cout << "	** ERROR: xyz file " << posFile << " could not be opened, ending." << endl;
		return 1;
	}


	// check del parameter
	if (del < 0.0 || del > 1.0){
		cout << "** ERROR: del = " << del << ", should be within 0 to 1. Ending here. " << endl;
		return 1;
	}


	// output opening statement to console
	cout << "=======================================================" << endl << endl;
	cout << "		dimerGel.cpp 									" << endl;
	cout << "		Athermal gelation of densely-packed dimers 		" << endl;
	cout << "		N 			= " << N << "						" << endl;
	cout << "		dr 			= " << dr << " 						" << endl;
	cout << "		dphi 		= " << dphi << " 					" << endl;
	cout << "		dg 			= " << dg << "						" << endl;
	cout << "		del 		= " << del << "						" << endl;
	cout << "		l2 			= " << l2 << " 						" << endl;
	cout << "		Ftol 		= " << Ftol << " 					" << endl;
	cout << "		seed 		= " << seed << "					" << endl;
	cout << "		pos file 	= " << posFile << "					" << endl;
	cout << endl;
	cout << "=======================================================" << endl << endl;


	// seed random number generator
	srand48(seed);


	/* * * * * * * * * * * * * * * * * * 

				INITIALIZE

					PARTICLES

	 * * * * * * * * * * * * * * * * * */

	// initialization variables
	double r1, r2, g1, g2, vp, vcap, rmax, L0;

	// particle radii and positions
	vector<double> radii(NMTOT,0.0);
	vector<double> l0(N,0.0);
	vector<double> mass(NMTOT,0.0);
	vector<double> pos(NDIM*NMTOT,0.0);

	// gaussian-distributed polydispersity
	vp = 0.0;
	rmax = 0.0;
	for (i=0; i<N; i++){
		// generate random numbers
		r1 = drand48();
		r2 = drand48();

		// calculate gaussian random variable using Box-Muller transform
		g1 = sqrt(-2.0*log(r1))*cos(2*PI*r2);
		g2 = sqrt(-2.0*log(r1))*sin(2*PI*r2);

		// get radii from g1
		r1 = g1*dr + 1.0;
		radii[2*i] = r1;
		radii[2*i + 1] = r1;

		// get l0 from g2
		l0[i] = (g2*dl0 + 2.0)*r1;
		if (l0[i] < 0)
			l0[i] = l0min;

		// compute mass
		mass[2*i] = (4.0/3.0)*PI*pow(r1,3.0);
		mass[2*i + 1] = mass[2*i];

		// add to particle volume
		vp += (8.0/3.0)*PI*pow(r1,3.0);
		if (l0[i] < 2.0*r1)
			vp -= (PI/12.0)*(4.0*r1 + l0[i])*pow(2.0*r1 - l0[i],2.0);

		// determine max radius
		if (radii[i] > rmax)
			rmax = radii[i];
	}


	// box lengths (initially a cube)
	L0 = pow(vp/(phi0),1.0/3.0);
	vector<double> L(NDIM,L0);
	double Lmin;
	
	// initialize particle positions randomly throughout box
	for (i=0; i<N; i++){
		// first monomer = random
		pos[2*NDIM*i] = L[0]*drand48();
		pos[2*NDIM*i + 1] = L[1]*drand48();
		pos[2*NDIM*i + 2] = L[2]*drand48();

		// second monomer = along x
		pos[2*NDIM*i + 3] = pos[2*NDIM*i] + l0[i];
		pos[2*NDIM*i + 4] = pos[2*NDIM*i + 1];
		pos[2*NDIM*i + 5] = pos[2*NDIM*i + 2];
	}


	// determine fundamental MD time unit
	double dtMD, dt0, dt;

	dtMD 	= 1.0;
	dt0 	= timeStepMag*dtMD;
	dt 		= dt0;




	/* * * * * * * * * * * * * * * * * * 

			INITIALIZE

				CELL LINKED LIST

	 * * * * * * * * * * * * * * * * * */

	// min box size
	double boxsize = 2.1*rmax;
	bool resetCLL;
	int NCELLS = 1;

	// number of boxes and lengths in each direction
	vector<int> sc(NDIM,0);
	vector<double> lc(NDIM,0.0);
	vector< vector<int> > nn;

	// use function to initialize cell linked list and box neighbors
	cellLinkedList(sc,lc,nn,NCELLS,L,boxsize);

	// linked-list variables
	vector<int> head(NCELLS,0);
	vector<int> last(NCELLS,0);
	vector<int> list(NMTOT+1,0);





	/* * * * * * * * * * * * * * * * * * 

			INITIAL FIRE

					MINIMZATION

	 * * * * * * * * * * * * * * * * * */


	// initialize velocity and force information
	vector<double> v(NDIM*NMTOT,0.0);
	vector<double> a(NDIM*NMTOT,0.0);
	vector<double> F(NDIM*NMTOT,0.0);

	// FIRE VARIABLES
	double P  		= 0;	
	double fnorm 	= 0;
	double vnorm 	= 0;
	double alpha   	= alpha0;

	double dtmax   	= 10*dt0;
	double dtmin   	= 1e-2*dt0;
	dt 				= dt0;

	int npPos      	= 0;
	int npNeg      	= 0;
	int npPMin      = 0;

	int fireit    	= 0;
	double fcheck  	= 10*Ftol;
	double U 		= 0;

	// cell linked-list variables
	int cellid, ci, cj, pi, pj, sctmp;

	// dimer variables
	double lx, ly, lz, l, ux, uy, uz, flx, fly, flz;

	// loop until force relaxes
	while ((fcheck > Ftol || npPMin < NMIN) && fireit < itmax){

		// VELOCITY-VERLET UPDATE 1: POSITIONS
		for (i=0; i<NMTOT; i++){
			for(d=0; d<NDIM; d++){
				// dof index
				ind = NDIM*i + d;

				// update positions
				pos[ind] += dt*v[ind] + 0.5*dt*dt*a[ind];

				// put back into box
				if (pos[ind] > L[d])
					pos[ind] -= L[d];
				else if (pos[ind] < 0)
					pos[ind] += L[d];

				// reset forces for force computation
				F[ind] = 0.0;
			}

			// reset linked list
			list[i] = 0;
		}
		list[NMTOT] = 0;


		// reset linked list head
		for (i=0; i<NCELLS; i++){
			head[i] = 0;
			last[i] = 0;
		}

		// sort particles into linked list
		for (i=0; i<NMTOT; i++){
			// 1. get cell id of current particle position
			cellid = 0;
			sctmp = 1;
			for (d=0; d<NDIM; d++){
				// add d index to 1d list
				cellid += floor(pos[NDIM*i + d]/lc[d])*sctmp;

				// increment dimensional factor
				sctmp *= sc[d];
			}

			// 2. add to head list or link within list
			// NOTE: particle ids are labelled starting from 1, pting to 0 means end of linked list
			if (head[cellid] == 0){
				head[cellid] = i + 1;
				last[cellid] = i + 1;
			}
			else{
				list[last[cellid]] = i + 1;
				last[cellid] = i + 1;
			}
		}


		// force computation using cell linked-lists
		U = 0;
		for (ci=0; ci<NCELLS; ci++){

			// get start of list of particles
			pi = head[ci];

			// loop over linked list
			while(pi > 0){
				// real particle index
				i = pi - 1;

				// next in list
				pj = list[pi];

				// loop down neighbors of pi - 1 in same cell
				while(pj > 0){
					// real index of pj
					j = pj - 1;	

					// check if dimer pair
					if ((i % 2 == 0 && j == i+1) || (j % 2 == 0 && i == j+1) ){
						// update pj and continue
						pj = list[pj];
						continue;
					}

					// contact distance
					sij = radii[i] + radii[j];

					// particle distance
					rij = 0.0;
					for (d=0; d<NDIM; d++){
						// distance across box
						dp = pos[NDIM*j + d] - pos[NDIM*i + d];

						// image distance
						dp -= L[d]*round(dp/L[d]);

						// add to total distance
						rij += dp*dp;

						// save vector component
						vij[d] = dp;
					}

					// check for overlap
					if (rij < sij*sij){
						// true distance
						rij = sqrt(rij);

						// force scale
						ftmp = (1.0 - (rij/sij))/sij;

						// add to forces
						for (d=0; d<NDIM; d++){
							F[NDIM*i + d] -= ftmp*(vij[d]/rij);
							F[NDIM*j + d] += ftmp*(vij[d]/rij);
						}

						// increase potential energy
						U += 0.5*pow(1 - (rij/sij),2.0);
					}

					// update pj
					pj = list[pj];
				}

				// test overlaps with forward neighboring cells
				for (cj=0; cj<NNN; cj++){
					// get first particle in neighboring cell
					pj = head[nn[ci][cj]];

					// loop over cells in neighboring cells
					while(pj > 0){
						// real index of pj
						j = pj - 1;	

						// check if dimer pair
						if ((i % 2 == 0 && j == i+1) || (j % 2 == 0 && i == j+1) ){
							// update pj and continue
							pj = list[pj];
							continue;
						}

						// contact distance
						sij = radii[i] + radii[j];

						// particle distance
						rij = 0.0;
						for (d=0; d<NDIM; d++){
							// distance across box
							dp = pos[NDIM*j + d] - pos[NDIM*i + d];

							// image distance
							dp -= L[d]*round(dp/L[d]);

							// add to total distance
							rij += dp*dp;

							// save vector component
							vij[d] = dp;
						}

						// check for overlap
						if (rij < sij*sij){
							// true distance
							rij = sqrt(rij);

							// force scale
							ftmp = (1.0 - (rij/sij))/sij;

							// add to forces
							for (d=0; d<NDIM; d++){
								F[NDIM*i + d] -= ftmp*(vij[d]/rij);
								F[NDIM*j + d] += ftmp*(vij[d]/rij);
							}

							// increase potential energy
							U += 0.5*pow(1 - (rij/sij),2.0);
						}

						// update pj
						pj = list[pj];
					}
				}

				// update pi
				pi = list[pi];
			}
		}


		// Add dimer bond force
		for (i=0; i<N; i++){
			// distance between pair images
			lx = pos[2*NDIM*i + 3] - pos[2*NDIM*i];
			lx -= L[0]*round(lx/L[0]);

			ly = pos[2*NDIM*i + 4] - pos[2*NDIM*i + 1];
			ly -= L[1]*round(ly/L[1]);

			lz = pos[2*NDIM*i + 5] - pos[2*NDIM*i + 2];
			lz -= L[2]*round(lz/L[2]);

			// bond length
			l = lx*lx + ly*ly + lz*lz;

			// unit vector
			ux = lx/l;
			uy = ly/l;
			uz = lz/l;

			// force on 1 due to 2
			ftmp = kl*(l - l0[i]);
			flx = ftmp*ux;
			fly = ftmp*uy;
			flz = ftmp*uz;

			// bonded force on monomer 1
			F[2*NDIM*i] += flx;
			F[2*NDIM*i + 1] += fly;
			F[2*NDIM*i + 2] += flz;

			// bonded force on monomer 2
			F[2*NDIM*i + 3] -= flx;
			F[2*NDIM*i + 4] -= fly;
			F[2*NDIM*i + 5] -= flz;

			// add to potential energy
			U += 0.5*kl*pow(l - l0[i],2.0);
		}


		// VELOCITY-VERLET UPDATE 2: VELOCITIES AND ACCELERATION
		for (i=0; i<NMTOT; i++){
			for(d=0; d<NDIM; d++){
				// dof index
				ind = NDIM*i + d;

				// new acceleration
				atmp = F[ind]/mass[i];

				// update velocities
				v[ind] += 0.5*(atmp + a[ind])*dt;

				// update acceleration
				a[ind] = atmp;
			}
		}

		// compute fnorm, vnorm and P
		fnorm = 0.0;
		vnorm = 0.0;
		P = 0.0;
		for (i=0; i<NMTOT; i++){
			for (d=0; d<NDIM; d++){
				ind = NDIM*i + d;
				fnorm 	+= F[ind]*F[ind];
				vnorm 	+= v[ind]*v[ind];
				P 		+= v[ind]*F[ind];
			}
		}
		fnorm = sqrt(fnorm);
		vnorm = sqrt(vnorm);

		// update fcheck based on fnorm (= force per degree of freedom)
		fcheck = fnorm/(NDIM*NMTOT);

		// update npPMin
		if (fcheck < Ftol)
			npPMin++;
		else
			npPMin = 0;

		// print to console
		if (fireit % NSKIP == 0){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << "		I N I T I A L 				" << endl;
			cout << " 	F I R E 						" << endl;
			cout << "		M I N I M I Z A T I O N 	" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** fireit = " << fireit << endl;
			cout << "	** fcheck = " << fcheck << endl;
			cout << "	** vnorm = " << vnorm << endl;
			cout << "	** dt = " << dt << endl;
			cout << "	** P = " << P << endl;
			cout << "	** Pdir = " << P/(fnorm*vnorm) << endl;
			cout << "	** alpha = " << alpha << endl;
			cout << "	** U = " << U << endl << endl;
		}

		// Step 1. adjust simulation based on net motion of degrees of freedom
		if (P > 0){
			// increase positive counter
			npPos++;

			// reset negative counter
			npNeg = 0;

			// alter simulation if enough positive steps have been taken
			if (npPos > NMIN){
				// change time step
				if (dt*finc < dtmax)
					dt *= finc;

				// decrease alpha
				alpha *= falpha;
			}
		}
		else{
			// reset positive counter
			npPos = 0;

			// increase negative counter
			npNeg++;

			// check if simulation is stuck
			if (npNeg > NNEGMAX){
				cout << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
				return 1;
			}

			// take half step backwards, reset velocities
			for (i=0; i<NMTOT; i++){
				for (d=0; d<NDIM; d++){
					// dof index
					ind = NDIM*i + d;

					// take half step backwards
					pos[ind] -= 0.5*dt*v[ind];

					// reset velocities
					v[ind] = 0;
				}
			}

			// decrease time step if past initial delay
			if (fireit > NDELAY){
				// decrease time step 
				if (dt*fdec > dtmin)
					dt *= fdec;

				// reset alpha
				alpha = alpha0;
			}
		}


		// update velocities (s.d. vs inertial dynamics) only if forces are acting
		if (fnorm > 0){
			for (i=0; i<NMTOT; i++){
				for (d=0; d<NDIM; d++){
					// dof index
					ind = NDIM*i + d;

					// update velocity
					v[ind] = (1 - alpha)*v[ind] + alpha*(vnorm/fnorm)*F[ind];
				}
			}
		}

		// update iterator
		fireit++;
	}
	// check if FIRE converged
	if (fireit == itmax){
		cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
		return 1;
	}
	else{
		cout << endl << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===========================================" << endl;
		cout << " 	F I R E 						" << endl;
		cout << "		M I N I M I Z A T I O N 	" << endl;
		cout << "	 C O N V E R G E D! 			" << endl;
		cout << "===========================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	** fireit = " << fireit << endl;
		cout << "	** fcheck = " << fcheck << endl;
		cout << "	** vnorm = " << vnorm << endl;
		cout << "	** dt = " << dt << endl;
		cout << "	** P = " << P << endl;
		cout << "	** Pdir = " << P/(fnorm*vnorm) << endl;
		cout << "	** alpha = " << alpha << endl;
		cout << "	** U = " << U << endl << endl;
	}


	/* * * * * * * * * * * * * * * * * * 

			DECOMPRESSION WITH

				GELATION & MINIMZATION

	 * * * * * * * * * * * * * * * * * */

	// initialize contact network
	int NPW = NMTOT*(NMTOT-1)/2;
	vector<bool> cij(NPW,0);
	vector<int> z(NMTOT,0);


	// initialize stress tensor
	vector<double> S(6,0.0);

	// instantaneous phi
	double phi = 0.0;
	vp = 0.0;
	for (i=0; i<N; i++){
		r1 = radii[2*i];
		vp += (8.0/3.0)*PI*pow(r1,3.0);
		if (l0[i] < 2.0*r1)
			vp -= (PI/12.0)*(4.0*r1 + l0[i])*pow(2.0*r1 - l0[i],2.0);
	}
	phi = vp/(L[0]*L[1]*L[2]);

	// decompression iterator
	int it = 0;
	double rscale;

	// attraction parameters
	double p1, u1, u2, h, lij, fdirtmp;

	// strain parameters
	double dgx, dgy, dgz;
	dg *= dphi;
	dgx = 1.0 + dg*(1.0 - del);
	dgy = 1.0 - (dg/(1.0 + dg))*(1 - del);
	dgz = 1.0 + dg*del;

	// determine frames to skip between xyz output
	int NSTEPS = round((phi - phimin)/dphi);
	int NXYZSTEPS = round((phi - phimin)/phiskip);
	int NPHISKIP = NSTEPS/NXYZSTEPS;

	// loop until phimin found
	while (phi > phimin && it < itmax){
		// reset FIRE variables that make have changed
		fireit 		= 0;
		fcheck 		= 10*Ftol;
		npPos       = 0;
		npNeg       = 0;
    	npPMin      = 0;
    	dt 			= dt0;
    	alpha 		= alpha0;

    	// loop over FIRE minimzation until force relaxes
		while ((fcheck > Ftol || npPMin < NMIN) && fireit < itmax){

			// VELOCITY-VERLET UPDATE 1: POSITIONS
			for (i=0; i<NMTOT; i++){
				for(d=0; d<NDIM; d++){
					// dof index
					ind = NDIM*i + d;

					// update positions
					pos[ind] += dt*v[ind] + 0.5*dt*dt*a[ind];

					// put back into box
					if (pos[ind] > L[d])
						pos[ind] -= L[d];
					else if (pos[ind] < 0)
						pos[ind] += L[d];

					// reset forces for force computation
					F[ind] = 0.0;
				}

				// reset contacts
				z[i] = 0;

				// reset linked list
				list[i] = 0;
			}
			list[NMTOT] = 0;


			// reset linked list head
			for (i=0; i<NCELLS; i++){
				head[i] = 0;
				last[i] = 0;
			}

			// sort particles into linked list
			for (i=0; i<NMTOT; i++){
				// 1. get cell id of current particle position
				cellid = 0;
				sctmp = 1;
				for (d=0; d<NDIM; d++){
					// add d index to 1d list
					cellid += floor(pos[NDIM*i + d]/lc[d])*sctmp;

					// increment dimensional factor
					sctmp *= sc[d];
				}

				// 2. add to head list or link within list
				// NOTE: particle ids are labelled starting from 1, pting to 0 means end of linked list
				if (head[cellid] == 0){
					head[cellid] = i + 1;
					last[cellid] = i + 1;
				}
				else{
					list[last[cellid]] = i + 1;
					last[cellid] = i + 1;
				}
			}


			// reset initial contact network
			for (i=0; i<NPW; i++)
				cij[i] = 0;

			// check new contacts
			for (ci=0; ci<NCELLS; ci++){

				// get start of list of particles
				pi = head[ci];

				// loop over linked list
				while(pi > 0){
					// real particle index
					i = pi - 1;

					// next in list
					pj = list[pi];

					// loop down neighbors of pi - 1 in same cell
					while(pj > 0){
						// real index of pj
						j = pj - 1;

						// check if dimer pair
						if ((i % 2 == 0 && j == i+1) || (j % 2 == 0 && i == j+1) ){
							// update pj and continue
							pj = list[pj];
							continue;
						}

						// contact distance
						sij = radii[i] + radii[j];
						sij *= 1 + l2;

						// particle distance
						rij = 0.0;
						for (d=0; d<NDIM; d++){
							// distance across box
							dp = pos[NDIM*j + d] - pos[NDIM*i + d];

							// image distance
							dp -= L[d]*round(dp/L[d]);

							// add to total distance
							rij += dp*dp;

							// save vector component
							vij[d] = dp;
						}

						// check for overlap, add to contacts
						if (rij < sij*sij)
							addContacts(i,j,cij,z,NMTOT);

						// update pj
						pj = list[pj];
					}

					// test overlaps with forward neighboring cells
					for (cj=0; cj<NNN; cj++){
						// get first particle in neighboring cell
						pj = head[nn[ci][cj]];

						// loop over cells in neighboring cells
						while(pj > 0){
							// real index of pj
							j = pj - 1;	

							// check if dimer pair
							if ((i % 2 == 0 && j == i+1) || (j % 2 == 0 && i == j+1) ){
								// update pj and continue
								pj = list[pj];
								continue;
							}

							// contact distance
							sij = radii[i] + radii[j];
							sij *= 1 + l2;

							// particle distance
							rij = 0.0;
							for (d=0; d<NDIM; d++){
								// distance across box
								dp = pos[NDIM*j + d] - pos[NDIM*i + d];

								// image distance
								dp -= L[d]*round(dp/L[d]);

								// add to total distance
								rij += dp*dp;

								// save vector component
								vij[d] = dp;
							}

							// check for overlap, add to contacts
							if (rij < sij*sij)
								addContacts(i,j,cij,z,NMTOT);

							// update pj
							pj = list[pj];
						}
					}

					// update pi
					pi = list[pi];
				}
			} 


			// force computation using cell linked-lists
			U = 0;
			fill(S.begin(),S.end(),0.0);
			for (ci=0; ci<NCELLS; ci++){

				// get start of list of particles
				pi = head[ci];

				// loop over linked list
				while(pi > 0){
					// real particle index
					i = pi - 1;

					// next in list
					pj = list[pi];

					// loop down neighbors of pi - 1 in same cell
					while(pj > 0){
						// real index of pj
						j = pj - 1;	

						// check if dimer pair
						if ((i % 2 == 0 && j == i+1) || (j % 2 == 0 && i == j+1) ){
							// update pj and continue
							pj = list[pj];
							continue;
						}

						// compute forces if previous overlapping
						if (cij[cmindex(i,j,NMTOT)]){
							// contact distance
							sij = radii[i] + radii[j];

							// particle distance
							rij = 0.0;
							for (d=0; d<NDIM; d++){
								// distance across box
								dp = pos[NDIM*j + d] - pos[NDIM*i + d];

								// image distance
								dp -= L[d]*round(dp/L[d]);

								// add to total distance
								rij += dp*dp;

								// save vector component
								vij[d] = dp;
							}
							rij = sqrt(rij);

							// dimensionless distance
		                    h = rij/sij;
		                    
		                    // lij: effective adhesion strength based on contacts
		                    // lij = 0.25*l2*((1/z[i]) + (1/z[j]));
		                    lij = 0.75*l2;
		                    p1 = 1.0 + lij;
		                    u1 = lij/(l2 - lij);
		                    u2 = l2*lij;
		                    
		                    // attractive force
		                    if (h < p1){
		                        // scalar part of force
		                        ftmp = (1 - h)/sij;

		                        // potential
		                        U = U + (0.5*pow(1 - h,2.0) - u2);
		                    }
		                    else{
		                       	// scalar part of force
		                        ftmp = u1*(h - 1 - l2)/sij;

		                        // potential
		                        U = U - 0.5*u1*pow(h - 1 - l2,2.0);
		                    }

							// add to forces & normal stresses
							for (d=0; d<NDIM; d++){
								// force on i due to j
								fdirtmp = -ftmp*(vij[d]/rij);

								// forces (equal and opposit)
								F[NDIM*i + d] += fdirtmp;
								F[NDIM*j + d] -= fdirtmp;

								// normal stresses
								S[d] -= fdirtmp*vij[d];
							}

							// off diagonal stresses (sorted XY, XZ, YZ)
							S[3] += vij[0]*((ftmp*vij[1])/rij);
							S[4] += vij[0]*((ftmp*vij[2])/rij);
							S[5] += vij[1]*((ftmp*vij[2])/rij);
						}

						// update pj
						pj = list[pj];
					}

					// test overlaps with forward neighboring cells
					for (cj=0; cj<NNN; cj++){
						// get first particle in neighboring cell
						pj = head[nn[ci][cj]];

						// loop over cells in neighboring cells
						while(pj > 0){
							// real index of pj
							j = pj - 1;	

							// check if dimer pair
							if ((i % 2 == 0 && j == i+1) || (j % 2 == 0 && i == j+1) ){
								// update pj and continue
								pj = list[pj];
								continue;
							}

							// compute forces if previous overlapping
							if (cij[cmindex(i,j,NMTOT)]){
								// contact distance
								sij = radii[i] + radii[j];

								// particle distance
								rij = 0.0;
								for (d=0; d<NDIM; d++){
									// distance across box
									dp = pos[NDIM*j + d] - pos[NDIM*i + d];

									// image distance
									dp -= L[d]*round(dp/L[d]);

									// add to total distance
									rij += dp*dp;

									// save vector component
									vij[d] = dp;
								}
								rij = sqrt(rij);

								// dimensionless distance
			                    h = rij/sij;
			                    
			                    // lij: effective adhesion strength based on contacts
			                    // lij = 0.25*l2*((1/z[i]) + (1/z[j]));
			                    lij = 0.75*l2;
			                    p1 = 1.0 + lij;
			                    u1 = lij/(l2 - lij);
			                    u2 = l2*lij;
			                    
			                    // attractive force
			                    if (h < p1){
			                        // scalar part of force
			                        ftmp = (1 - h)/sij;

			                        // potential
			                        U = U + (0.5*pow(1 - h,2.0) - u2);
			                    }
			                    else{
			                       	// scalar part of force
			                        ftmp = u1*(h - 1 - l2)/sij;

			                        // potential
			                        U = U - 0.5*u1*pow(h - 1 - l2,2.0);
			                    }

								// add to forces & normal stresses
								for (d=0; d<NDIM; d++){
									// force on i due to j
									fdirtmp = -ftmp*(vij[d]/rij);

									// forces (equal and opposit)
									F[NDIM*i + d] += fdirtmp;
									F[NDIM*j + d] -= fdirtmp;

									// normal stresses
									S[d] -= fdirtmp*vij[d];
								}

								// off diagonal stresses (sorted XY, XZ, YZ)
								S[3] += vij[0]*((ftmp*vij[1])/rij);
								S[4] += vij[0]*((ftmp*vij[2])/rij);
								S[5] += vij[1]*((ftmp*vij[2])/rij);
							}

							// update pj
							pj = list[pj];
						}
					}

					// update pi
					pi = list[pi];
				}
			}

			// Add dimer bond force
			for (i=0; i<N; i++){
				// distance between pair images
				lx = pos[2*NDIM*i + 3] - pos[2*NDIM*i];
				lx -= L[0]*round(lx/L[0]);

				ly = pos[2*NDIM*i + 4] - pos[2*NDIM*i + 1];
				ly -= L[1]*round(ly/L[1]);

				lz = pos[2*NDIM*i + 5] - pos[2*NDIM*i + 2];
				lz -= L[2]*round(lz/L[2]);

				// bond length
				l = lx*lx + ly*ly + lz*lz;

				// unit vector
				ux = lx/l;
				uy = ly/l;
				uz = lz/l;

				// force on 1 due to 2
				ftmp = kl*(l - l0[i]);
				flx = ftmp*ux;
				fly = ftmp*uy;
				flz = ftmp*uz;

				// bonded force on monomer 1
				F[2*NDIM*i] += flx;
				F[2*NDIM*i + 1] += fly;
				F[2*NDIM*i + 2] += flz;

				// bonded force on monomer 2
				F[2*NDIM*i + 3] -= flx;
				F[2*NDIM*i + 4] -= fly;
				F[2*NDIM*i + 5] -= flz;

				// add to potential energy
				U += 0.5*kl*pow(l - l0[i],2.0);

				// // add to stress tensor (note minus sign, keeps things positive)
				// S[0] -= lx*flx;
				// S[1] -= ly*fly;
				// S[2] -= lz*flz;

				// // sorted XY, XZ, YZ
				// S[3] -= lx*fly;
				// S[4] -= ly*flz;
				// S[5] -= ly*flz;
			}	

			// VELOCITY-VERLET UPDATE 2: VELOCITIES AND ACCELERATION
			for (i=0; i<NMTOT; i++){
				for(d=0; d<NDIM; d++){
					// dof index
					ind = NDIM*i + d;

					// new acceleration
					atmp = F[ind]/mass[i];

					// update velocities
					v[ind] += 0.5*(atmp + a[ind])*dt;

					// update acceleration
					a[ind] = atmp;
				}
			}

			// compute fnorm, vnorm and P
			fnorm = 0.0;
			vnorm = 0.0;
			P = 0.0;
			for (i=0; i<NMTOT; i++){
				for (d=0; d<NDIM; d++){
					ind = NDIM*i + d;
					fnorm 	+= F[ind]*F[ind];
					vnorm 	+= v[ind]*v[ind];
					P 		+= v[ind]*F[ind];
				}
			}
			fnorm = sqrt(fnorm);
			vnorm = sqrt(vnorm);

			// update fcheck based on fnorm (= force per degree of freedom)
			fcheck = fnorm/(NDIM*NMTOT);

			// update npPMin
			if (fcheck < Ftol)
				npPMin++;
			else
				npPMin = 0;

			// print to console
			if (fireit % NSKIP == 0){
				cout << endl << endl;
				cout << "===========================================" << endl;
				cout << "		G E L A T I O N 			" << endl;
				cout << " 	F I R E 						" << endl;
				cout << "		M I N I M I Z A T I O N 	" << endl;
				cout << "===========================================" << endl;
				cout << endl;
				cout << "	** fireit = " << fireit << endl;
				cout << "	** it = " << it << endl;
				cout << "	** fcheck = " << fcheck << endl;
				cout << "	** vnorm = " << vnorm << endl;
				cout << "	** dt = " << dt << endl;
				cout << "	** P = " << P << endl;
				cout << "	** Pdir = " << P/(fnorm*vnorm) << endl;
				cout << "	** alpha = " << alpha << endl;
				cout << "	** U = " << U << endl;
				cout << "	** phi = " << phi << endl << endl;
			}

			// Step 1. adjust simulation based on net motion of degrees of freedom
			if (P > 0){
				// increase positive counter
				npPos++;

				// reset negative counter
				npNeg = 0;

				// alter simulation if enough positive steps have been taken
				if (npPos > NMIN){
					// change time step
					if (dt*finc < dtmax)
						dt *= finc;

					// decrease alpha
					alpha *= falpha;
				}
			}
			else{
				// reset positive counter
				npPos = 0;

				// increase negative counter
				npNeg++;

				// check if simulation is stuck
				if (npNeg > NNEGMAX){
					cout << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
					return 1;
				}

				// take half step backwards, reset velocities
				for (i=0; i<NMTOT; i++){
					for (d=0; d<NDIM; d++){
						// dof index
						ind = NDIM*i + d;

						// take half step backwards
						pos[ind] -= 0.5*dt*v[ind];

						// reset velocities
						v[ind] = 0;
					}
				}

				// decrease time step if past initial delay
				if (fireit > NMIN){
					// decrease time step 
					if (dt*fdec > dtmin)
						dt *= fdec;

					// reset alpha
					alpha = alpha0;
				}
			}


			// update velocities (s.d. vs inertial dynamics) only if forces are acting
			if (fnorm > 0){
				for (i=0; i<NMTOT; i++){
					for (d=0; d<NDIM; d++){
						// dof index
						ind = NDIM*i + d;

						// update velocity
						v[ind] = (1 - alpha)*v[ind] + alpha*(vnorm/fnorm)*F[ind];
					}
				}
			}

			// update iterator
			fireit++;
		}

		// print + continue regardless of convergence
		cout << endl << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===========================================" << endl;
		cout << " 	F I R E 						" << endl;
		cout << "		M I N I M I Z A T I O N 	" << endl;
		cout << "	 C O N V E R G E D! 			" << endl;
		cout << "===========================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	** fireit = " << fireit << endl;
		cout << "	** it = " << it << endl;
		cout << "	** fcheck = " << fcheck << endl;
		cout << "	** vnorm = " << vnorm << endl;
		cout << "	** dt = " << dt << endl;
		cout << "	** P = " << P << endl;
		cout << "	** Pdir = " << P/(fnorm*vnorm) << endl;
		cout << "	** alpha = " << alpha << endl;
		cout << "	** U = " << U << endl;
		cout << "	** phi = " << phi << endl << endl;

		// print positions to xyz file during initial minimization
		if (it % NPHISKIP == 0)
			printPos(posout,pos,radii,L,S,z,NMTOT);

		// incremement iterator
		it++;

		// alternate volumetric strain and decompression
		if (it % 2 == 0){
			// even: deform boundary
			L[0] *= dgx;
			L[1] *= dgy;
			L[2] *= dgz;

			// affine displacements
			for (i=0; i<NMTOT; i++){
				pos[i*NDIM] *= dgx;
				pos[i*NDIM + 1] *= dgy;
				pos[i*NDIM + 2] *= dgz;
			}

			// update linked list cell geometry
			for (d=0; d<NDIM; d++)
				lc[d] = L[d]/sc[d];

			// if any lc is smaller than what box size should be, reset 
			rmax = 0.0;
			for (i=0; i<N; i++){
				if (radii[2*i] > rmax)
					rmax = radii[2*i];
			}

			resetCLL = 0;
			boxsize = 2.1*rmax;
			for (d=0; d<NDIM; d++){
				if (lc[d] < boxsize){
					resetCLL = 1;
					break;
				}
			}

			// check if L is too small
			Lmin = 1e6;
			for (d=0; d<NDIM; d++){
				if (L[d] < Lmin)
					Lmin = L[d];
			}

			// if too small, stop box deformation
			if (Lmin < 0.5*L0){
				dgx = 1.0;
				dgy = 1.0;
				dgz = 1.0;
			}

			if (resetCLL){
				cout << "** Need to reset box size, lc < boxsize" << endl;
				cellLinkedList(sc,lc,nn,NCELLS,L,boxsize);
				cout << "sc[0]=" << sc[0] << ", ";
				cout << "sc[1]=" << sc[1] << ", ";
				cout << "sc[2]=" << sc[2] << ", ";
				cout << endl;
				cout << "lc[0]=" << lc[0] << ", ";
				cout << "lc[1]=" << lc[1] << ", ";
				cout << "lc[2]=" << lc[2] << ", ";
				cout << endl;
				head.resize(NCELLS);
				last.resize(NCELLS);
			}
		}
		else{
			// odd: shrink particles
			rscale = pow((phi - dphi)/phi,1.0/NDIM);
			for (i=0; i<NMTOT; i++){
				mass[i] *= pow(rscale,NDIM);
				radii[i] *= rscale;
			}
			for (i=0; i<N; i++)
				l0[i] *= rscale;
		}

		// recompute phi
		vp = 0.0;
		for (i=0; i<N; i++){
			r1 = radii[2*i];
			vp += (8.0/3.0)*PI*pow(r1,3.0);
			if (l0[i] < 2.0*r1)
				vp -= (PI/12.0)*(4.0*r1 + l0[i])*pow(2.0*r1 - l0[i],2.0);
		}
		phi = vp/(L[0]*L[1]*L[2]);
	}


	// close objects
	posout.close();

	// print to console, return
	cout << "	** FINISHED MAIN FOR sphereGel.cpp, ENDING." << endl << endl << endl;
	return 0;
}









/* 

	&&&&&&&&&&&&&&&&&&&&&&& FUNCTION DEFINITIONS &&&&&&&&&&&&&&&&&&&&&

	FUNCTIONS DEFINED

	printPos 		: output particle positions to .pos file for processing and visualization
	cmindex 		: map (i,j) pairs to contact vector
	addContacts 	: assign contacts to contact matrix and z vector

	&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 

*/



void printPos(ofstream &posout, vector<double> &pos, vector<double> &radii, vector<double> &L, vector<double> &S, vector<int> &z, int N){
	// local variables
	int i, d;

	// check to see if file is open
	if (!posout.is_open()){
		cout << "	** ERROR: IN outputxyz(), input xyzout has no open file, ending." << endl;
		exit(1);
	}
	else
		cout << "	** Printing xyz information to file." << endl;

	// output xyz info
	posout << N << endl;
	posout << L[0] << "  " << L[1] << "  " << L[2] << endl;
	posout << S[0] << "  " << S[1] << "  " << S[2] << endl;
	posout << S[3] << "  " << S[4] << "  " << S[5] << endl;

	// loop over particles, print positions to file
	for (i=0; i<N; i++){
		// output particle positions
		for (d=0; d<NDIM; d++)
			posout << setw(wnum) << setprecision(pnum) << pos[NDIM*i + d];

		// output particle radius
		posout << setw(wnum) << setprecision(pnum) << radii[i];

		// output number of contacts
		posout << setw(wnum) << z[i] << endl;
	}
}

int cmindex(int i, int j, int N){
	if (i > j)
		return N*j + i - (j+1)*(j+2)/2;
	else
		return N*i + j - (i+1)*(i+2)/2;
}


void addContacts(int i, int j, vector<bool> &cij, vector<int> &z, int N){
	// local variables
	int ci;

	// determine contact index
	if (i > j)
		ci = N*j + i - (j+1)*(j+2)/2;
	else
		ci = N*i + j - (i+1)*(i+2)/2;

	// add to matrix
	cij[ci] = 1;

	// add to particle contact count
	z[i]++;
	z[j]++;
}




// function to set / reset cell linked list based on particle size
void cellLinkedList(vector<int> &sc, vector<double> &lc, vector< vector<int> > &nn, int &NCELLS, vector<double> &L, double boxsize){
	// local variables
	int d, i;

	// cell box lengths in each direction
	NCELLS = 1;
	for (d=0; d<NDIM; d++){
		// determine number of cells along given dimension by rmax
		sc[d] = round(L[d]/boxsize);

		// just in case, if < 3, change to 3 so cell neighbor checking will work
		if (sc[d] < 3)
			sc[d] = 3;

		// determine cell box length by number of cells
		lc[d] = L[d]/sc[d];

		// count total number of cells
		NCELLS *= sc[d];
	}

	// neighbor cells for each cell (26 neighbors / cell)
	int zi, scx, scy, scxy;
	int nntmp0, nntmp1, nntmp2, nntmp3, nntmp4;
	nn.clear();
	nn.resize(NCELLS);

	scx = sc[0];
	scy = sc[1];
	scxy = scx*scy;

	// loop over cells, save neighbor for each cell index
	for (i=0; i<NCELLS; i++){
		// reshape entry
		nn[i].resize(NNN);

		// z sclice
		zi = i/scxy;		
		
		// faces
		nntmp0 		= (i+NCELLS-1) % NCELLS; 			// left neighbor (i-1)
		nntmp1 		= i-scx;							// bottom neighbor (j-1)
		nntmp2 		= (i+NCELLS-scxy) % NCELLS;			// backward neighbor (k-1)
		nn[i][0] 	= (i+1) % NCELLS; 					// right neighbor (i+1)
		nn[i][1] 	= i+scx;							// top neighbor (j+1)
		nn[i][2] 	= (i+scxy) % NCELLS;				// forward neighbor (k+1)		

		// y-direction bc
		if (i-zi*scxy < scx)							// if on bottom row, look up (y direction)
			nntmp1 = i-scx+scxy;	
		if ((i+scx)/(scxy) > zi)						// if on top row, look down (y direction)
			nn[i][1] = i-scxy+scx;	

		// edges

		// * xy plane
		nn[i][3] = nntmp1 + 1;							// i+1, j-1
		nn[i][4] = nn[i][1] + 1;						// i+1, j+1

		// * xz plane
		nn[i][5] = (nntmp2 + 1) % NCELLS;				// i+1, k-1
		nn[i][6] = (nn[i][2] + 1) % NCELLS;				// i+1, k+1

		// * yz plane
		nntmp3 		= nntmp2 - scx;						// j-1, k-1
		nn[i][7] 	= nntmp2 + scx;						// j+1, k-1
		nntmp4  	= nn[i][2] - scx;					// j-1, k+1
		nn[i][8] 	= nn[i][2] + scx; 					// j+1, k+1

		// y-direction bc
		if (i-zi*scxy < scx){							// if on bottom row, look up (y direction)
			nntmp3 = nntmp2-scx+scxy;
			nntmp4 = nn[i][2]-scx+scxy;	
		}
		if ((i+scx)/scxy > zi){							// if on top row, look down (y direction)
			nn[i][7] = nntmp2-scxy+scx;
			nn[i][8] = nn[i][2]-scxy+scx;			
		}
		


		// cubic vertices
		nn[i][9] = (nntmp3 + 1) % NCELLS; 				// i+1, j-1, k-1
		nn[i][10] = (nn[i][7] + 1) % NCELLS; 			// i+1, j+1, k-1
		nn[i][11] = (nntmp4 + 1) % NCELLS; 				// i+1, j-1, k+1
		nn[i][12] = (nn[i][8] + 1) % NCELLS; 			// i+1, j+1, k+1



		// right BC
		if ((i+1) % scx == 0){
			nn[i][0] = i-scx+1;
			nn[i][3] = nntmp1-scx+1;
			nn[i][4] = nn[i][1]-scx+1;
			nn[i][5] = nntmp2-scx+1;
			nn[i][6] = nn[i][2]-scx+1;
			nn[i][9] = nntmp3-scx+1;
			nn[i][10] = nn[i][7]-scx+1;
			nn[i][11] = nntmp4-scx+1;
			nn[i][12] = nn[i][8]-scx+1;
		}
	} 
}


