/*

	Jack Treado
	06/2020
	sphereGel.cpp

	Generate gel of spheres via athermal extension

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
const int pnum 				= 14;

const double phi0 			= 0.9;
const double phimin 		= 0.1;
const double timeStepMag 	= 0.01;

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
void printXYZ(ofstream& xyzout, vector<double>& pos, vector<double>& radii, vector<double>& L, int N);
void printCM(ofstream& cmout, vector<bool>& cij, int NPW);
int cmindex(int i, int j, int N);
void addContacts(int i, int j, vector<bool>& cij, vector<int>& z, int N);

// MAIN
int main(int argc, char const *argv[])
{
	// variables used throughout main
	int i, j, d, ind;
	double sij, rij, dp, ftmp, atmp;
	vector<double> vij(NDIM,0.0);

	// parameters to be read in 
	int N, seed;
	double l2, dr, dphi, dlz, Ftol;

	// read in parameters from command line input
	string N_str 		= argv[1];
	string dr_str 		= argv[2];
	string dphi_str 	= argv[3];
	string dlz_str 		= argv[4];
	string l2_str 		= argv[5];
	string Ftol_str 	= argv[6];
	string seed_str 	= argv[7];
	string xyzFile 		= argv[8];
	string cmFile 		= argv[9];

	stringstream Nss(N_str);
	stringstream drss(dr_str);
	stringstream dphiss(dphi_str);
	stringstream dlzss(dlz_str);
	stringstream l2ss(l2_str);
	stringstream Ftolss(Ftol_str);
	stringstream seedss(seed_str);

	Nss >> N;
	drss >> dr;
	dphiss >> dphi;
	dlzss >> dlz;
	l2ss >> l2;
	Ftolss >> Ftol;
	seedss >> seed;

	// open xyz file
	ofstream xyzout;
	xyzout.open(xyzFile.c_str());
	if (!xyzout.is_open()){
		cout << "	** ERROR: xyz file " << xyzFile << " could not be opened, ending." << endl;
		return 1;
	}

	// open contact matrix file
	ofstream cmout;
	cmout.open(cmFile.c_str());
	if (!cmout.is_open()){
		cout << "	** ERROR: cm file " << cmFile << " could not be opened, ending." << endl;
		return 1;
	}


	// output opening statement to console
	cout << "=======================================================" << endl << endl;
	cout << "		sphereGel.cpp 									" << endl;
	cout << "		Jack Treado, 2020   							" << endl << endl;
	cout << "		Athermal gelation of densely-packed spheres 	" << endl;
	cout << "		N 			= " << N << "						" << endl;
	cout << "		dr 			= " << dr << " 						" << endl;
	cout << "		dphi 		= " << dphi << " 					" << endl;
	cout << "		dlz 		= " << dlz << "						" << endl;
	cout << "		l2 			= " << l2 << " 						" << endl;
	cout << "		Ftol 		= " << Ftol << " 					" << endl;
	cout << "		seed 		= " << seed << "					" << endl;
	cout << "		xyz file 	= " << xyzFile << "					" << endl;
	cout << "		cm file 	= " << cmFile << " 					" << endl;
	cout << endl;
	cout << "=======================================================" << endl << endl;


	// seed random number generator
	srand48(seed);


	/* * * * * * * * * * * * * * * * * * 

				INITIALIZE

					PARTICLES

	 * * * * * * * * * * * * * * * * * */

	// initialization variables
	double r1, r2, grv, rsum3, rmax, L0;

	// particle radii and positions
	vector<double> radii(N,0.0);
	vector<double> mass(N,0.0);
	vector<double> pos(NDIM*N,0.0);

	// gaussian-distributed polydispersity
	rsum3 = 0.0;
	rmax = 0.0;
	for (i=0; i<N; i++){
		// generate random numbers
		r1 = drand48();
		r2 = drand48();

		// calculate gaussian random variable using Box-Muller transform
		grv = sqrt(-2.0*log(r1))*cos(2*PI*r2);

		// get radius
		radii[i] = grv*dr + 1.0;

		// compute mass
		mass[i] = (4.0/3.0)*PI*pow(radii[i],3.0);

		// add to cubed rsum
		rsum3 += pow(radii.at(i),3.0);

		// determine max radius
		if (radii[i] > rmax)
			rmax = radii[i];
	}


	// box lengths (initially a cube)
	L0 = pow(4.0*PI*rsum3/(3.0*phi0),1.0/3.0);
	vector<double> L(NDIM,L0);

	// initialize particle positions randomly throughout box
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++)
			pos[NDIM*i + d] = L[d]*drand48();
	}


	// print positions to xyz file
	printXYZ(xyzout,pos,radii,L,N);


	// determine fundamental MD time unit
	double dtMD, dt0, dt;

	dtMD 	= 1.0;
	dt0 	= timeStepMag*dtMD;
	dt 		= dt0;







	/* * * * * * * * * * * * * * * * * * 

			INITIALIZE

				CELL LINKED LIST

	 * * * * * * * * * * * * * * * * * */

	// Cell linked-list variables

	// cell box lengths in each direction
	vector<int> sc(NDIM,0);
	vector<double> lc(NDIM,0.0);
	int NCELLS = 1;
	for (d=0; d<NDIM; d++){
		// determine number of cells along given dimension by rmax
		sc[d] = round(L[d]/(2.1*rmax));

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
	vector< vector<int> > nn;
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

	// linked-list variables
	vector<int> head(NCELLS,0);
	vector<int> last(NCELLS,0);
	vector<int> list(N+1,0);




	/* * * * * * * * * * * * * * * * * * 

			INITIAL FIRE

					MINIMZATION

	 * * * * * * * * * * * * * * * * * */


	// initialize velocity and force information
	vector<double> v(NDIM*N,0.0);
	vector<double> a(NDIM*N,0.0);
	vector<double> F(NDIM*N,0.0);

	// FIRE VARIABLES
	double P  		= 0;	
	double fnorm 	= 0;
	double vnorm 	= 0;
	double alpha   	= alpha0;

	double dtmax   	= 10*dt0;
	double dtmin   	= 1e-8*dt0;
	dt 				= dt0;

	int npPos      	= 0;
	int npNeg      	= 0;
	int npPMin      = 0;

	int fireit    	= 0;
	double fcheck  	= 10*Ftol;
	double U 		= 0;

	// cell linked-list variables
	int cellid, ci, cj, pi, pj, sctmp;

	// loop until force relaxes
	while ((fcheck > Ftol || npPMin < NMIN) && fireit < itmax){

		// VELOCITY-VERLET UPDATE 1: POSITIONS
		for (i=0; i<N; i++){
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
		list[N] = 0;


		// reset linked list head
		for (i=0; i<NCELLS; i++){
			head[i] = 0;
			last[i] = 0;
		}

		// sort particles into linked list
		for (i=0; i<N; i++){
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


		// VELOCITY-VERLET UPDATE 2: VELOCITIES AND ACCELERATION
		for (i=0; i<N; i++){
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
		for (i=0; i<N; i++){
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
		fcheck = fnorm/(NDIM*N);

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
			for (i=0; i<N; i++){
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
			for (i=0; i<N; i++){
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
	int NPW = N*(N-1)/2;
	vector<bool> cij(NPW,0);
	vector<int> z(N,0);

	// instantaneous phi
	double phi = 0.0;
	for (i=0; i<N; i++)
		phi += mass[i];
	phi /= L[0]*L[1]*L[2];

	// decompression iterator
	int it = 0;

	// attraction parameters
	double p1, u1, u2, h, lij;

	// determine change in Lz
	int NSTEPS = round((phi - phimin)/dphi);
	double lzscale, rscale;

	// determine frames to skip between xyz output
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
			for (i=0; i<N; i++){
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
			list[N] = 0;


			// reset linked list head
			for (i=0; i<NCELLS; i++){
				head[i] = 0;
				last[i] = 0;
			}

			// sort particles into linked list
			for (i=0; i<N; i++){
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
							addContacts(i,j,cij,z,N);

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
								addContacts(i,j,cij,z,N);

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

						// compute forces if previous overlapping
						if (cij[cmindex(i,j,N)]){
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
		                    lij = 0.25*l2*((1/z[i]) + (1/z[j]));
		                    p1 = 1.0 + lij;
		                    u1 = lij/(l2 - lij);
		                    u2 = l2*lij;
		                    
		                    // attractive force
		                    if (h < p1){
		                        // scalar part of force
		                        ftmp = 1 - h;

		                        // potential
		                        U = U + (0.5*sij*pow(1 - h,2.0) - u2);
		                    }
		                    else{
		                       	// scalar part of force
		                        ftmp = u1*(h - 1 - l2);

		                        // potential
		                        U = U - 0.5*sij*u1*pow(h - 1 - l2,2.0);
		                    }

							// add to forces
							for (d=0; d<NDIM; d++){
								F[NDIM*i + d] -= ftmp*(vij[d]/rij);
								F[NDIM*j + d] += ftmp*(vij[d]/rij);
							}
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

							// compute forces if previous overlapping
							if (cij[cmindex(i,j,N)]){
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
			                    lij = 0.25*l2*((1/z[i]) + (1/z[j]));
			                    p1 = 1.0 + lij;
			                    u1 = lij/(l2 - lij);
			                    u2 = l2*lij;
			                    
			                    // attractive force
			                    if (h < p1){
			                        // scalar part of force
			                        ftmp = 1 - h;

			                        // potential
			                        U = U + 0.5*sij*pow(1 - h,2.0) - u2;
			                    }
			                    else{
			                       	// scalar part of force
			                        ftmp = u1*(h - 1 - l2);

			                        // potential
			                        U = U - 0.5*sij*u1*pow(h - 1 - l2,2.0);
			                    }

								// add to forces
								for (d=0; d<NDIM; d++){
									F[NDIM*i + d] -= ftmp*(vij[d]/rij);
									F[NDIM*j + d] += ftmp*(vij[d]/rij);
								}
							}

							// update pj
							pj = list[pj];
						}
					}

					// update pi
					pi = list[pi];
				}
			}


			// VELOCITY-VERLET UPDATE 2: VELOCITIES AND ACCELERATION
			for (i=0; i<N; i++){
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
			for (i=0; i<N; i++){
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
			fcheck = fnorm/(NDIM*N);

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
				for (i=0; i<N; i++){
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
				for (i=0; i<N; i++){
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
			if (it % NPHISKIP == 0){
				printXYZ(xyzout,pos,radii,L,N);
				printCM(cmout,cij,NPW);
			}
		}

		// incremement iterator
		it++;

		// scale particles and box lengths based on delta phi
		lzscale = 1.0 + (dlz - 1.0)*(L0/L[2])/NSTEPS;
		L[2] *= lzscale;
		lc[2] *= lzscale;

		rscale = pow(lzscale*(1.0 - (dphi/phi)),1.0/NDIM);
		for (i=0; i<N; i++){
			mass[i] *= pow(rscale,NDIM);
			radii[i] *= rscale;
		}

		// recompute phi
		phi = 0.0;
		for (i=0; i<N; i++)
			phi += mass[i];
		phi /= L[0]*L[1]*L[2]; 


		// apply affine strain
		for (i=0; i<N; i++)
			pos[i*NDIM + 2] *= lzscale;
	}

	


	// close objects
	xyzout.close();
	cmout.close();

	// print to console, return
	cout << "	** FINISHED MAIN FOR sphereGel.cpp, ENDING." << endl << endl << endl;
	return 0;
}









/* 

	&&&&&&&&&&&&&&&&&&&&&&& FUNCTION DEFINITIONS &&&&&&&&&&&&&&&&&&&&&

	FUNCTIONS DEFINED

	printXYZ 		: output particle positions to .xyz file for processing and visualization
	printCM			: output contact matrix in a row vector / frame
	cmindex 		: map (i,j) pairs to contact vector
	addContacts 	: assign contacts to contact matrix and z vector

	&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 

*/



void printXYZ(ofstream& xyzout, vector<double>& pos, vector<double>& radii, vector<double>& L, int N){
	// local variables
	int i, d;
	char atom = 'C';

	// check to see if file is open
	if (!xyzout.is_open()){
		cout << "	** ERROR: IN outputxyz(), input xyzout has no open file, ending." << endl;
		exit(1);
	}
	else
		cout << "	** Printing xyz information to file." << endl;

	// output xyz info
	xyzout << N << endl;
	xyzout << "Lattice=\"" << L[0] << " 0.0 0.0 0.0 " << L[1] << " 0.0 0.0 0.0 " << L[2] << "\" " << '\t';
	xyzout << "Properties=species:S:1:pos:R:3:radius:R:1" << endl;

	// loop over particles, print positions to file
	for (i=0; i<N; i++){
		// output particle type
		xyzout << setw(w) << atom;

		// output particle positions
		for (d=0; d<NDIM; d++)
			xyzout << setw(wnum) << setprecision(pnum) << pos[NDIM*i + d];

		// output particle radius
		xyzout << setw(wnum) << setprecision(pnum) << radii[i] << endl;
	}
}


void printCM(ofstream& cmout, vector<bool>& cij, int NPW){
	// local variables
	int ci;

	for (ci=0; ci<NPW; ci++)
		cmout << setw(w) << cij[ci];
	cmout << endl;
}


int cmindex(int i, int j, int N){
	if (i > j)
		return N*j + i - (j+1)*(j+2)/2;
	else
		return N*i + j - (i+1)*(i+2)/2;
}


void addContacts(int i, int j, vector<bool>& cij, vector<int>& z, int N){
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



