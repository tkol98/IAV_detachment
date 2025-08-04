#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdbool.h>
#include <chrono>
#include <random>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <iterator>

#define MATRIX_INT std::vector<std::vector<int>>
#define MATRIX_DOUBLE std::vector<std::vector<double>>

//VALEURS POUR RAN2
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

// ------------------------------------------------------------------------------------------------------------------
// -- SYSTEM PARAMETERS
const double Taille_Systeme= 1000.0;	// nm, system size, should be increased soon
const double roh_SA= 0.03; 				// /(nm^2)
//int N_R= 200000;						// number of receptors, => determined by density roh_SA
const int N_R= roh_SA * Taille_Systeme * Taille_Systeme;		// number of receptors

// -- PARTICLE PARAMETERS
const double L= 300.0;					// nm, length of virus (for non-spherical)
//const double diameter= 120.0; 			// nm, diameter of virus
constexpr double r_sph= 60.0;	 				// nm, radius of spherical virus
// INTERACTIONS:
constexpr double k_on= 0.56 * 1e3;				// 1/s, binding rate HA
constexpr double k_off= 9. * 1e3;				// 1/s, unbinding rate HA
//constexpr double k_d= 2.0 * 1e3;				// 1/s, NA cleaving rate
double k_d;	//1/s, NA cleaving rate,  for external input

// -- SIMULATION PARAMETERS
constexpr double D= 7.5; 						// nm, maximal interaction distance (λ)
constexpr double R_Verlet= 10.0; 				// nm, verlet radius
constexpr double R_close= 1.5* r_sph; 			// nm, Radius of the circle containing all the "close" receptors

constexpr double DT= 1.0 * 1e-6;				// seconds -> 1 µs, time step size
constexpr int N_dT= 25000;						// number of time steps

// -- DIFFUSION PARAMETERS
constexpr double D_receptor= 373333.;						// (nm^2)/s,
constexpr double D_sphere= 373333.;							// (nm^2)/s, diffusion constant D* of spherical particle	
constexpr double D_r= (3./4.) * D_sphere/(r_sph*r_sph);		// (rad^2)/s, rotational diffusion constant
constexpr double DT_Brown= 1.0 * 1e-8;					// s, step size of Brownian diffusion
const int N_Brown= int(DT/DT_Brown);					// number of Brownian steps between reactions
// ---------------
int rejections=0;
int total_tests=0;
// ------------------------------------------------------------------------------------------------------------------

// ---------------- HELPER FUNCTIONS -------------------------------
bool isodd(int num) 
{// Verifies whether a number is odd
	if(num % 2 == 0)
		return false;
	else
		return true;
}

double distance2(std::vector<double> r1,std::vector<double> r2)
{// squared distance between two vectors 
	double dist2;

	dist2=pow((r1[0]-r2[0]),2)+pow((r1[1]-r2[1]),2);

	return dist2;
}

double ran2(long *idum)
{// a random number between 0 and 1 generator
	// idum should be negative
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    
    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

double RGV(long *g)
{/*	Returns a normally distributed random number between -1 and 1 using the 
	* Box-Muller Method 
	*/
	
	std::vector<double> random_numbers;
	double p;
	
	// get two (not norm. distr.) random numbers between 0 and 1
    for(int j=0; j<2; j++)
    {
        p = ran2(g);
        random_numbers.push_back(p);
    }
 	double r1 = random_numbers[0];
	double r2 = random_numbers[1];

  	double r3 = sqrt(-2*(log(1-r1)));	// Box-Muller "radius"
  	double r = cos(r2*2*M_PI)*r3;		// real part of the complex Box-muller number 

  	return r;
}


double Get_time_reaction(double a_tot, long *idum)
{// sample a reaction time (Gillespie algorithm)
    double t;
    t =  (-1)*(log(1 - ran2(idum)))/a_tot;
    
    return t;
}

int Cell_to_Index(int N_x, int N_y, int ix, int iy) 
{//gets cell coordinates as input and returns cell index
	int ixp;
	int iyp;
	int nc;

	ixp= ix;
	iyp = iy;

	// ifs provide periodic boundary conditions
	if (ix<1)
	{
		ixp = N_x;
	}
	if (ix>N_x)
	{
		ixp = 1;
	}
	if (iy<1)
	{
		iyp = N_y;
	}
	if (iy>N_y)
	{
		iyp = 1;
	}
	nc = (iyp - 1)*N_x + ixp;

	return nc;
}

std::vector<int> Index_to_Cell(int N_x, int N_y, int nc) 
{// gets cell index as input and returns cell coordinates
	std::vector<int> ixy;
	int ix;
	int iy;

	iy = (nc-1)/N_x + 1;		// coordinates start at "1" (not "0")
	ix = nc - (iy - 1)*N_x;

	ixy.push_back(ix);
	ixy.push_back(iy);

	return ixy;
}

int Get_Id(std::vector<int> L_R_line, int i)
{// find the receptor "i" in the Ligand list "L_R_line"/get the corresponding index "id"
    int id = 1;
    for(int w = 0, end= L_R_line.size(); w < end; w++)
    {
        if(L_R_line[w + 1] == i)
        {
            return id;
        }
        else
        {
            id += 1;   
        }
    }
    return -1; // if the receptor is not found, return -1
}

std::vector<double> get_position_linked_ligand(std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<double> *X_CM_x, std::vector<double> *X_CM_y, std::vector<double> *X_CM_theta, std::vector<int> *L_to_R, int receptor, double theta_mod)
{/* Returns the position (in the Labframe of reference) of the ligand which is bound to
	* the given (index) receptor "receptor". Uses the given list "L_to_R".
	* If no bound ligand is found, return CM position of the virion
	*/
	std::vector<double> position;
	position.clear();
	int ligand = -1;
	double theta;

	for(int j = 0, end= L_to_R->size(); j < end; j++)
	{
		if(L_to_R->at(j) == receptor)
			ligand = j + 1;
	}

	if(ligand == -1)
	{// if no ligand is found, return position of virion
		position.push_back(X_CM_x->at(X_CM_x->size() - 1));
		position.push_back(X_CM_y->at(X_CM_y->size() - 1));

		return position;
	}
	// calculate the position of the ligand in the Labframe of reference
	theta = X_CM_theta->at((X_CM_theta->size()) - 1) + theta_mod;

	position.push_back((X_l_x->at(ligand-1))*cos(theta) - (X_l_y->at(ligand-1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1));
	position.push_back((X_l_x->at(ligand-1))*sin(theta) + (X_l_y->at(ligand-1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1));

	return position;
}

std::vector<double> get_new_position_linked_ligand(std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<double> r_trial, std::vector<double> *X_CM_theta, std::vector<int> *L_to_R, int receptor)
{/* Returns the new position (in the Labframe of reference) of the ligand which is bound to
	* the given (index) receptor "receptor", meaning that we add the given diffusion vector r_trial
	* to the coords. 
	* If no bound ligand is found, return coords of r_trial and print a warning.
	*/
	std::vector<double> position;
	position.clear();
	int ligand = -1;
	double theta;

	// find the ligand index
	for(int j = 0, end= L_to_R->size(); j < end; j++)
	{
		if(L_to_R->at(j) == receptor)
			ligand = j + 1;
	}

	if(ligand == -1)	// if the ligand is not found
	{
		std::cout << "WARNING get_new_position_linked_ligand" << "\n";
		position.push_back(r_trial[0]);
		position.push_back(r_trial[1]);

		return position;			// return coords of trial vector
	}

	theta = X_CM_theta->at(X_CM_theta->size() - 1);
	// add r_trial vector to get a new position of the ligand (in the labframe of reference
	position.push_back((X_l_x->at(ligand-1))*cos(theta) - (X_l_y->at(ligand-1))*sin(theta) + r_trial[0]);
	position.push_back((X_l_x->at(ligand-1))*sin(theta) + (X_l_y->at(ligand-1))*cos(theta) + r_trial[1]);

	return position;
}

//______________________________________________________________________________________________________________________________
// INITIALIZE
void initialize_sphericalIAV(std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<int> *T_l, std::vector<int> *L_to_R, long *idum, bool clustered)
{/* Initializes a triangular ligand grid on a spherical virium. 
	 * 
	 * T_l denotes type of ligand: 0-> NA, 1->HA
	 *  
	 * if clustered==false: we randomly distribute the NA	
	 * 	with a ratio of NA:HA 1:5.
	 * if clustered== true, we create random patches of ~3 ligands of NA, 
	 * 	alternating through the quadrants.
	 * 
	*/
	double dx= 2* D* sqrt(3);		// horizontal distance between two ligands in the lattice
	double dy= D ;	// the vertical distance of the grid lines
	
	double offset;
	int row=0, i=1;
	double end;
	double xmax2= pow(r_sph,2); 		// square of the maximal x value for a lattice line
	
	while (xmax2 >= 0){		
		if (row==0){
			// the very center
			X_l_x->push_back(0);
			X_l_y->push_back(0);
			// the horizontal center line
			while (i<= r_sph/dx){
				X_l_x->push_back(dx*i);
				X_l_y->push_back(0);
				
				X_l_x->push_back(-dx*i);
				X_l_y->push_back(0);
				i++;
			}
		}
		else{
			// offset for every second line
			if(isodd(row)==true){offset= dx/2, i=0;}
			else{
				offset= 0;
				// place the vertical center line explicitly to prevent placing two ligands at the same position
				i= 1;
				X_l_x->push_back(0);
				X_l_x->push_back(0);
				X_l_y->push_back(row*dy);
				X_l_y->push_back(-row*dy);
			}
			end= (sqrt(xmax2)/dx - offset/dx);

			while (i<= end){
				// upper semicircle
				X_l_x->push_back(dx*i + offset);
				X_l_y->push_back(row*dy);
				X_l_x->push_back(-dx*i - offset);
				X_l_y->push_back(row*dy);
			
				// lower semicircle
				X_l_x->push_back(dx*i + offset);
				X_l_y->push_back(-row*dy);
				X_l_x->push_back(-dx*i - offset);
				X_l_y->push_back(-row*dy);
				i++;
			}
		}
		row++;
		xmax2= pow(r_sph,2) - pow(row*dy,2); 
	}
		// --------------------------------------- PLACE NA
	// needs to be adjusted with size of virion
	int N_l= X_l_y->size();
	int n_NA= int(N_l/6);
	int n_patches= int(n_NA/3);
	int n=0;
	
	// "random" patches
	if (clustered){
		std::vector<double> cx_NA;		// centers of the NA patches
		std::vector<double> cy_NA;		// centers of the NA patches
		std::vector<int> n_c(n_patches, 0);			// number of NA of each patch
		int defined=0;
		// we alternate through the quadrants of the circle and get random patch centers
		for (int patch=0; patch< n_patches; patch++){
			cx_NA.push_back(ran2(idum)	* pow(-1, patch) * (r_sph-D));		// "-D" to avoid patches at the rim
			cy_NA.push_back(ran2(idum)	* pow(-1, patch/2) * (r_sph-D));
		}	
		// assign ligands to patch centers
		for (int l = 0; l < N_l; l++){
			defined=0; 
			if(n < n_NA){
				// iterate through patches 
				for (int i_p=0; i_p< n_patches; i_p++){
					// look if ligand is close to patch center
					if (pow(X_l_x->at(l) - cx_NA.at(i_p) ,2) + pow(X_l_y->at(l) - cy_NA.at(i_p) ,2)  < pow(2*D,2) && n_c.at(i_p) < 3){
						T_l->push_back(0); // NA
						n_c.at(i_p) +=1;
						n++;
						defined=1;
						break;}
				}
			}
			if(not defined){T_l->push_back(1);} // HA
				
			L_to_R->push_back(-1);	// -1= no receptor is bound to ligand
		}		
		// if we did not spread all available NA yet, make patches bigger until we did:
		double widen_radius=0;
		while(n < n_NA){
			widen_radius+= 0.1;
			for (int l = 0; (l < N_l && n < n_NA); l++){
				for (int i_p=0; i_p< n_patches; i_p++){
					// look if ligand is close to patch center
					if (T_l->at(l)!= 0 && pow(X_l_x->at(l) - cx_NA.at(i_p) ,2) + pow(X_l_y->at(l) - cy_NA.at(i_p) ,2)  < pow(2*D,2) + widen_radius){
						T_l->at(l)= 0; // NA
						n++;
						break;}
				}
			}
		}
	}
	if (clustered== false){
		// innit all ligands as HA
		for (int l=0; l< N_l; l++){ 
			T_l->push_back(1); // HA
			L_to_R->push_back(-1);	// -1= no receptor is bound to ligand
		}
		std::cout << "T_l->size(): " << T_l->size() << std::endl;
		std::cout << "N_l: " << N_l << std::endl;
		std::cout << "n_NA: " << n_NA << std::endl;

		// distribute NA
		while(n < n_NA){
			int lig= ran2(idum) * N_l; // gets casted to int between [0, N_l-1]
			if (T_l->at(lig) != 0){T_l->at(lig)= 0; n++;}
			else{;}// skip
		}
	}
	else{
		std::cout << "Please distinguish whether virion is polarized or not. \nExit Program.";
		exit(0);		
	}
	// ---------------------------------------
	// Print and Visualize Grid
    std::ofstream fou;
    fou.open("ligands.xyz");
	fou << N_l << '\n';	
	fou << "A Spherical (2D) virus: \tType (0=NA, 1= HA) \t X[nm] \t Y [nm] \n"; // Header
	for (int j = 0;j < N_l; j++){
		if (T_l->at(j)==1){
			fou << '1' << '\t' << X_l_x->at(j) << '\t' << X_l_y->at(j)  << '\t' << '\n'; }
		else if (T_l->at(j)==0){
			fou << '0' << '\t' << X_l_x->at(j) << '\t' << X_l_y->at(j)  << '\t' << '\n'; }
	}

	std::cout<< "Number of Ligands: " << N_l << '\n';
	std::cout<< "count_NA, n_patches " << n_NA << '\t' << n_patches  << '\n';

	fou.flush();
	fou.close();
}

void initialize_IAV(std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<int> *T_l, std::vector<int> *L_to_R, long *idum)
{/* Initializes a triangular ligand grid on a bacillaform virium.
	 * 9 columns (bottom to top) of ligands on the virium:  (D=λ)
	 * y0 | y1 | y2 | y1 | y2 | y1 | y2 | y1 | y0
	 * y0 on the edges a bit shorter, then a grid of y1 and y2 with offset of sqrt(3)D
	 * x_coords from -D to +D
	 * y coords from -L/2 to + L/2 
	 * 		[y0 from -L/2 + D*sqrt(3) 	to 	L/2 - D*sqrt(3)]
	 * 		[y1 from 			-L/2  			to 		L/2]
	 * 		[y2 from -L/2 - D*sqrt(3) 	to 	L/2 + D*sqrt(3)]
	 * 
	 * x-spacing:	 D
	 * y-spacing:	2 * sqrt(3) * D
	 * 
	 * T_l denotes type of ligand: 0-> NA, 1->HA
	 * Virium polarized so that the lower fifth is populated by NA
	*/
	double b = D*2;
	double d = (sqrt(3)*b)/2;
	std::vector<double> y0;
	std::vector<double> y1;
	std::vector<double> y2;


	
	for (int i = 0; i <= int(L/(2*d))-1; i++)
	{
		y0.push_back(-L/2 + d + i*2*d);
	}

	for (int j = 0; j <= int(L/(2*d)); j++)
	{
		y1.push_back(-L/2 + j*2*d);
	}

	for (int k = 0; k <= int(L/(2*d))+1; k++)
	{
		y2.push_back(-L/2 - d + k*2*d);
	}

	double adjust = abs((y2[0] + y2[y2.size()-1])/2);  // = 0 (??)
    

	int y0_size= y0.size();
	int y1_size= y1.size();
	int y2_size= y2.size();
	
	for (int a = 0; a < 9; a++)
	{
		if(a==0)
		{
			for (int l = 0; l < y0_size; l++)
			{
				X_l_x->push_back(-2*b);
				X_l_y->push_back(y0.at(l) + adjust);
			}
		}
		else if(a==8)
		{
			for (int m = 0; m < y0_size; m++)
			{
				X_l_x->push_back(2*b);
				X_l_y->push_back(y0.at(m) + adjust);
			}
		}
		else if(isodd(a)==true)
		{
			for (int n = 0; n < y1_size; n++)
			{
				X_l_x->push_back((a-4)*b/2);
				X_l_y->push_back(y1.at(n) + adjust);
			}
		}
		else if(isodd(a)==false)
		{
			for (int o = 0; o < y2_size; o++)
			{
				X_l_x->push_back((a-4)*b/2);
				X_l_y->push_back(y2.at(o) + adjust);
			}
		}	
	}

	// Save the resulting Ligandgrid: (compatible for ovito)
	int N_l =X_l_y->size();  // number of ligands
	int count_NA=0;
	
	std::ofstream fou;
    fou.open("ligands.xyz");
	fou << N_l << '\n';	
	fou << "A Bacillaform virus: \tType (0=NA, 1= HA) \t X[nm] \t Y [nm] \n"; // Header
	for (int p = 0; p < N_l; p++)
	{	
		//POLARIZED:		lower sixth is "0"= NA , [-L/2 + L/6= -L/3]
		if(X_l_y->at(p) < (-L/3))	{
			T_l->push_back(0); // NA
			count_NA +=1;}
		else
			{T_l->push_back(1); } // HA
		
		fou << T_l->at(p)  << '\t' << X_l_x->at(p) << '\t' << X_l_y->at(p)  << '\t' <<  '\n'; 	
	    L_to_R->push_back(-1);	// -1= no receptor is bound to ligand
	}
	// Maybe funnel the ligand info into run_info later (using structures??)
	std::cout << "Ligand Grid created.\n";
	std::cout << "Number of Ligands: " << N_l << '\n';
	std::cout << "Number of NA: " << count_NA << '\n';

	fou.flush();
	fou.close();
}

void initialize_SA(std::vector<double> *X_r_x, std::vector<double> *X_r_y, std::vector<int> *T_r, long *m, bool clustered)
{/* Randomly sample the coordinates of N_R receptors and fill them in X_r
	* x-and y-coords.: -SystemSize/2 to + Systemsize/2
	* 
	* m: random seed
	* T_r: assign each receptor a "1" (meaning "active"/"not cleaved")
	*/
	double x;
	double y;
	
 	x = ran2(m)*Taille_Systeme - Taille_Systeme/2;
	y = ran2(m)*Taille_Systeme - Taille_Systeme/2;		// squared system
		
	X_r_x->push_back(x);
    X_r_y->push_back(y);	

	if (clustered==1){
		// Create Receptor Cluster:
		std::ofstream ClusteredReceptors;
		ClusteredReceptors.open("clusteredReceptors.txt");
		
		int N_cl= 6;		// number of clusters
		int r_cl= 100;		// radius of clusters (nm)
		double cl_x, cl_y; 
		
		double roh_bg= 0.2; 	// receptor background density
		double roh_cl= 7 *roh_bg;		// the receptor density in the cluster (with the background receptors adds up to 8-fold the background density)
		int nr_cl= static_cast<int>(roh_cl * M_PI * r_cl*r_cl);		// number of receptors in cluster, 8 times background density * area_of_cluster
		
		// Initialize the receptor background
		unsigned int Rec_bg= N_R - nr_cl*N_cl;	
		while(X_r_x->size() < Rec_bg) {
			x = ran2(m)*Taille_Systeme - Taille_Systeme/2;
			y = ran2(m)*Taille_Systeme - Taille_Systeme/2;

			X_r_x->push_back(x);
			X_r_y->push_back(y);  
		}//	
		
		int sign_x, sign_y;	// used for gaussian distr. positions of the rec.
		
		std::cout << "Cluster Centers:   ";
		for (int row=0; row< 2; row++){
			for (int col=0; col< 3; col++){
					
				// add two cluster centers per column
				cl_x = 	 Taille_Systeme/6 + col * Taille_Systeme/3 - Taille_Systeme/2;
				cl_y = 	 Taille_Systeme/4 + row * Taille_Systeme/2 - Taille_Systeme/2;
			
				std::cout << '[' << cl_x << ','<<cl_y<<"]  ";

				for(int i=0; i< nr_cl; i++){
					// get gaussian distr. around cluster center, maybe rather uniformly?
					if (ran2(m)<0.5) {sign_x= 1;} else {sign_x= -1;}
					if (ran2(m)<0.5) {sign_y= 1;} else {sign_y= -1;}
					
					x = (cl_x + sign_x* RGV(m)/3* r_cl) ;		// -> 3σ= r
					y = (cl_y + sign_y* RGV(m)/3* r_cl);
					
					// periodic boundary conditions/wraparound
					if (x> Taille_Systeme/2){x-= Taille_Systeme;}
					if (x< -Taille_Systeme/2){x+= Taille_Systeme;}
					if (y> Taille_Systeme/2){y-= Taille_Systeme;}
					if (y< -Taille_Systeme/2){y+= Taille_Systeme;}
									
					X_r_x->push_back(x);
					X_r_y->push_back(y);  
					
					ClusteredReceptors << x << '\t' << y << '\t' << '\n'; 	
				}
			}
		}
		
		// Prints
		std::cout << "\nRadius of cluster: " << r_cl	<< " nm"<<	'\n';
		std::cout << "Number of clusters N_cl: "<< N_cl	<<	'\n';
		std::cout << "Number of receptors per cluster nr_cl: "<< nr_cl	<<	'\n';
		std::cout << "Number of Receptor in uniform background:  "<< Rec_bg	<<	'\n';

		//std::cout << "roh_rec:  "<< roh_rec	<<	'\n';

		ClusteredReceptors.flush();
		ClusteredReceptors.close();
	}
	else {
		std::cout << "\nReceptors are not Clustered.\n";
		// if not Clustered just create a receptor field
		while(X_r_x->size() < N_R) {
			x = ran2(m)*Taille_Systeme - Taille_Systeme/2;
			y = ran2(m)*Taille_Systeme - Taille_Systeme/2;

			X_r_x->push_back(x);
			X_r_y->push_back(y);  
		}
	}
	
	for (int j = 0, end= X_r_x->size(); j < end; j++)
	{
		T_r->push_back(1);			// 1: "active"/not cleaved
	}
}

void Initialize_L_R(MATRIX_INT *L_R, std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<double> X_r_x, std::vector<double> X_r_y, std::vector<double> X_CM_x, std::vector<double> X_CM_y, std::vector<double> X_CM_theta, MATRIX_INT Neigh, std::vector<int> hoc_R, std::vector<int> ll_R, std::vector<int> T_l, std::vector<int> L_to_R, int N_x, int N_y, double R_Cell, int N_L, int step)
{/* Initialize L_R, using the neighbour list
	* L_R assigns each Ligand a vector listing the interacting receptors (receptors in range).
	* The functions checks for all neighbouring simulation cells whether a particle is in range.
	*/
    double theta;
    double x, y;
	for (int alpha = 0; alpha < N_L; ++alpha)
	{ 
		// Calculate ligand coordinates in Labframe of ref.:
		theta = X_CM_theta.at((X_CM_theta.size()) - 1);
	    x = (X_l_x->at(alpha))*cos(theta) - (X_l_y->at(alpha))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
		y = (X_l_x->at(alpha))*sin(theta) + (X_l_y->at(alpha))*cos(theta) + X_CM_y[X_CM_y.size() - 1];

		// Get Grid-coords. and cell index of current ligand
		int ix = int((x + Taille_Systeme/2)/R_Cell) + 1;
		int iy = int((y + Taille_Systeme/2)/R_Cell) + 1;
		int index = Cell_to_Index(N_x, N_y, ix, iy);
        
        // iterate through neighbour cells
		for (int j = 0; j < 9; ++j)
		{
			int i = Neigh[index-1][j];
			
		// if neighbour cell is not empty:
			if(hoc_R[i] != 0)
			{
		// check if receptor is in range
	            if(pow(x - X_r_x[hoc_R[i] - 1],2) + pow(y - X_r_y[hoc_R[i] - 1],2) < R_Verlet*R_Verlet)
    				(L_R->at(alpha)).push_back(hoc_R[i]); 
    	// go to next receptor/ligand in the cell
				int g = ll_R[hoc_R[i]];

		// repeat while there is a next receptor
				while(g != 0)
				{
	                if(pow(x - X_r_x[g - 1],2) + pow(y - X_r_y[g - 1],2) < R_Verlet*R_Verlet)
					    (L_R->at(alpha)).push_back(g);
					int k = g;	
					g = ll_R[k];
				}
			}
		}
		
		// If there is a ligand in a bond with its bonded receptor outside the verlet range:
		//	add it as well to L_R
		if(T_l[alpha] == 2 && pow(x - X_r_x[L_to_R[alpha] - 1],2) + pow(y - X_r_y[L_to_R[alpha] - 1],2) > R_Verlet*R_Verlet){
			std::cout << "WARNING: Receptor " << L_to_R[alpha] << " and Ligand " << alpha << " form a bond with excessive length " << pow(x - X_r_x[L_to_R[alpha] - 1],2) + pow(y - X_r_y[L_to_R[alpha] - 1],2)<< '\n'; 
			(L_R->at(alpha)).push_back(L_to_R[alpha]);			
		}
		// insert L_Rsize (=N_L?) at the start of L_Rlist
		int N = (L_R->at(alpha)).size();
		(L_R->at(alpha)).insert((L_R->at(alpha)).begin(), N);
	}
}

void Initialize_R_L(MATRIX_INT *R_L, std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<double> X_r_x, std::vector<double> X_r_y, std::vector<double> X_CM_x, std::vector<double> X_CM_y, std::vector<double> X_CM_theta, MATRIX_INT Neigh, std::vector<int> hoc_L, std::vector<int> ll_L, std::vector<int> T_r, std::vector<int> T_l, std::vector<int> L_to_R, int N_x, int N_y, double R_Cell, int N_L, int step)
{/* Initialize R_L, using the neighbour list
	* R_L assigns each Receptor a vector listing the interacting ligands (ligands in range)
	* The functions checks for all neighbouring simulation cells whether a particle is in range.
	*/
    double theta;
    double x, y;
    double count = 0;

	for (int alpha = 0; alpha < N_R; ++alpha)
	{ 
		// Calculate receptor coordinates in the Labframe of ref.:
		int ix = int((X_r_x[alpha] + Taille_Systeme/2)/R_Cell) + 1;
		int iy = int((X_r_y[alpha] + Taille_Systeme/2)/R_Cell) + 1;
		int index = Cell_to_Index(N_x, N_y, ix, iy);
        // iterate through neighbour cells
		for (int j = 0; j < 9; ++j)
		{
			int i = Neigh[index - 1][j];			
			if(hoc_L[i] != 0) 		// if neighbour cell is not empty:
			{
				// get the ligand (lab frame of reference coordinates)
				theta = X_CM_theta.at((X_CM_theta.size()) - 1);
	    		x = (X_l_x->at(hoc_L[i] - 1))*cos(theta) - (X_l_y->at(hoc_L[i] - 1))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
				y = (X_l_x->at(hoc_L[i] - 1))*sin(theta) + (X_l_y->at(hoc_L[i] - 1))*cos(theta) + X_CM_y[X_CM_y.size() - 1];
				// check whether ligand is in range
	            if(pow(x - X_r_x[alpha],2) + pow(y - X_r_y[alpha],2) < R_Verlet*R_Verlet)
    				{(R_L->at(alpha)).push_back(hoc_L[i]);}
		    					
				int g = ll_L[hoc_L[i]];
				count = 0;
				while(g != 0) 		// repeat while there is a next ligand
				{
					theta = X_CM_theta.at((X_CM_theta.size()) - 1);
				    x = (X_l_x->at(g-1))*cos(theta) - (X_l_y->at(g-1))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
					y = (X_l_x->at(g-1))*sin(theta) + (X_l_y->at(g-1))*cos(theta) + X_CM_y[X_CM_y.size() - 1];

	                if(pow(x - X_r_x[alpha],2) + pow(y - X_r_y[alpha],2) < R_Verlet*R_Verlet)
					    (R_L->at(alpha)).push_back(g);
					int k = g;	
					g = ll_L[k];
					count += 1;
				}
			}
		}	
		// if a receptor is in a bridge, the bonded ligand gets added even if its outside of the range (not sure if this ever happens) 
		if(T_r[alpha] == 2)
		{// search through the ligand
			for(int z = 0; z < N_L; z++)
			{
				if(T_l[z] == 2 && L_to_R[z] == alpha + 1) // I found the bonded ligand
				{
					theta = X_CM_theta.at((X_CM_theta.size()) - 1);
					x = (X_l_x->at(z))*cos(theta) - (X_l_y->at(z))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
					y = (X_l_x->at(z))*sin(theta) + (X_l_y->at(z))*cos(theta) + X_CM_y[X_CM_y.size() - 1];
					// (if it is not already there) add it as well to R_L 
					if(pow(x - X_r_x[alpha],2) + pow(y - X_r_y[alpha],2) > R_Verlet*R_Verlet){
						std::cout << "WARNING: Receptor " << alpha << " and Ligand " << z << " form a bond with excessive length " << pow(x - X_r_x[alpha],2) + pow(y - X_r_y[alpha],2) << '\n'; 
						(R_L->at(alpha)).push_back(z+1);
					}
				}
			}
		}

		int N = (R_L->at(alpha)).size();
		(R_L->at(alpha)).insert((R_L->at(alpha)).begin(), N);
	}
}

void Initialize_Affinity(MATRIX_DOUBLE *a, MATRIX_INT *L_R, std::vector<double> *X_CM_x, std::vector<double> *X_CM_y, std::vector<double> *X_CM_theta, std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<double> *X_r_x, std::vector<double> *X_r_y, std::vector<int> *L_to_R, std::vector<int> T_l, std::vector<int> T_r, double k, int N_L, int t1, int t2)
{/* Intializes the affinities for all ligands for a given interaction t1-t2
	* For each Ligand we fill a vector which denotes either "0" or "k" for 
	* 	each corresponding receptor in Verlet range and which leads with the
	* 	corresponding sum.
	* 
	* a: [[5k, 0, k, ...,k], [0], [3k, 0, k, 0, k ,0 , 0, k], ....]
	* 			|			  |					|	
	* 		|ligand1|	  |ligand2|			|ligand3|		....
	*/
	double a_sum = 0;
	double theta, x, y;

// iterate through ligands
	for (int alpha = 0; alpha < N_L; ++alpha)
	{
		// Get ligands coords. in the Labframe of ref.
		theta = X_CM_theta->at((X_CM_theta->size()) - 1);
	    x = (X_l_x->at(alpha))*cos(theta) - (X_l_y->at(alpha))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
		y = (X_l_x->at(alpha))*sin(theta) + (X_l_y->at(alpha))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);
	
		a_sum = 0;
		
// iterate through receptors (in verlet range)
		for (int i = 1, n_rec=(L_R->at(alpha)).size(); i <n_rec ; ++i)
		{

// check whether in interaction range
			if((pow(x - X_r_x->at(L_R->at(alpha)[i] - 1),2) + pow(y - X_r_y->at(L_R->at(alpha)[i] - 1),2) < D*D))
			{
// check for correct types of ligands/receptors
		        if(T_l[alpha] == t1 && T_r[((L_R->at(alpha)).at(i)) - 1] == t2) 
		        {
	    	       	if(t1 == 2 && t2 == 2) //in case of a bridge
	    	       	{
	    	       	    if(L_to_R->at(alpha) == L_R->at(alpha)[i])	// check whether they are bound to each other 
	    	       	    {
	    	       	        a_sum += k;
			                (a->at(alpha)).push_back(k);
	    	       	    }
	    	       	    else  // both are bound in bridges, but not with each other
	    	       	    {
		                    (a->at(alpha)).push_back(0);
	    	       	    }
	    	       	}
	    	       	else // no bridge
	    	       	{
			            a_sum += k;
			            (a->at(alpha)).push_back(k);
	    	       	}
		        }
		        else // ligand/receptor types do not correspond to the interaction
		        {
		            (a->at(alpha)).push_back(0);
			    }
			}
			else // not in range
			{
				(a->at(alpha)).push_back(0);
			}    
		}
// insert a_sum as first element for each ligand
		(a->at(alpha)).insert((a->at(alpha)).begin(), a_sum);
	}
}

MATRIX_INT Set_Neigh(int N_Cell, int N_x, int N_y) 
{// creates matrix where every cell gets its 8 neighbour indeces assigned (Neighbour-List)
	int NNeigh = 9;
	std::vector<int> ixy(2);

	MATRIX_INT Neigh;
	Neigh.resize(N_Cell, std::vector<int>(NNeigh));

	for (int i = 0; i < N_Cell; ++i)
	{
		ixy = Index_to_Cell(N_x, N_y, i+1);
		int ix = ixy[0];
		int iy = ixy[1];
		int cnt = 0;
		// iterate over neighbour coords.
		for (int icx = ix-1; icx < ix + 2; ++icx)
		{
			for (int icy = iy-1 ; icy < iy + 2; ++icy)
			{// add neighbour indeces to matrix
				int index = Cell_to_Index(N_x,N_y,icx,icy);
				Neigh[i][cnt] = index;
				cnt=cnt+1;
			}
		}
	}
	return Neigh;
}

std::vector<int> Get_bridges(std::vector<int> *T_l, std::vector<int> *L_to_R)
{/* returns "list_bridged_receptors",
	*  a flat list with all receptors (their index) that are forming currently a bridge
	*/
	std::vector<int> list_bridged_receptors;

	for(int i = 0, end= T_l->size(); i<end; i++){
		if(T_l->at(i) == 2)
			list_bridged_receptors.push_back(L_to_R->at(i));
	}

	return list_bridged_receptors;
}


MATRIX_INT Init_CellLists(double R_Cell, int N_Cell, int N_part, int N_x, int N_y, std::vector<double> *X_part_x, std::vector<double> *X_part_y)
{/* Creates the lists hoc, ll, bl
	* hoc[i]: (last added) Ligand/Receptor in the Cell i
	* ll[j]: points to the Ligand/Receptor which is in the same cell as j and was added before
	* 			0-> there are no additional L/R in the same cell/ j is the first one which got added
	* bl: same list as ll, but pointing in the opposite direction
	* 		i.e.: ll[j]=k, bl[k]=j
	*/
	std::vector<int> hoc(N_Cell + 1);
	std::vector<int> ll(N_part + 1);
	std::vector<int> bl(N_part + 1);
	MATRIX_INT v;

	for (int icell = 1; icell < N_Cell + 1; ++icell)
	{
		hoc[icell] = 0;
	}
	for (int i = 1; i < N_part + 1; ++i)
	{
		int ix = int((X_part_x->at(i-1) + Taille_Systeme/2)/R_Cell) + 1;
		int iy = int((X_part_y->at(i-1) + Taille_Systeme/2)/R_Cell) + 1;
		int nc = Cell_to_Index(N_x,N_y,ix,iy);

		ll[i] = hoc[nc];
		if (hoc[nc] != 0)
		{
			bl[hoc[nc]] = i;
		}
		hoc[nc] = i;
	}
	
	v.push_back(ll);
	v.push_back(bl);
	v.push_back(hoc);
	
	return v;
}

// ----------------------------------- AFFINITIES: ------------------------------------------------
std::vector<double> Get_a_X_tot(MATRIX_DOUBLE *a_X, int N_L)  
{/* Creating a flat vector a which lists for each ligand the total affinity/sum
	* and leads with the total sum of these sums "a_X_tot"
	*/
	std::vector<double> a_X_tot;
	a_X_tot.clear();
	double a_X_sum = 0.0;

	for (int i = 0; i < N_L; ++i)
	{
		a_X_tot.push_back((a_X->at(i)).at(0));	// takes a_sum of each ligand
		a_X_sum += (a_X->at(i)).at(0);
	}
	a_X_tot.insert(a_X_tot.begin(), a_X_sum);	// a_X_tot leads with total sum of ligand affinities
	return a_X_tot;
}

std::vector<double> Get_a_tot(std::vector<double> *a_d_tot, std::vector<double> *a_on_tot, std::vector<double> *a_off_tot, int N_L) 
{/* Puts the total affinities in a vector "a", which is used in order to choose a reaction
	* a: [a_tot, a_d, a_on, a_off]
	*/
	double a_tot;
	double a_d = 0.0;
	double a_on = 0.0;
	double a_off = 0.0;
	double epsilon = 0.0000001;
	std::vector<double> a;
	a.clear();

	for (int i = 0; i < N_L; ++i)
	{
		a_d += a_d_tot->at(i+1);
		if(a_d < epsilon)
		{
		    a_d = 0.0;
		}
	}

	for (int j = 0; j < N_L; ++j)
	{
		a_on += a_on_tot->at(j+1);
		if(a_on < epsilon)
		{
		    a_on = 0.0;
		}
	}

	for (int k = 0; k < N_L; ++k)
	{
		a_off += a_off_tot->at(k+1);
		if(a_off < epsilon)
		{
		    a_off = 0.0;
		}
	}

	a_tot = a_d + a_on + a_off;
	a.push_back(a_tot);
	a.push_back(a_d);
	a.push_back(a_on);
	a.push_back(a_off);

	return a;
}

// -----------------------------------------------------------------------------------------------------------------------------
// CHECK BREAKAGE
bool check_breakage(std::vector<double> *X_CM_x, std::vector<double> *X_CM_y, std::vector<double> *X_CM_theta, std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<double> *X_r_x, std::vector<double> *X_r_y, std::vector<int> *T_l, std::vector<int> *L_to_R, double DX, double DY, double DTHETA, bool verbose= true)
{/* Verifies that the diffusion does not break a bridge (looks at all bridges) .
 * Gets all the current Ligand and Receptor positions as well as the diffusion
 * steps DX, DY and rotation DTheta
 * Returns 	bk= true	-> Bond breakage
 * 				bk= false	-> no bond breakage	
*/
	bool bk = false;
	int receptor = 0;
    double theta;
    double x;
    double y;

    theta = DTHETA + X_CM_theta->at((X_CM_theta->size()) - 1);

	for (int i = 0, end= X_l_x->size(); i < end; ++i)		// iterate through ligands
	{
		if (T_l->at(i) == 2)					// those, that form bridges
		{
			receptor = L_to_R->at(i);
		    
		    x = (X_l_x->at(i))*cos(theta) - (X_l_y->at(i))*sin(theta) + X_CM_x->at((X_CM_x->size()) -1) ;
		    y = (X_l_x->at(i))*sin(theta) + (X_l_y->at(i))*cos(theta) + X_CM_y->at((X_CM_y->size()) -1) ;
		    
		    // check the lengths of the bond
			if (pow(x + DX - X_r_x->at(receptor - 1),2) + pow(y + DY - X_r_y->at(receptor - 1),2) > D*D)
			{
				bk = true;
				if (verbose){
				std::cout << "Breaking of a normal bridge. Squared bond length: " << pow(x + DX - X_r_x->at(receptor - 1),2) + pow(y + DY - X_r_y->at(receptor - 1),2) << '\n';
				std::cout << "Lig " << i<<  ", X_lig:  " << x << '\t' << y << '\n';	
				std::cout << "Rec " << receptor -1 << ", X_rec:  " << X_r_x->at(receptor - 1) << '\t' << X_r_y->at(receptor - 1) << '\n';	
			//	std::cout << "Rec-1 " << receptor -2 << ", X_rec:  " << X_r_x->at(receptor - 2) << '\t' << X_r_y->at(receptor - 2) << '\n';	
			//	std::cout << "Rec+1 " << receptor  << ", X_rec:  " << X_r_x->at(receptor ) << '\t' << X_r_y->at(receptor ) << '\n';	
				
				}
			}
		}
	}
	
	return bk;
}


bool check_breakage_receptors(std::vector<double> *X_CM_x, std::vector<double> *X_CM_y, std::vector<double> *X_CM_theta, std::vector<double> *X_l_x, std::vector<double> *X_l_y, double r_x, double r_y, std::vector<int> *T_l, std::vector<int> *L_to_R, std::vector<int> *list_bridged_receptors, double DX, double DY, int rec)
{//	
	bool bk_rec;
	int receptor = 0;
    double theta;
    double x;
    double y;
    int i = -1;
    
    theta = X_CM_theta->at((X_CM_theta->size()) - 1);
	//std::cout << "overextended bridges : ";

	while(receptor != rec + 1)
	{
		if (T_l->at(i + 1) == 2)
		{
			receptor = L_to_R->at(i + 1);
			//std::cout << receptor << '\n';
		}
	
		i += 1;

	}	

    x = (X_l_x->at(i))*cos(theta) - (X_l_y->at(i))*sin(theta) + X_CM_x->at((X_CM_x->size()) -1) ;
    y = (X_l_x->at(i))*sin(theta) + (X_l_y->at(i))*cos(theta) + X_CM_y->at((X_CM_y->size()) -1) ;
    
	if (pow(x - (r_x + DX),2) + pow(y - (r_y + DY),2) > D*D)
	{
		bk_rec = true;
	}

	return bk_rec;
}

// -----------------------------------------------------------------------------------------------------------------------------
// REACTIONS:
int Choose(std::vector<double> V, long *g) 
{/* Returns an int i which defines either the type of reaction (1-3) or the index of the reacting ligand/receptor 
	* (depending on what we choose as an input 'V':
	* - Using the total affinites "V= a_tot", we choose a reaction with the random function "ran2".
	* 						i: 1= destroy SA, 2= forming a bridge, 3= breaking a bridge
	* 
	* - Using the total affinities of the ligands (i.e. "V= a_d_tot"),
	* 	 we choose the index i of the reacting ligand
	* - Using the affinities of a specific ligand (i.e. "V= a_d.at(alpha)"),
	*  	 we choose the index of the reacting receptor i (with ligand "alpha")
	* 
	*/
	double r = ran2(g); 
	double sumpar = 0.0;
	
    if(V[0] == 0)   {return -1;} // V[0]= sum of affinities=0 -> no reaction 
        
  	int i = 0;							
    int V_size=V.size();
	while(sumpar <= r && i<V_size) 
        { 
		i += 1;		// we skip V[0] since its the sum
		sumpar += (V[i]/V[0]);	
	}
	
	return i;
}

void Make_reaction(MATRIX_DOUBLE *a_d, MATRIX_DOUBLE *a_on, MATRIX_DOUBLE *a_off, std::vector<double> a_d_tot, std::vector<double> a_on_tot, std::vector<double> a_off_tot, MATRIX_INT L_R, MATRIX_INT R_L, std::vector<double> *X_CM_x, std::vector<double> *X_CM_y, std::vector<double> *X_CM_theta, std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<double> *X_r_x, std::vector<double> *X_r_y, std::vector<int> *T_l, std::vector<int> *T_r, std::vector<int> *L_to_R, std::vector<int> *list_bridged_receptors, std::vector<int> *list_CB, int r, long *idum, int step, std::ofstream *depleted)
{// Perform the reaction chosen by the input 'r'
	//choisit le ligand et le rÃ©cepteur et effectue la rÃ©action choisie auparavant 
    double a_sum;
    double theta, x, y;
    std::vector<int> new_list_CB;

	if (r == 1) //destruction SA
	{
		int lig = Choose(a_d_tot, idum); 			// choose reacting ligand (total index)
		int ir = Choose((a_d->at(lig-1)), idum); 	// choose reacting receptor (index in ligand list)
		int rec = L_R[lig-1][ir]; 					// get the receptors real index 

		*depleted << step << "," << rec << "\n"; 		// print the depleted receptor, numbered from 1 to N_R (!!should rather change that)

		// Update affinities: remove k_d
		(a_d->at(lig-1))[0] -= (a_d->at(lig-1))[ir]; // subtract the affinity of this reaction from the sum
		(a_d->at(lig-1))[ir] = 0;					// remove affinity to the same reaction again

		// remove the receptor from the affinities of all ligands in verlet range
		for (int j = 1, end= R_L[rec-1].size(); j < end; j++) 
		{
			// get Ligand Labframe coords.
			theta = X_CM_theta->at((X_CM_theta->size()) - 1);
		    x = (X_l_x->at(R_L[rec-1][j] - 1))*cos(theta) - (X_l_y->at(R_L[rec-1][j] - 1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
			y = (X_l_x->at(R_L[rec-1][j] - 1))*sin(theta) + (X_l_y->at(R_L[rec-1][j] - 1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);

			if((pow(x - X_r_x->at(rec-1),2) + pow(y - X_r_y->at(rec-1),2) < D*D)) // if Ligand in interaction range of receptor i
			{
				int ir2 = Get_Id(L_R[(R_L[rec-1][j]) - 1], rec);				// index of the receptor i in the R_L of the current ligand 

			    if(T_l->at(R_L[rec-1][j] - 1) == 0)	// if Ligand is NA, deactivate receptor for a_d
				{
					(a_d->at(R_L[rec-1][j] - 1))[0] -= (a_d->at(R_L[rec-1][j] - 1))[ir2];
				    (a_d->at(R_L[rec-1][j] - 1))[ir2] = 0.0;
				}
				if(T_l->at(R_L[rec-1][j] - 1) == 1)	// if Ligand is HA, deactivate receptor for a_on
				{
					(a_on->at(R_L[rec-1][j] - 1))[0] -= (a_on->at(R_L[rec-1][j] - 1))[ir2];
					(a_on->at(R_L[rec-1][j] - 1))[ir2] = 0.0;
				}
			}
		}
	    T_r->at(rec-1) = 0;								// adjust type of receptor to "off"
	}


	else if(r == 2) // construct a bridge
	{
		int lig = Choose(a_on_tot, idum); 			// choose reacting ligand (total index)
		int ir = Choose((a_on->at(lig - 1)), idum); 	// choose reacting receptor (index in ligand list)
		int rec = L_R[lig - 1][ir];						// get the receptors real index 

		// Update affinities: add k_off
		(a_off->at(lig - 1))[0] += k_off;
		(a_off->at(lig - 1))[ir] = k_off; 
		
		// adjust the affinities of all ligands in verlet range	
		for (int k = 1, end=R_L[rec-1].size(); k < end; k++)
		{
			theta = X_CM_theta->at((X_CM_theta->size()) - 1);
		    x = (X_l_x->at(R_L[rec-1][k] - 1))*cos(theta) - (X_l_y->at(R_L[rec-1][k] - 1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
			y = (X_l_x->at(R_L[rec-1][k] - 1))*sin(theta) + (X_l_y->at(R_L[rec-1][k] - 1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);
			// check if in interaction range
			if((pow(x - X_r_x->at(rec-1),2) + pow(y - X_r_y->at(rec-1),2) < D*D))
			{
				int ir2 = Get_Id(L_R[(R_L[rec-1][k]) - 1], rec);
							    
				if (T_l->at(R_L[rec-1][k] - 1) == 1)		 	// if Ligand is HA, deactivate receptor for a_on
				{
					(a_on->at(R_L[rec-1][k] - 1))[0] -= (a_on->at(R_L[rec-1][k] - 1))[ir2];
					(a_on->at(R_L[rec-1][k] - 1))[ir2] = 0.0;				
				}
				else if (T_l->at(R_L[rec-1][k] - 1) == 0)		// if Ligand is NA, deactivate receptor for a_d
				{
					(a_d->at(R_L[rec-1][k] - 1))[0] -= (a_d->at(R_L[rec-1][k] - 1))[ir2];
					(a_d->at(R_L[rec-1][k] - 1))[ir2] = 0.0;						
				}
			}
		}
		// adjust affinities between bound ligand and all other receptors -> 0.0
		for (int k = 0, end= (a_on->at(lig-1)).size(); k < end; k++)
		{
		    a_on->at(lig-1).at(k) = 0.0;
		}
		// list them as bridges
		T_l->at(lig - 1) = 2;								// bridge-type
		T_r->at(rec - 1) = 2;								// bridge-type
		L_to_R->at(lig - 1) = rec;	

		// update list_bridged_receptors and list_CB
		list_bridged_receptors->push_back(rec);
	}


	else if(r == 3) // destroy a bridge
	{
		int lig = Choose(a_off_tot, idum); 			// choose reacting ligand (total index)
		int ir = Choose((a_off->at(lig-1)), idum); 	// choose reacting receptor (index in ligand list)
		int rec = L_R[lig-1][ir]; 						// get the receptors real index 
		
		// Update
		(a_off->at(lig-1))[0] = 0.0;
		(a_off->at(lig-1))[ir] = 0.0; 
		// Update types
		T_l->at(lig-1) = 1;
		T_r->at(rec - 1) = 1;

		// adjust the affinities of all ligands in verlet range	of receptor "i"
		for (int m = 1, end= R_L[rec-1].size(); m < end; m++)
		{

			theta = X_CM_theta->at((X_CM_theta->size()) - 1);
		    x = (X_l_x->at(R_L[rec-1][m] - 1))*cos(theta) - (X_l_y->at(R_L[rec-1][m] - 1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
			y = (X_l_x->at(R_L[rec-1][m] - 1))*sin(theta) + (X_l_y->at(R_L[rec-1][m] - 1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);
			// check if in interaction range
			if((pow(x - X_r_x->at(rec-1),2) + pow(y - X_r_y->at(rec-1),2) < D*D))
			{
				int ir2 = Get_Id(L_R[(R_L[rec-1][m]) - 1], rec);
	
				// if ligand HA and receptor not the current one: increase a_on
				if (T_l->at(R_L[rec-1][m] - 1) == 1 && R_L[rec-1][m] != lig)	
				{
					(a_on->at(R_L[rec-1][m] - 1))[0] += k_on;
					(a_on->at(R_L[rec-1][m] - 1))[ir2] = k_on;				
				}
				else if (T_l->at(R_L[rec-1][m] - 1) == 0)	// if ligand NA, increase a_d
				{
					(a_d->at(R_L[rec-1][m] - 1))[0] += k_d;
					(a_d->at(R_L[rec-1][m] - 1))[ir2] = k_d;						
				}
			}
		}
		
		// adjust the affinities of all receptors (which are on "1") in verlet range of ligand "alpha"
		a_sum = 0;
		for (int k = 1, end= L_R[lig-1].size(); k < end; k++)
		{
			theta = X_CM_theta->at((X_CM_theta->size()) - 1);
		    x = (X_l_x->at(lig-1))*cos(theta) - (X_l_y->at(lig-1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
			y = (X_l_x->at(lig-1))*sin(theta) + (X_l_y->at(lig-1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);
		
		    if ((T_r->at(L_R[lig-1][k] - 1) == 1)  && (pow(x - X_r_x->at(L_R[lig-1][k] - 1),2) + pow(y - X_r_y->at(L_R[lig-1][k] - 1),2) < D*D))
		    {
		        a_on->at(lig-1)[k] += k_on;
		        a_sum += k_on;
		    }
		    else
		    {
		        a_on->at(lig-1)[k] = 0.0;
		    }
		}
		a_on->at(lig-1)[0] = a_sum; 
		
		L_to_R->at(lig - 1) = -1;		// remove bridge from L_to_R


		//update list_bridged_receptors
		int idb = Get_Id(*list_bridged_receptors, rec);

		// if rec is at first index:
		if(idb == -1 && list_bridged_receptors->at(0) == rec)
			idb = 0;
		int size_rec_list= list_bridged_receptors->size();
		if(idb == size_rec_list && list_bridged_receptors->at(0) == rec)	// some wraparound??
			idb = 0;
		if(rec != list_bridged_receptors->at(idb))
			std::cout << "WARNING MKREACT" << "\n";
		list_bridged_receptors->erase(list_bridged_receptors->begin()+idb); // remove the bridge from index
	}    
	else
	{
	    std::cout << "error in choice of reaction: no reaction has been chosen" << '\n';
	}
}


//-----------------------------------------------------------------------------------------------------------------------------
// DIFFUSION

void diffusion(std::vector<double> *X_CM_x, std::vector<double> *X_CM_y, std::vector<double> *X_CM_theta, std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<double> *X_r_x, std::vector<double> *X_r_y, std::vector<int> *T_l, std::vector<int> *L_to_R, long *idum, double delta_t, double D_II, double D_ortho, double D_rot)
{// diffuse the Virus and the connected receptors
    std::vector<double> random_numbers;

    for(int k = 0; k<3; k++)
		{random_numbers.push_back(RGV(idum));}
    
    double old_theta = X_CM_theta->at((X_CM_theta->size()) - 1);
    
	double DX_Brwn = (-1)*sqrt(2*D_II*delta_t)*sin(old_theta)*random_numbers[0] + sqrt(2*D_ortho*delta_t)*cos(old_theta)*random_numbers[1];
	double DY_Brwn = sqrt(2*D_II*delta_t)*cos(old_theta)*random_numbers[0] + sqrt(2*D_ortho*delta_t)*sin(old_theta)*random_numbers[1]; 
	double DTHETA_Brwn = sqrt(2*D_r*delta_t)*random_numbers[2];
	
	// diffuse the Virion
	X_CM_x->push_back(X_CM_x->at((X_CM_x->size()) -1) + DX_Brwn);
	X_CM_y->push_back(X_CM_y->at((X_CM_y->size()) -1) + DY_Brwn);
	double new_theta= old_theta + DTHETA_Brwn;
	X_CM_theta->push_back(new_theta);
	
	// diffuse the connected receptors
	for (int lig = 0; lig < X_l_x->size(); ++lig){
		if (T_l->at(lig) == 2){// get lig-receptor pair
			int rec = L_to_R->at(lig) -1;
		    
		    // get ligand translation
		    double dx = (X_l_x->at(lig))*cos(new_theta) - (X_l_y->at(lig))*sin(new_theta) - 
			(X_l_x->at(lig))*cos(old_theta) + (X_l_y->at(lig))*sin(old_theta) + DX_Brwn;
		    double dy = (X_l_x->at(lig))*sin(new_theta) + (X_l_y->at(lig))*cos(new_theta) -
		     (X_l_x->at(lig))*sin(old_theta) - (X_l_y->at(lig))*cos(old_theta) + DY_Brwn;
		    
		    // apply the same translation to the bound receptor
			X_r_x->at(rec) = X_r_x->at(rec) + dx;
	       	X_r_y->at(rec) = X_r_y->at(rec) + dy;
	       	
			// periodic boundary conditions:
	       	double size = Taille_Systeme;
			if(X_r_x->at(rec) > size/2)
				X_r_x->at(rec) = X_r_x->at(rec) - size;
			if(X_r_x->at(rec) < -1*size/2)
				X_r_x->at(rec) = X_r_x->at(rec) + size;

			if(X_r_y->at(rec) > size/2)
				X_r_y->at(rec) = X_r_y->at(rec) - size;
			if(X_r_y->at(rec) < -1*size/2)
				X_r_y->at(rec) = X_r_y->at(rec) + size;
		}
	}
}

// DIFFUSE BOUND RECEPTORS
void diffusion_bound_receptors(std::vector<double> *X_CM_x, std::vector<double> *X_CM_y, std::vector<double> *X_CM_theta, std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<double> *X_r_x, std::vector<double> *X_r_y, std::vector<int> *T_l, std::vector<int> *L_to_R, std::vector<int> *list_bridged_receptors, long *idum, double delta_t)
{/* Moves all the bound receptor by a brownian motion:
 *  D_Brwn = sqrt(2*D_receptor*delta_t)*random_number
 * 
 * Bridged receptors which would break because of the diffusion do not 
 * move.
 */
    std::vector<double> random_numbers;
    double DX_Brwn, DY_Brwn;
    
    for(int k = 0; k<2; k++)
		{random_numbers.push_back(0.0);}
		
    // cycle through the bridged receptors
    for(int j = 0, end= list_bridged_receptors->size(); j < end; j++){
		int rec= list_bridged_receptors->at(j) - 1;
    	
		// move bridged receptors
		int trial=0;
		bool bk_rec = true;
		while (bk_rec ==true){ // try until we find diffusion that does not break bonds
			if (trial== 100)
			{std::cout << "WARNING!\n No suitable receptor diffusion after 100 trials for rec " << rec << ".\n Last diffusion vector: [" << DX_Brwn << ',' << DY_Brwn << "]" << std::endl; exit(0);}
			
			random_numbers[0] = RGV(idum);
			random_numbers[1] = RGV(idum);

			DX_Brwn = sqrt(2*D_receptor*delta_t)*random_numbers[0];
			DY_Brwn = sqrt(2*D_receptor*delta_t)*random_numbers[1];

			bk_rec = check_breakage_receptors(X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x->at(rec), X_r_y->at(rec), T_l, L_to_R, list_bridged_receptors, DX_Brwn, DY_Brwn, rec);
			++trial;
			}
		total_tests+= trial;
		rejections+= trial -1;
		
		// Diffuse the receptor
		X_r_x->at(rec) = X_r_x->at(rec) + DX_Brwn;
		X_r_y->at(rec) = X_r_y->at(rec) + DY_Brwn;
		 
		// periodic boundary conditions/wraparound
		double size = Taille_Systeme;
		if(X_r_x->at(rec) > size/2)
			X_r_x->at(rec) = X_r_x->at(rec) - size;
		if(X_r_x->at(rec) < -1*size/2)
			X_r_x->at(rec) = X_r_x->at(rec) + size;
		// for y
		if(X_r_y->at(rec) > size/2)
			X_r_y->at(rec) = X_r_y->at(rec) - size;
		if(X_r_y->at(rec) < -1*size/2)
			X_r_y->at(rec) = X_r_y->at(rec) + size;
	}
}

// DIFFUSE UNBOUND RECEPTORS
void diffusion_loose_receptors(std::vector<double> *X_r_x, std::vector<double> *X_r_y, std::vector<int> *list_bridged_receptors, long *idum, double delta_t)
{/* Moves the receptor which are not bound by a brownian motion:
 *  D_Brwn = sqrt(2*D_receptor*delta_t)*random_number
 */
    std::vector<double> random_numbers;
    
    for(int k = 0; k<2; k++)
		{random_numbers.push_back(0.0);}
	// cycle thorugh all receptors
    for (int rec = 0; rec < N_R; ++rec){
		// find bridged receptors
    	bool bridged = false;
    	for(int j = 0, end= list_bridged_receptors->size(); j < end; j++){
    		if(list_bridged_receptors->at(j) - 1 == rec)
				{bridged = true;}
    	}
		// skip bridged receptors
    	if(bridged == true){;}
		// diffuse the rest
		else{
			random_numbers[0] = RGV(idum);
			random_numbers[1] = RGV(idum);
			double DX_Brwn = sqrt(2*D_receptor*delta_t)*random_numbers[0];
			double DY_Brwn = sqrt(2*D_receptor*delta_t)*random_numbers[1];
	
			// Diffuse the receptor
			X_r_x->at(rec) = X_r_x->at(rec) + DX_Brwn;
			X_r_y->at(rec) = X_r_y->at(rec) + DY_Brwn;
			 
			// periodic boundary conditions/wraparound
			double size = Taille_Systeme;
			if(X_r_x->at(rec) > size/2)
				X_r_x->at(rec) = X_r_x->at(rec) - size;
			if(X_r_x->at(rec) < -1*size/2)
				X_r_x->at(rec) = X_r_x->at(rec) + size;
			// for y
			if(X_r_y->at(rec) > size/2)
				X_r_y->at(rec) = X_r_y->at(rec) - size;
			if(X_r_y->at(rec) < -1*size/2)
				X_r_y->at(rec) = X_r_y->at(rec) + size;
		}
	}
}


// -----------------------------------------------------------------------------------------------------------------------------
// UPDATE
void UpdateCell_Lists(int Old_Cell, int New_Cell, int ipar, std::vector<int> *ll, std::vector<int> *bl, std::vector<int> *hoc)
{/* Update the simulation cell lists for the particle "ipar" which changed its cell
	* ipar: index of the current particle
	* hoc is the last particle in the cell, ll points at the index of the previous particle until 0, bl points at the next particle until 0.
	* Can be used for either ll_L or ll_R, depending on which get provided in the function arguments.
	*/
	
	// REMOVE ipar out of the Old_Cell
	if (ll->at(ipar) != 0)
	{
		if (bl->at(ipar) != 0) // ipar is somewhere in the MIDDLE of the cell (index-wise)
		{
			ll->at(bl->at(ipar)) = ll->at(ipar);
			bl->at(ll->at(ipar)) = bl->at(ipar);
		}
		else // ipar is the FIRST particle in the cell (the last one added)
		{
			bl->at(ll->at(ipar)) = 0;
			hoc->at(Old_Cell) = ll->at(ipar);
		}
	}
	else
	{
		if (bl->at(ipar) != 0) // ipar is the LAST particle in cell (linked list)
		{
			ll->at(bl->at(ipar)) = 0;
		}
		else // ipar is the ONLY particle in the cell
		{
			hoc->at(Old_Cell) = 0;
		}
	}


	// INSERT ipar in the NEW_Cell (as last added)
	bl->at(ipar) = 0;
	ll->at(ipar) = hoc->at(New_Cell);
	if (hoc->at(New_Cell) !=0)
	{
		bl->at(hoc->at(New_Cell)) = ipar;
	}
	hoc->at(New_Cell) = ipar;
}

void Update_Cell_L(std::vector<int> *ll_L, std::vector<int> *bl_L, std::vector<int> *hoc_L, std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<double> *X_CM_x, std::vector<double> *X_CM_y, std::vector<double> *X_CM_theta, int N_L, int N_x, int N_y, double R_Cell, double X_CM_x_init, double X_CM_y_init, double X_CM_theta_init)
{/* Update the lists denoting the ligands in each simulation cell
	* Checks for each Ligand whether it left its Sim. cell. and,
	* If so, updates then ll, bl and hoc, using Update_CellList
	*/
	double theta = X_CM_theta->at(X_CM_theta->size() - 1);
	double otheta = X_CM_theta_init;

	for (int i = 1; i < N_L + 1; i++)
	{
				// ligand coord. in virus framce of ref.
				double x = (X_l_x->at(i-1))*cos(theta) - (X_l_y->at(i-1))*sin(theta);
				double y = (X_l_x->at(i-1))*sin(theta) + (X_l_y->at(i-1))*cos(theta);
				// get Labframe coord. of the ligands
				double x_p = x + X_CM_x->at(X_CM_x->size() - 1);
				double y_p = y + X_CM_y->at(X_CM_y->size() - 1);
				
				// get index of the cooresponding simulation cell
				int ix = int((x_p + Taille_Systeme/2)/R_Cell) + 1;
				int iy = int((y_p + Taille_Systeme/2)/R_Cell) + 1;
				int nc = Cell_to_Index(N_x,N_y,ix,iy);

				// the same for the old virus position
				double ox = (X_l_x->at(i-1))*cos(otheta) - (X_l_y->at(i-1))*sin(otheta);
				double oy = (X_l_x->at(i-1))*sin(otheta) + (X_l_y->at(i-1))*cos(otheta);
				double ox_p = ox + X_CM_x_init;
				double oy_p = oy + X_CM_y_init;
				// old sim. Cell
				int oix = int((ox_p + Taille_Systeme/2)/R_Cell) + 1;
				int oiy = int((oy_p + Taille_Systeme/2)/R_Cell) + 1;
				int oc = Cell_to_Index(N_x,N_y,oix,oiy);

				if(oc != nc)		// Diffusion leads to ligands changing their simulation cell
				{
					UpdateCell_Lists(oc, nc, i, ll_L, bl_L, hoc_L);
				}
	}
}

void Update_Cell_R(std::vector<int> *ll_R, std::vector<int> *bl_R, std::vector<int> *hoc_R, std::vector<double> *X_r_x, std::vector<double> *X_r_y, std::vector<double> *x_init_R, std::vector<double> *y_init_R, int N_x, int N_y, double R_Cell)
{/*Update the lists denoting the receptor in each simulation cell
	* Checks for each Receptor whether it left its Sim. cell. and,
	* If so, updates then ll, bl and hoc, using Update_CellList
	* [one should ll,bl, hoc]
	*/
	for (int i = 1; i < N_R + 1; i++)
	{
				// Get Receptor positions
				double x = (X_r_x->at(i-1));
				double y = (X_r_y->at(i-1));
				double x_p = x;
				double y_p = y;
				// get index of the cooresponding simulation cell
				int ix = int((x_p + Taille_Systeme/2)/R_Cell) + 1;
				int iy = int((y_p + Taille_Systeme/2)/R_Cell) + 1;
				int nc = Cell_to_Index(N_x,N_y,ix,iy);
				
				// the same for the old position
				double ox_p = x_init_R->at(i-1);
				double oy_p = y_init_R->at(i-1);
				// old sim. Cell
				int oix = int((ox_p + Taille_Systeme/2)/R_Cell) + 1;
				int oiy = int((oy_p + Taille_Systeme/2)/R_Cell) + 1;
				int oc = Cell_to_Index(N_x,N_y,oix,oiy);

				if(oc != nc)			// check if particles left their sim. cell
				{
					UpdateCell_Lists(oc, nc, i, ll_R, bl_R, hoc_R);
				}	
	}

}

void Reinitialize_Lists(MATRIX_INT *L_R, MATRIX_INT *R_L, MATRIX_DOUBLE *a_d, MATRIX_DOUBLE *a_on, MATRIX_DOUBLE *a_off, std::vector<double> *a_d_tot, std::vector<double> *a_on_tot, std::vector<double> *a_off_tot, std::vector<double> *affinity, MATRIX_INT Neigh, std::vector<int> hoc_R, std::vector<int> hoc_L, std::vector<int> ll_L, std::vector<int> ll_R, std::vector<double> X_l_x, std::vector<double> X_l_y, std::vector<double> X_r_x, std::vector<double> X_r_y, std::vector<double> X_CM_x, std::vector<double> X_CM_y, std::vector<double> X_CM_theta, int N_x, int N_y, double R_Cell, int N_L, std::vector<int> T_l, std::vector<int> T_r, std::vector<int> *L_to_R, int step)
{/* Update L_R, R_L and the affinities a_d, a_on, a_off, by emptying/clearing them and
	* reinitilizing all again.
	*/
	// CLEAR the lists
    for(int i=0; i<N_L; i++)
    {
        (L_R->at(i)).clear();
        (a_d->at(i)).clear();
        (a_on->at(i)).clear();
        (a_off->at(i)).clear();
    }
    for(int i=0; i<N_R; i++)
    {
        (R_L->at(i)).clear();
    }
    a_d_tot->clear();
    a_on_tot->clear();
    a_off_tot->clear();
    affinity->clear();
    
    // Initialize everything again 
	Initialize_L_R(L_R, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_R, ll_R, T_l, *L_to_R, N_x, N_y, R_Cell, N_L, step);
	Initialize_R_L(R_L, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_L, ll_L, T_r, T_l, *L_to_R, N_x, N_y, R_Cell, N_L, step);
	Initialize_Affinity(a_d, L_R, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, L_to_R, T_l, T_r, k_d, N_L, 0, 1);
	Initialize_Affinity(a_on, L_R, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, L_to_R, T_l, T_r, k_on, N_L, 1, 1);
    Initialize_Affinity(a_off, L_R, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, L_to_R, T_l, T_r, k_off, N_L, 2, 2);
    *a_d_tot = Get_a_X_tot(a_d, N_L);
	*a_on_tot = Get_a_X_tot(a_on, N_L);
	*a_off_tot = Get_a_X_tot(a_off, N_L);
	*affinity = Get_a_tot(a_d_tot, a_on_tot, a_off_tot, N_L);
}

void Reinitiliazise_Affinity(MATRIX_DOUBLE *a_on, MATRIX_DOUBLE *a_off, MATRIX_DOUBLE *a_d, std::vector<double> *a_d_tot, std::vector<double> *a_on_tot, std::vector<double> *a_off_tot, std::vector<double> *affinity, MATRIX_INT *L_R, std::vector<double> *X_CM_x, std::vector<double> *X_CM_y, std::vector<double> *X_CM_theta, std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<double> *X_r_x, std::vector<double> *X_r_y, std::vector<int> *L_to_R, std::vector<int> T_l, std::vector<int> T_r, double N_L)
{/* Update the affinities, by emptying/clearing the lists and
	* reinitilizing all again.
	*/
	// CLEAR everything
	for(int i=0; i<N_L; i++)
    {
        (a_d->at(i)).clear();
        (a_on->at(i)).clear();
        (a_off->at(i)).clear();
    }

    a_d_tot->clear();
    a_on_tot->clear();
    a_off_tot->clear();
    affinity->clear();

	// INIT again
	Initialize_Affinity(a_d, L_R, X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x, X_r_y, L_to_R, T_l, T_r, k_d, N_L, 0, 1);
	Initialize_Affinity(a_on, L_R, X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x, X_r_y, L_to_R, T_l, T_r, k_on, N_L, 1, 1);
    Initialize_Affinity(a_off, L_R, X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x, X_r_y, L_to_R, T_l, T_r, k_off, N_L, 2, 2);

    *a_d_tot = Get_a_X_tot(a_d, N_L);
	*a_on_tot = Get_a_X_tot(a_on, N_L);
	*a_off_tot = Get_a_X_tot(a_off, N_L);

	*affinity = Get_a_tot(a_d_tot, a_on_tot, a_off_tot, N_L);
}



//___________________________________________________________________________________________________________________________
// CONTINUATION
using namespace std; 

std::vector<double> load_CM(string dir){
	/* Returns a std:vector of four elements, listing the CM_x, CM_y, CM_theta and N_b 
	 * of the last step of the previous sim.
	 * Needs to be given the correct directory to find the file 'trajectory.txt'.
	 * 
	 */
	
	std::vector<double> CM;
	CM.clear();
	
	 // Open the input file named "trajectory.txt" 
    ifstream inputFile(dir+ "trajectory.txt"); 
    string line;
    // skip to the last line
	while (inputFile >> std::ws && getline(inputFile, line)){;} 
    inputFile.close();     // Close the file 
	
	int pos = 0;	
	for(int i= 0; i<4; i++){
		pos= line.find(',');
		CM.push_back(stof(line.substr(0, pos)));
		line.erase(0, pos+1);
	}
	return CM;
}

void load_ligands(string dir, std::vector<double> *X_l_x, std::vector<double> *X_l_y, std::vector<int> *T_l, string file= "ligands.xyz" ){
	/* Loads the coordinates of the "sphericalVirus.xyz" file 
	 * (or any other given file) from the given directory as well
	 * as their types, and packs them in the vectors X_r_x, X_r_y and T_l.
	 * 
	 * In T_l: "1" denote HA ligands
	 * 			"0" denote NA ligands
	 */
	
	// Open the input file named "receptor.txt" 
    ifstream inputFile(dir+ file); 
	
	// extract the Last and second to last line of the file
    string line;
    getline(inputFile, line);	// skip the two header lines	
	getline(inputFile, line);
    
	// push the values to the vectors
	int pos = 0;		
	while(getline(inputFile, line)){
		T_l->push_back(stoi(line.substr(0, 1)));
		line.erase(0, 2);

		pos= line.find('\t');
		X_l_x->push_back(stof(line.substr(0, pos)));
		line.erase(0, pos+1);
		X_l_y->push_back(stof(line));
	}
    inputFile.close();     // Close the file 
 
	std::cout << "Loaded all ligands. Number of ligands: X_l_x.size()= " << X_l_x->size() << '\n';
}

void load_receptors(string dir, std::vector<double> *X_r_x, std::vector<double> *X_r_y){
	/* Loads the coordinates of the "receptor.txt" file from the given 
	 * directory and packs them in the vectors X_r_x and X_r_y.
	 */
	
	// Open the input file named "receptor.txt" 
    ifstream inputFile(dir+ "receptor.txt"); 
    ifstream copy(dir+ "receptor.txt"); 
	
	// extract the Last and second to last line of the file
    string x_line, y_line;
    getline(inputFile, y_line);		// 'x_line' lags one line behind 'y_line'
    while (inputFile >> std::ws && getline(inputFile, y_line)){
		getline(copy, x_line);
	} 
    inputFile.close();     // Close the file 
    
	// skip 'X_r_x = [' -> index 9 (might have to get changed, depending on save format))
	x_line= x_line.substr(9, x_line.size()-9-1);
	y_line= y_line.substr(9, y_line.size()-9-1);
	
	// push the values to the vectors
	int x_pos, y_pos = 0;		
	while(x_line.find(',') != string::npos){
		x_pos= x_line.find(',');
		X_r_x->push_back(stof(x_line.substr(0, x_pos)));
		x_line.erase(0, x_pos+1);

		y_pos= y_line.find(',');
		X_r_y->push_back(stof(y_line.substr(0, y_pos)));
		y_line.erase(0, y_pos+1);	
	}
	
	if (y_line.find(',') != string::npos){
		std::cout << "We have more Receptor Y-coords, than receptor X-coords. Something went wrong saving the data." << '\n';
		std::cout << "Components of Y vector not loaded: " << y_line << '\n';
	} else {std::cout << "Loaded all receptors. Number of receptors: X_r_x.size()= " << X_r_x->size() << std::endl;}
}

void load_SA(string dir, std::vector<int> *T_r){
	/* Use the data in depleted_SA.txt to switch the SA in the list T_r
	 * inactive ("0").
	 *
	 * Keep in mind the index shift in T_r:
	 * i.e.: depleted_SA.txt : [500, 503, 520, ....] -> T_r[499]= 0, T_r[502]= 0, T_r[519]= 0
	 */
	
	 // Open the input file named "trajectory.txt" 
    ifstream f_SA(dir+ "depleted_SA.txt"); 
    string line;
    int pos;

    getline(f_SA, line);						// skip first line (Header)
	while (getline(f_SA, line)){				//f_SA >> std::ws && 
		pos= line.find(',');							// substr(pos+1)= index of receptor
		T_r->at(stoi(line.substr(pos+1))  - 1)=0; 		// adjust type of receptor to "off"/ cleaved
	} 
    f_SA.close();     // Close the file 
    
    std::cout << "Loaded depleted SA from: " << dir << std::endl;
}

void load_bridges(string dir, std::vector<int> *L_to_R, std::vector<int> *T_l, std::vector<int> *T_r){
	/* Use the data in "L_to_R.txt" in order to fill L_to_R which
	 * denotes to each ligand the index of the bound receptor 
	 * or "-1" if none is bound.
	 * Updates the vectors T_l and T_r by denoting the corresponding
	 * bound particles (with a "2").
	 * 
	 * Keep in mind the index shift in T_r:
	 * 	i.e.: L_to_R= [..., 100, ..., 125, ...] -> T_r[99]=2, T_r[124]=2
	 */
	 
	 // Open the input file named "L_to_R.txt" 
    ifstream inputFile(dir+ "L_to_R.txt"); 
    string line;
    // skip to the last line
	while (inputFile >> std::ws && getline(inputFile, line)){;} 
    inputFile.close();     // Close the file 
	
	int pos = 0;	
	while(line.find(',')!= string::npos){
		pos= line.find(',');
		L_to_R->push_back(stof(line.substr(0, pos)));
		line.erase(0, pos+1);
	}
	std::cout << "Loaded L_to_R." << std::endl;
	
	int n_br=0;
	for (int l= 0, end= L_to_R->size(); l < end; l++){
		if (L_to_R->at(l)!= -1)
			{n_br +=1;
			T_l->at(l)= 2;
			T_r->at(L_to_R->at(l) -1)=2;}		//L_to_R->at(l)= index_rec
	}
	
	std::cout << "Loaded bridges also in T_l and T_r. Number of bridges: " << n_br << std::endl;
}

void load_CBs(string dir, std::vector<int> *list_CB){
	/* Use the data in "CB.txt" in order to fill list_CB which
	 * denotes the indexes of the receptors of the constraining
	 * bridges.
	 */
	 
	 // Open the input file named "CB.txt" 
    ifstream inputFile(dir+ "CB.txt"); 
    string line;
    // skip to the last line
	while (inputFile >> std::ws && getline(inputFile, line)){;} 
    inputFile.close();     // Close the file 
	
	int pos = 0;	
	while(line.find(',')!= string::npos){
		pos= line.find(',');
		list_CB->push_back(stof(line.substr(0, pos)));
		line.erase(0, pos+1);
	}
	if (list_CB->size()==0){std::cout << "\nWARNING!! List_CB is empty." << std::endl;}
	else{std::cout << "Loaded list_CB. Number of constraing bridges: " <<  list_CB->size() << std::endl;}
}

// --------------------------- --------------------------- --------------------------- --------------------------- --------------------------- --------------------------- --------------------------- --------------------------- 
//																						MAIN
// --------------------------- --------------------------- --------------------------- --------------------------- --------------------------- --------------------------- --------------------------- --------------------------- 
int main(int argc, char* argv[])
{ /*
	* Use in the following format:
	* ./Run <K_d>  [<Run Iter>] 
	* ("Iter" gets one out of 10 random seeds assigned.)
	*/
	std::cout << "Run " << argv[0] << std::endl;
	
	long idum;			// the random seed    (has to be < 0)
	// ------------------------------------ PARSE INPUT --------------------------- --------------------------- --------------------------- --------------------------- 
	if (argc == 3){
		std::cout << "1st Input: " << argv[1] << std::endl;
		k_d= atof(argv[1]);
		int iter= std::stoi(argv[2]);
		if (iter==1){idum= -111142451;}
		else if (iter==2){idum= -211142450;}
		else if (iter==3){idum= -311190231;}
		else if (iter==4){idum= -487658742;}
		else if (iter==5){idum= -523741243;}
		else if (iter==6){idum= -672237412;}
		else if (iter==7){idum= -777777772;}
		else if (iter==8){idum= -818237586;}
		else if (iter==9){idum= -923245123;}
		else if (iter==10){idum= -102355612;}
		else {std::cout << "Iteration <" << iter << "> has not been specified yet. Exit programm."<< std::endl;return 0;}
	}else if (argc == 2){
		k_d= atof(argv[1]);
		srand (time(NULL));
		idum =  -1*(rand() % 888888888 + 111111111); //-111142450;-111121278;-111121278;-111140767; -111121051; 
	}
	else{
		std::cout << "Wrong amount of arguments.\n\
		use Format '/Run <K_d> [<Run Iter>]'.\nExit program."<< std::endl;
		return 0;}
	// ------------------------------------ PARSE INPUT --------------------------- --------------------------- --------------------------- --------------------------- 
	std::cout << "seed du GNPA: " << idum << '\n';
	std::cout << " ---------------------------------------------- " << std::endl;
    
	// If the simulation is no continuation: init type of particle and polarization
	// if it is a continuation: provide path to ligand file and ALL depleted_SA .txt files
	bool diffusing_receptors= true;
	bool continuation=false;
	string dir_cont= "./";
	int time_cont=0;

	std::cout << "DT_Brown: " << DT_Brown <<  '\n';
	std::cout << "DT: " << DT <<  '\n';
	std::cout << "N_Brown: " << N_Brown << '\n';

	std::cout << "D_receptor: " << D_receptor  << "\t nm^2 /s" << '\n';
	std::cout << "D_sphere*: " << D_sphere  << "\t nm^2 /s" << '\n';
	std::cout << "D_r*: " << D_r  << "\t nm^2 /s" << '\n';
	
	
	std::cout << "Simulation of a 2D system with diffusing receptors" << "\n\n\n" << std::flush;
	clock_t startTime = clock();
	
	std::string annex= "";
	if (continuation){std::string annex= ""; 
		std::cout << "Start time of continuation [ms]: " <<  time_cont << '\n';}  //std::string annex= "_continued"
	
	// init output files
    std::ofstream trajectory;
    std::ofstream receptor;
    std::ofstream run_info;
    std::ofstream depleted;
    std::ofstream CB;
    std::ofstream ltr;
    std::ofstream bl;

    trajectory.open ("trajectory" + annex + ".txt");
    receptor.open ("receptor" + annex + ".txt");
    run_info.open ("run_info" + annex + ".txt");
    depleted.open("depleted_SA" + annex + ".txt");
    ltr.open("L_to_R" + annex + ".txt");
    bl.open("bridges" + annex + ".txt");

    trajectory << "TIME,X_CM,Y_CM,THETA_CM,N_b" << '\n';
    depleted << "STEP,SA" << '\n';

    run_info << "idum: " << idum << '\n';
    run_info << "DT: " << DT << "  N_dT: " << N_dT <<  '\n';
    run_info << "System Size: " << Taille_Systeme << "  R_Verlet: " << R_Verlet << '\n'; //"  Background Receptors: " << N_R << 
    run_info << "k_d: " << k_d << "  k_off: " << k_off << "  k_on: " << k_on<< '\n';
    
    // ----- GET REC DENSITY BELOW PARTICLE:
    std::ofstream RecDensity;   
    RecDensity.open ("RecDensity" + annex + ".txt");
	RecDensity << "step,RecDensity,ActiveDensity" << std::endl;
    // -----
    
	// init variables
	std::vector<double> a_d_tot;
	std::vector<double> a_on_tot;
	std::vector<double> a_off_tot;

	std::vector<double> affinity;

	// Lists to access the Ligands/Receptors in each simulation Cell
	MATRIX_INT v_L;
	std::vector<int> ll_L;
	std::vector<int> bl_L;
	std::vector<int> hoc_L;
	MATRIX_INT v_R;
	std::vector<int> ll_R;
	std::vector<int> bl_R;
	std::vector<int> hoc_R;

	double R_Cell;
	int N_x;
	int N_y;
	int N_Cell;

	MATRIX_INT Neigh;

	// X_l: coordinates of the Ligands on the Virium (virial frame of reference)
	std::vector<double> X_l_x;  
	std::vector<double> X_l_y;
	
	// X_CM: coordinates of the Ligands in the lab frame of reference 
	std::vector<double> X_CM_x; 
	std::vector<double> X_CM_y;
	std::vector<double> X_CM_theta; // orientation of the virium (angle between x axes?)
	
	// X_r: coordinates of the Receptors (SA) on the cell (Labframe of ref.)
	std::vector<double> X_r_x;
	std::vector<double> X_r_y; 
	
	std::vector<int> T_l;	// Ligand type: 0=NA, 1=SA, 2=bound/bridge
	std::vector<int> T_r;	/* Receptor type: 0= cleaved, 1= active, 2=bound/bridge, 
										* !!Attention: the indexes are shifted in T_r!!
										* The indeces denoted in L_to_R, list_bridged_receptors, etc... are shifted by +1 as compared to the 
										* T_r vector. 
										* i.e.: L_to_R= [..., 100, ..., 125, ...] -> T_r[99]=2, T_r[124]=2
										* 		depleted_SA.txt : [500, 503, 520, ....] -> T_r[499]= 0, T_r[502]= 0, T_r[519]= 0
										*/
	
	std::vector<int> L_to_R; 					// lists for each ligand the index of the bonded receptor, "-1" when unbound
	std::vector<int> list_bridged_receptors;	// lists the indexes of all receptors which currently forming a bridge
	std::vector<int> list_CB;					// the indexes of all constraining bridges (indexes of the involved receptors)
	
	double time_react;
	double time_tot;
	std::vector<double> time_vect;
	time_vect.push_back(0);
	
	int N_b;
    std::vector<int> bridge_vect;
    bridge_vect.push_back(0.0);

	// ----------------------------------------------------------------------------------------------
	//  --------------------------------- INITIALIZE: -----------------------------------------------
	// ----------------------------------------------------------------------------------------------
	if (!continuation){
		// ________
		// IAV: initialize Ligand grid on IAV
		// initialize_IAV(&X_l_x, &X_l_y, &T_l, &L_to_R, &idum);			// Bacilla Virus
		// initialize_sphericalIAV(&X_l_x, &X_l_y, &T_l, &L_to_R, &idum, true);	// Spherical clustered
		// initialize_sphericalIAV(&X_l_x, &X_l_y, &T_l, &L_to_R, &idum, false);	// Spherical unclustered

		// Load given Ligand grid
		string f= "./ligands_equally.xyz";
		load_ligands(dir_cont, &X_l_x, &X_l_y, &T_l, f);
		
		for (int i=0, end= T_l.size(); i < end; ++i){
			L_to_R.push_back(-1);	// -1= no receptor is bound to ligand
		}
		// ________
		// SA: Initialize Receptors (random distr. on a squared grid + clusters[if indicated])
		initialize_SA(&X_r_x, &X_r_y, &T_r, &idum, false);	// not clustered
		// initialize_SA(&X_r_x, &X_r_y, &T_r, &idum, true);	// clustered
		
		//_________
		// Init CM:
		X_CM_x.push_back(0.0);
		X_CM_y.push_back(0.0);
		X_CM_theta.push_back(0.0);
	}
	else{// We Continue a Stopped simulation
		// IAV: initialize Ligand grid on IAV
		string f= "ligands.xyz";
		load_ligands(dir_cont, &X_l_x, &X_l_y, &T_l, f);
		
		// SA: Initialize Receptors
		load_receptors(dir_cont, &X_r_x, &X_r_y);
		
		// init T_r vector with 1s ("active"/not cleaved)
		for (int j = 0, end= N_R; j < end; j++)
			{T_r.push_back(1);}
		// load the inactive SA
		load_SA(dir_cont, &T_r);
		
		// Init CM:
		std::vector<double> old_CM;
		old_CM= load_CM(dir_cont);
		
		// old_CM[0] = end time
		X_CM_x.push_back(old_CM[1]);
		X_CM_y.push_back(old_CM[2]);
		X_CM_theta.push_back(old_CM[3]);
}
	// ----------------------------------------------------------------------------------------------

   	int N_L = X_l_x.size();		// Number of Ligands
	// print initial receptors
	receptor << "X_r_x = [";
	for(int i = 0, end= X_r_x.size(); i<end; i++)
		{receptor << X_r_x[i] << ',';}
	receptor << "]" << '\n' << "X_r_y = [";
	for(int i = 0, end=X_r_y.size(); i<end; i++)
		{receptor << X_r_y[i] << ',';}	
	receptor << "]";
	receptor.flush();

	// Test number of Receptors
	if (N_R != X_r_x.size()){
		std::cout << "WARNING: While Initializing Receptors. Number of Receptors does not correspond to the desired value N_R.\n\n";
		return 0;}
	else if (X_r_x.size() != X_r_y.size()){
		std::cout << "WARNING: While Initializing Receptors. Number of Receptors X coordinates does not correspond to the Y coordinates.\n\n";
		return 0;}
    run_info << "N_Ligands: " << N_L << "  N_Receptors: " << N_R << "  N_NA: See .out-file" << '\n';
    run_info.flush();

	
	// L_R: listing the interacting(those in range) receptors for each ligands, R_L: listing the interacting ligands for each receptor 
	MATRIX_INT L_R; 
	L_R.resize(N_L, std::vector<int>(0));
	MATRIX_INT R_L;
	R_L.resize(N_R, std::vector<int>(0));		// N_R: Number of Receptors, defined as a constant
	MATRIX_INT L_R_old; 
	L_R_old.resize(N_L, std::vector<int>(0));

	// Affinities
	MATRIX_DOUBLE a_d;							
	a_d.resize(N_L, std::vector<double>(0));
	MATRIX_DOUBLE a_on;
	a_on.resize(N_L, std::vector<double>(0));
	MATRIX_DOUBLE a_off;
	a_off.resize(N_L, std::vector<double>(0));


    
	// SETUP SYSTEM-GRID
	// Divide System in Cells of size R_verlet
	N_x = int(Taille_Systeme/R_Verlet); //D->R_Verlet
	N_y = N_x;
	N_Cell = N_x*N_y;

	R_Cell = Taille_Systeme/N_x;	//= R_verlet

	// create a neighbour list
	Neigh = Set_Neigh(N_Cell, N_x, N_y);

	// INIT LL,BL,HOC:
	// lists the elements present in the simulation cells. hoc assigns each sim.cell the last
	// element, ll[and bl] allows to access the rest
	// The Ligands:
	// The Ligands (we need ligands in the labframe of reference):
	std::vector<double> x_h;
	std::vector<double> y_h;
	double th;
	th = X_CM_theta[X_CM_theta.size() - 1];
	for (int alpha = 0; alpha < N_L; ++alpha)
	{ 
		// Calculate ligand coordinates in Labframe of ref.:
	    x_h.push_back((X_l_x[alpha])*cos(th) - (X_l_y[alpha])*sin(th) + X_CM_x[X_CM_x.size() - 1]);
		y_h.push_back((X_l_x[alpha])*sin(th) + (X_l_y[alpha])*cos(th) + X_CM_y[X_CM_y.size() - 1]);
	}
	v_L = Init_CellLists(R_Cell, N_Cell, N_L, N_x, N_y, &x_h, &y_h);

	x_h.clear();
	y_h.clear();
	
	ll_L = v_L[0];
	bl_L = v_L[1];
	hoc_L = v_L[2];
	// The Receptors
	v_R = Init_CellLists(R_Cell, N_Cell, N_R, N_x, N_y, &X_r_x, &X_r_y);
	ll_R = v_R[0];
	bl_R = v_R[1];
	hoc_R = v_R[2];

	// Initialize Verlet Lists: L_R and R_L
	Initialize_L_R(&L_R, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_R, ll_R, T_l, L_to_R, N_x, N_y, R_Cell, N_L, 0);
	Initialize_R_L(&R_L, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_L, ll_L, T_r, T_l, L_to_R, N_x, N_y, R_Cell, N_L, 0); 
	
	//_____________________________________________________________________
	//_____________________________________________________________________	
	if (continuation){
		// get L_to_R:
		load_bridges(dir_cont, &L_to_R, &T_l, &T_r);
		// get bridges_list:
		list_bridged_receptors.clear();
    	list_bridged_receptors = Get_bridges(&T_l, &L_to_R);
    	std::cout <<"\n\n\n" << std::endl;
	}
	//_____________________________________________________________________
	//___________________________________________________________________	

	// Initialize affinities (for each ligand there is a list of affinities for all receptors in verlet range)
	Initialize_Affinity(&a_d, &L_R, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &L_to_R, T_l, T_r, k_d, N_L, 0, 1);
	Initialize_Affinity(&a_on, &L_R, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &L_to_R, T_l, T_r, k_on, N_L, 1, 1);
    Initialize_Affinity(&a_off, &L_R, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &L_to_R, T_l, T_r, k_off, N_L, 2, 2);

	// get affinity vectors (a_tot for each vector)
    a_d_tot = Get_a_X_tot(&a_d, N_L);
	a_on_tot = Get_a_X_tot(&a_on, N_L);
	a_off_tot = Get_a_X_tot(&a_off, N_L);
	// a_tot to choose the type of reaction	
	affinity = Get_a_tot(&a_d_tot, &a_on_tot, &a_off_tot, N_L);

	// Init Diffusion variables (Necessary to calculate maximum displacement for 
	// ligands and receptors to check whether its necessary to update the Cell lists)
  	std::vector<double> displacement_L;
  	std::vector<double> x_L;				// ligands after time step
  	std::vector<double> y_L;
  	std::vector<double> x_init_L;		// ligands before time step
  	std::vector<double> y_init_L;
  	double max_L = 0.0;
  	for(int w=0; w<N_L;w++)
  	{
  		x_init_L.push_back(0.0);
  		y_init_L.push_back(0.0);
  		x_L.push_back(0.0);
  		y_L.push_back(0.0);
  		displacement_L.push_back(0.0);
  	}

	double X_CM_x_init = X_CM_x[X_CM_x.size() - 1];
	double X_CM_y_init = X_CM_y[X_CM_y.size() - 1];
	double X_CM_theta_init = X_CM_theta[X_CM_theta.size() - 1];
	
  	double theta;

  	std::vector<double> displacement_R;
  	std::vector<double> displacement_close_R;
  	std::vector<double> x_init_R;			// receptors before time step
  	std::vector<double> y_init_R;
  	double max_R = 0.0;
  	for(int w=0; w<N_R;w++)
  	{
  		x_init_R.push_back(0.0);
  		y_init_R.push_back(0.0);
  		displacement_R.push_back(0.0);
  	}

	// initialize labframe coordinates of the ligands and receptors
  	for(int w=0; w < N_L; w++){
		theta = X_CM_theta_init;
		x_init_L[w] = X_l_x[w]*cos(theta) - X_l_y[w]*sin(theta) + X_CM_x_init;
	    y_init_L[w] = X_l_x[w]*sin(theta) + X_l_y[w]*cos(theta) + X_CM_y_init;
	}
	for(int w=0;w < N_R; w++){
		x_init_R[w] = X_r_x[w];
		y_init_R[w] = X_r_y[w];
	}

	std::vector<double> new_rCM;
	std::vector<int> bridges_therm;
	int NBrow;

	int gen = 0;
	
	// Detachment Conditions
	int tau_threshold= 50;		// micro s| if the particle takes longer to attach, we stop the sim 
	bool Attached= 0; 			// flag that shows whether the particle is already attached to the cell surface
	bool Detached= 0;			// flag that shows whether the particle detached from the cell surface (after attaching), condition for exiting the 
	int time_detach, time_attach, current_time;
	int total_reaction=0;

	// small safeguard for small DT
	if (DT < 1. * 1e-6)
		{std::cout << "WARNING: DT is smaller than 1 µs. \n A lot of the output and timechecks are in µs, using smaller DT might lead to problems." 
		"\n At the very least make 'current_time-var' a float and try to proceed after. Exit Programm..." << std::endl; exit(0);}

	run_info << "\nInitialized everything (IAV, Receptors, Neighbourlist, etc.) in  " << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds.";
	run_info.flush();
	// ----------------------------------------------------------------------------------------------
	// ------------------------ START SIM: ----------------------------------------------------------
	// ----------------------------------------------------------------------------------------------
    for (int step_i = 1; step_i < N_dT + 1; ++step_i){
		time_tot = 0.0;
		int count_react= 0;
		current_time= time_cont + (step_i * (DT * 1e6)); 		// current time in µs
		//___________________________________________________________________________________________________________________________
		// REACTIONS: list/perform reactions of one time step, get constraining bridges
		while(time_tot < DT){
    		int r = Choose(affinity, &idum); // choose reaction, r: 1= destroy SA, 2= forming a bridge, 3= breaking a bridge
    	    time_react = Get_time_reaction(affinity[0], &idum);
			time_tot += time_react;
			count_react += 1;
    	    if(time_tot < DT){
				// do reacton and get new affinities    
        	    Make_reaction(&a_d, &a_on, &a_off, a_d_tot, a_on_tot, a_off_tot, L_R, R_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &T_l, &T_r, &L_to_R, &list_bridged_receptors, &list_CB, r, &idum, step_i, &depleted);
 				
                a_d_tot.clear();
                a_on_tot.clear();
                a_off_tot.clear();
                affinity.clear();
        	    
        	    a_d_tot = Get_a_X_tot(&a_d, N_L);
                a_on_tot = Get_a_X_tot(&a_on, N_L);
            	a_off_tot = Get_a_X_tot(&a_off, N_L);
            	affinity = Get_a_tot(&a_d_tot, &a_on_tot, &a_off_tot, N_L);

				// ------------------------
				// ATTACHMENT | DETACHMENT
				N_b = list_bridged_receptors.size();	

				if (Attached==0 and N_b >= 1) // ATTACHMENT
					{Attached= 1; time_attach= current_time; std::cout << "\nParticle attached at step " << step_i << ", t= " << current_time
					<< " µs, ('Substep': time_tot= " << time_tot << ", time_react= " << time_react << ")" << std::endl;}
				else if (Attached==1 and N_b == 0)// DETACHMENT (only if particle is already attached)
					{Detached= 1; time_detach= current_time; std::cout << "\nParticle detached at t= " << time_detach 
					<< " µs, ('Substep': time_tot= " << time_tot << ", time_react= " << time_react << ")" << std::endl; break;} // stop sim
			}
		}
		total_reaction += count_react;
    	depleted.flush();
		// ----------------------
		// Treating DETACHMENT:
		if (Detached==1){std::cout << "\n Stop simulation." << std::endl; break;} // stop sim
		else if (Attached==0 and (current_time >= tau_threshold))// Particle takes too long to attach | stop sim
			{std::cout << "\n Particle takes too long to attach. tau_theshold= "<< tau_threshold << ", Current Time=" << current_time
			 << "\n Exit Code " << std::endl; Detached= 1; time_detach= current_time; break;}


		// count number of bridges and save them in "bridge_vect"
		N_b = list_bridged_receptors.size();
		bridge_vect.push_back(N_b);	
			
		//----------------------------------------------------------------------------------------------
		//		DIFFUSION:
		double D_II, D_ortho, D_rot;
		// UPDATE the DIFFUSION CONSTANTS:
		if (Attached){
			D_II= D_sphere/N_b;				// D_parallel/N_b;
			D_ortho= D_sphere/N_b;			// D_orthogonal/N_b;
			D_rot= D_r/N_b;					// D_THETA/N_b;	
		} else //avoid division by 0 bridges 
			{D_II= D_sphere, D_ortho= D_sphere, D_rot= D_r;}
		
		// PRINT Diffusion parameters
	    if((step_i) % (N_dT /25) == 0){ // print 25 times
			std::cout << "\nt= " << current_time << " µs, step: " <<  step_i  << '\n';
			std::cout << "N_b: \t" <<  N_b << '\n';
			std::cout << "D_II: \t" <<  D_II << "\t nm^2 /s" << '\n';
			std::cout << "D_ortho: \t" <<  D_ortho  << "\t nm^2 /s" << '\n';
			std::cout << "D_rot: \t" <<  D_rot << "\t nm^2 /s" << '\n';
			std::cout << "Reactions in current step: \t" <<  count_react << '\n';
			
		}

		// diffuses the Virus Particle and the bound receptors
		diffusion(&X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &T_l, &L_to_R, &idum, DT, D_II, D_ortho, D_rot); // includes rotation

		//----------------------------------------------------------------------------------------------	
		// RECEPTOR DIFFUSION:
		if (diffusing_receptors){
			// diffusion substeps
			for(int iBrow=0; iBrow < N_Brown; iBrow++){
				// diffuses the bound receptors by Brownian motion (count rejections in "rejections")
				diffusion_bound_receptors(&X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &T_l, &L_to_R, &list_bridged_receptors, &idum, DT_Brown);
			}
			// diffuse the unbound receptors
			diffusion_loose_receptors(&X_r_x, &X_r_y, &list_bridged_receptors, &idum, DT);
		}
		if(X_CM_x[X_CM_x.size()-1] != X_CM_x[X_CM_x.size()-1])
			{std::cout << "\n\n Something is wrong. We should not enter here! I" 
			<< " \n Step= " << step_i << "t= " << current_time << " µs , X_CM_x[X_CM_x.size()-1]= " << X_CM_x[X_CM_x.size()-1] << std::endl; return 0;}
		if(X_CM_y[X_CM_y.size()-1] != X_CM_y[X_CM_y.size()-1])
			{std::cout << "\n\n Something is wrong. We should not enter here! II" 
			<< " \n Step= " << step_i << "t= " << current_time << " µs , X_CM_y[X_CM_y.size()-1]= " << X_CM_y[X_CM_y.size()-1] << std::endl; return 0;}
		if(X_CM_theta[X_CM_theta.size()-1] != X_CM_theta[X_CM_theta.size()-1])
			{std::cout << "\n\n Something is wrong. We should not enter here! III" 
			<< " \n Step= " << step_i << "t= " << current_time << " µs , X_CM_theta[X_CM_theta.size()-1]= " << X_CM_theta[X_CM_theta.size()-1] << std::endl; return 0;}


		//___________________________________________________________________________________________________________________________		
		//UPDATES:			
		new_rCM.clear();

		// check breakage after receptor diffusion
		bool check_bk = false;
		check_bk = check_breakage(&X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &T_l, &L_to_R, 0.0, 0.0, 0.0);
		if(check_bk == true)
		{
			std::cout << "\n" << "overextended bridge 2nd test in step: " << step_i << ", t= " << current_time << " µs \n";
			std::cout << "x_CM = [" << X_CM_x[X_CM_x.size()-1] << "," << X_CM_y[X_CM_y.size()-1] << "]" << "\n";
			return 0;
		}

		// get the displacements of the Ligands...
		for(int w=0; w < N_L; w++)
		{
			theta = X_CM_theta[X_CM_theta.size() - 1];
			x_L[w] = X_l_x[w]*cos(theta) - X_l_y[w]*sin(theta) + X_CM_x[X_CM_x.size() -1];
		    y_L[w] = X_l_x[w]*sin(theta) + X_l_y[w]*cos(theta) + X_CM_y[X_CM_y.size() -1];
			displacement_L[w] = pow(x_L[w]-x_init_L[w],2) + pow(y_L[w]-y_init_L[w],2);
		}
		// ... and the receptors
		for (int w=0; w < N_R; w++){// We only look at receptors "close" to the Virion
			if(pow(X_r_x[w]-X_CM_x[X_CM_x.size()-1],2)+pow(X_r_y[w]-X_CM_y[X_CM_y.size()-1],2) < R_close*R_close )
			{
				displacement_R[w] = pow(X_r_x[w]-x_init_R[w],2) + pow(X_r_y[w]-y_init_R[w],2);
				displacement_close_R.push_back(displacement_R[w]);
			}
		}
		// Get the largest displacements for R and L
		max_L = sqrt(*max_element(displacement_L.begin(), displacement_L.end()));
		max_R = sqrt(*max_element(displacement_close_R.begin(), displacement_close_R.end()));

		// if it is possible that something entered/exited a verlet list
		if(max_R + max_L > R_Verlet - D) 
	    {// update the Cells and Verlet-lists

	    	Update_Cell_L(&ll_L, &bl_L, &hoc_L, &X_l_x, &X_l_y, &X_CM_x, &X_CM_y, &X_CM_theta, N_L, N_x, N_y, R_Cell, X_CM_x_init, X_CM_y_init, X_CM_theta_init);
	    	Update_Cell_R(&ll_R, &bl_R, &hoc_R, &X_r_x, &X_r_y, &x_init_R, &y_init_R, N_x, N_y, R_Cell);
			Reinitialize_Lists(&L_R, &R_L, &a_d, &a_on, &a_off, &a_d_tot, &a_on_tot, &a_off_tot, &affinity, Neigh, hoc_R, hoc_L, ll_L, ll_R, X_l_x, X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, N_x, N_y, R_Cell, N_L, T_l, T_r, &L_to_R, step_i);

			// get new init_L, init_R and X_CM_init
			for(int w=0; w < N_L; w++)
			{	
				x_init_L[w] = x_L[w];
				y_init_L[w] = y_L[w];					
	    	}

	    	X_CM_x_init = X_CM_x[X_CM_x.size() - 1];
	    	X_CM_y_init = X_CM_y[X_CM_y.size() - 1];
	    	X_CM_theta_init = X_CM_theta[X_CM_theta.size() - 1];

	    	for(int w=0; w< N_R; w++)
	    	{
	    		x_init_R[w] = X_r_x[w];
	    		y_init_R[w] = X_r_y[w];
	    	}
	    	displacement_close_R.clear();	    
	    }

		// update the affinities
    	else
    	{Reinitiliazise_Affinity(&a_on, &a_off, &a_d, &a_d_tot, &a_on_tot, &a_off_tot, &affinity, &L_R, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &L_to_R, T_l, T_r, N_L);}

    	// -------- ---- PRINTS --- -------
		// if((step_i) % (N_dT /200) == 0)  // get 200 trajectory frames
		if((current_time) % (10) == 0)		// every 10 µs
		{ // ovito compliant (xyz)
	    	trajectory << current_time << ','  << X_CM_x.at(X_CM_x.size() - 1) << "," << X_CM_y.at(X_CM_y.size() - 1) << "," << X_CM_theta.at(X_CM_theta.size() - 1) << "," << bridge_vect.at(bridge_vect.size() - 1) << '\n';
	      	trajectory.flush();		
	      	//_____
			int close_Rec= 0;
			int closeActive_Rec= 0; 
			for (int w=0; w < N_R; w++){// Count the receptors "close" to the Virion
				if(pow(X_r_x[w]-X_CM_x[X_CM_x.size()-1],2)+pow(X_r_y[w]-X_CM_y[X_CM_y.size()-1],2) < R_close*R_close )
					{++close_Rec; if (T_l[w]==1){++closeActive_Rec;}}
			}
	      	RecDensity << current_time << ','  << close_Rec/(M_PI*R_close*R_close) << ','  << closeActive_Rec/(M_PI*R_close*R_close) << '\n';
	      	RecDensity.flush();	
			//_____
						
		    for(int k = 0, end= list_bridged_receptors.size(); k<end; k++)
				{bl << list_bridged_receptors[k] << ",";}
		    bl << '\n';

		    for(int k = 0, end=L_to_R.size(); k<end; k++)
				{ltr << L_to_R[k] << ",";}
		    ltr << '\n';
			
			if (diffusing_receptors and (current_time % 1000 == 0)){	// rec_pos every ms
				receptor << '\n' << "X_r_x = [";
				for(int i = 0, end=X_r_x.size(); i<end; i++)
					{receptor << X_r_x[i] << ',';}
				receptor << "]" << '\n' << "X_r_y = [";
				for(int i = 0, end=X_r_y.size(); i<end; i++)
					{receptor << X_r_y[i] << ',';}	
				receptor << "]";
				receptor.flush();
			}
		    bl.flush();
		    ltr.flush();
	    }
	}
	// -----------------------------------------------------------------------------------------------
	// SIMULATION FINISHED/ STOPPED
	int sim_end_time;
	if (Detached==0){
		sim_end_time= time_cont + (N_dT * (DT * 1e6)); 		// current time in µs;
		std::cout << "\n--------------------------------------------------------------------------\n";
		std::cout << "Finished Run without detachment." << std::endl;
		// write time data
		run_info << "\n\nRun finished, No detachment: ";
	}
	else{ // Detached ==1
		sim_end_time= time_detach;
		std::cout << "\n--------------------------------------------------------------------------\n";
		std::cout << "Finished Run with detachment at step= " << time_detach << std::endl;
		// write time data
		run_info << "\n\nRun Stopped after Detachment at step= " << time_detach << ": ";
	}

	// write Rejections data
	run_info << "\nDuration of run, real time [s]: " << double( clock() - startTime ) / (double)CLOCKS_PER_SEC << ",\t Simulation time: " << sim_end_time << '\n';
	run_info << "\nSimulation speed: " << (double)sim_end_time * 3.6 * (double)CLOCKS_PER_SEC/ (double)(clock() - startTime ) << " seconds/ hourCPU";
	run_info << "\nRejections: " << rejections << "\tTotalTrials: "<< total_tests << "\t-> rejection rate: " << (double)rejections/(double)total_tests << '\n';
	run_info << "Total number of Reactions: " << total_reaction << ",\t Reactions per step: " << (double)total_reaction/(double)sim_end_time << '\n';

	run_info.flush();
	depleted.close();
	trajectory.close();
	receptor.close();
	run_info.close();


	// get Mean number of bridges
	// int num_sim_steps= sim_end_time - time_cont /(DT*1e6);
	int num_sim_steps=bridge_vect.size();
	int sum_bridges=0;
	for (int t_i=0; t_i < num_sim_steps; ++t_i){sum_bridges+= bridge_vect[t_i];}
	double mean_Nbr= (double)sum_bridges/num_sim_steps;
	std::cout << "num_sim_steps: " << num_sim_steps << ", mean_Nbr: " << mean_Nbr << std::endl;
	// write Summary file
	std::ofstream outfile;
	outfile.open("../../Summary.txt", std::ios_base::app); // append instead of overwrite
	outfile << atof(argv[1])<< ',' << (mean_Nbr)<< ',' << (sim_end_time) << ',' << (time_attach) << std::endl; 
	outfile.close();

	// -------------------  write RecDensity File -------------------						
	int close_Rec= 0;
	int closeActive_Rec= 0; 
	for (int w=0; w < N_R; w++){// Count the receptors "close" to the Virion
		if(pow(X_r_x[w]-X_CM_x[X_CM_x.size()-1],2)+pow(X_r_y[w]-X_CM_y[X_CM_y.size()-1],2) < R_close*R_close )
			{++close_Rec; if (T_l[w]==1){++closeActive_Rec;}}
	}
	RecDensity << time_cont + (sim_end_time) << ','  << close_Rec/(M_PI*R_close*R_close) << ','  << closeActive_Rec/(M_PI*R_close*R_close) << '\n';
	RecDensity.flush();	
	// --------------------------------------			


    return 0;
}
