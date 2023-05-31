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
#define MATRIX_FLOAT std::vector<std::vector<float>>

#define k_d 0.56
//#define k_on 0.56 //608.478 //5*1000000/(0.62*(4/3)*M_PI*(pow(7.5,3)))
//#define k_off 1.0
//#define k_on 0.56
#define k_off 25.0
#define k_on 0.56

//#define NBrow 10


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



# define L 300.0 //longueur de virus
# define diameter 60.0 //diamÃ¨tre du virus
# define D 7.5 //distance maximale pour les rÃ©actions
// # define MU 0.65 //viscositÃ© de l'eau 
// # define KbT 1.38*pow(10,-23)*308 //constante de Boltzmann et tempÃ©rature (K)
# define D_parallel 3.38*pow(10,6) //coefficients de diffusion
# define D_orthogonal 1.75*pow(10,6)
# define D_THETA 275
# define D_receptor 1.0*pow(10,3) 
# define DT 0.001 //time step
# define N_BRWN 100001 //500000 //nb de steps
# define Taille_Systeme 1000.0
# define N_R 200000//100000 //nombre de rÃ©cepteurs

#define R_Verlet 10.0 //rayon pour la liste de Verlet 




bool isodd(int num) //vÃ©rifie la paritÃ© d'un entier
{
	if(num % 2 == 0)
		return false;
	else
		return true;
}



float ran2(long *idum) //idum doit être un nombre négatif
{
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    
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



float RGV(long *g) //trouve le nombre alÃ©atoire pour le mvt brownien
{
	std::vector<float> random_numbers;
	float p;

    for(int j=0; j<2; j++)
    {
        p = ran2(g);
        random_numbers.push_back(p);
    }
 	float r1 = random_numbers[0];
	float r2 = random_numbers[1];

  	float r3 = sqrt(-2*(log(1-r1)));
  	float r = cos(r2*2*M_PI)*r3;

  	return r;
}



void initialize_IAV(std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<int> *T_l, std::vector<int> *L_to_R, long *idum)
{//Initialise la position des ligands dans le repÃ¨re du CM du virus
	float b = D*2;
	float d = (sqrt(3)*b)/2;
	std::vector<float> y0;
	std::vector<float> y1;
	std::vector<float> y2;



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

	float adjust = abs((y2[0] + y2[y2.size()-1])/2);
   

	for (int a = 0; a < 9; a++)
	{
		if(a==0)
		{
			for (int l = 0; l < y0.size(); l++)
			{
				X_l_x->push_back(-2*b);
				X_l_y->push_back(y0.at(l) + adjust);
				//std::cout << y0.at(l) + adjust << "   " << y0.at(l) << "    " << X_l_y->at(X_l_y->size()-1) <<'\n';
			}
		}
		else if(a==8)
		{
			for (int m = 0; m < y0.size(); m++)
			{
				X_l_x->push_back(2*b);
				X_l_y->push_back(y0.at(m) + adjust);
			}
		}
		else if(isodd(a)==true)
		{
			for (int n = 0; n < y1.size(); n++)
			{
				X_l_x->push_back((a-4)*b/2);
				X_l_y->push_back(y1.at(n) + adjust);
			}
		}
		else if(isodd(a)==false)
		{
			for (int o = 0; o < y2.size(); o++)
			{
				X_l_x->push_back((a-4)*b/2);
				X_l_y->push_back(y2.at(o) + adjust);
			}
		}	
	}

	//POLARIZED
	for (int p = 0; p < X_l_y->size(); p++)
	{
		if(X_l_y->at(p) < (-2*L/5))	
			T_l->push_back(0);
		else
			T_l->push_back(1);
	
	    L_to_R->push_back(-1);
	}

	//RANDOM
	/*for (int p = 0; p < X_l_y->size(); p++)
	{
		T_l->push_back(1);
		L_to_R->push_back(-1);
	}
	int ligand;
	int count = 0;
	while(count < 18)
	{
		ligand = (int)(ran2(idum)*100);
		if(T_l->at(ligand) != 0)
		{
			T_l->at(ligand) = 0;
			count += 1;
		}				
	}*/

	//4 STRIPES
	/*for (int p = 0; p < X_l_y->size(); p++)
	{
		T_l->push_back(1);
		L_to_R->push_back(-1);
	}
	T_l->at(0) = 0;
	T_l->at(24) = 0;
	T_l->at(49) = 0;
	T_l->at(74) = 0;
	T_l->at(98) = 0;
	T_l->at(2) = 0;
	T_l->at(26) = 0;
	T_l->at(51) = 0;
	T_l->at(76) = 0;
	T_l->at(100) = 0;
	T_l->at(4) = 0;
	T_l->at(28) = 0;
	T_l->at(53) = 0;
	T_l->at(78) = 0;
	T_l->at(102) = 0;
	T_l->at(6) = 0;
	T_l->at(30) = 0;
	T_l->at(55) = 0;
	T_l->at(80) = 0;
	T_l->at(104) = 0;
	*/

	//2 STRIPES
	/*for (int p = 0; p < X_l_y->size(); p++)
	{
		T_l->push_back(1);
		L_to_R->push_back(-1);
	}
	T_l->at(11) = 0;
	T_l->at(23) = 0;
	T_l->at(36) = 0;
	T_l->at(48) = 0;
	T_l->at(61) = 0;
	T_l->at(73) = 0;
	T_l->at(86) = 0;
	T_l->at(22) = 0;
	T_l->at(35) = 0;
	T_l->at(47) = 0;
	T_l->at(60) = 0;
	T_l->at(72) = 0;
	T_l->at(85) = 0;
	T_l->at(97) = 0;
	*/

	//1 MIDDLE STRIPE
	/*for (int p = 0; p < X_l_y->size(); p++)
	{
		T_l->push_back(1);
		L_to_R->push_back(-1);
	}
	T_l->at(16) = 0;
	T_l->at(41) = 0;
	T_l->at(66) = 0;
	T_l->at(91) = 0;
	T_l->at(5) = 0;
	T_l->at(29) = 0;
	T_l->at(54) = 0;
	T_l->at(79) = 0;
	T_l->at(103) = 0;
	T_l->at(17) = 0;
	T_l->at(42) = 0;
	T_l->at(67) = 0;
	T_l->at(92) = 0;
	*/

	//3 PATCHES
	/*for (int p = 0; p < X_l_y->size(); p++)
	{
		T_l->push_back(1);
		L_to_R->push_back(-1);
	}
	T_l->at(12) = 0;
	T_l->at(17) = 0;
	T_l->at(24) = 0;
	T_l->at(25) = 0;
	T_l->at(29) = 0;
	T_l->at(30) = 0;
	T_l->at(37) = 0;
	T_l->at(42) = 0;
	T_l->at(77) = 0;
	T_l->at(89) = 0;
	T_l->at(90) = 0;
	T_l->at(101) = 0;
	*/


	//2 PATCHES
	/*for (int p = 0; p < X_l_y->size(); p++)
	{
		T_l->push_back(1);
		L_to_R->push_back(-1);
	}
	T_l->at(24) = 0;
	T_l->at(27) = 0;
	T_l->at(36) = 0;
	T_l->at(37) = 0;
	T_l->at(39) = 0;
	T_l->at(40) = 0;
	T_l->at(49) = 0;
	T_l->at(52) = 0;
	T_l->at(61) = 0;
	T_l->at(62) = 0;
	T_l->at(64) = 0;
	T_l->at(65) = 0;
	T_l->at(74) = 0;
	T_l->at(77) = 0;
	*/

	//4 PATCHES
	/*for (int p = 0; p < X_l_y->size(); p++)
	{
		T_l->push_back(1);
		L_to_R->push_back(-1);
	}
	16,24,28,29,36,37,39,41,49,51,52,56,64,68,69,81
	T_l->at(16) = 0;
	T_l->at(24) = 0;
	T_l->at(28) = 0;
	T_l->at(29) = 0;
	T_l->at(36) = 0;
	T_l->at(37) = 0;
	T_l->at(39) = 0;
	T_l->at(41) = 0;
	T_l->at(49) = 0;
	T_l->at(51) = 0;
	T_l->at(52) = 0;
	T_l->at(56) = 0;
	T_l->at(64) = 0;
	T_l->at(68) = 0;
	T_l->at(69) = 0;
	T_l->at(81) = 0;
	*/
}



void initialize_SA(std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *T_r, long *m)
{//initialise les rÃ©cpteurs de maniÃ¨re alÃ©atoire dans le plan
	float x;
	float y;
	int check;
	
 	x = ran2(m)*Taille_Systeme - Taille_Systeme/2;
	y = ran2(m)*Taille_Systeme - Taille_Systeme/2;
		
	X_r_x->push_back(x);
    X_r_y->push_back(y);	

	while(X_r_x->size() < N_R) 
	{
		x = ran2(m)*Taille_Systeme - Taille_Systeme/2;
		y = ran2(m)*Taille_Systeme - Taille_Systeme/2;
	    check = 0;

		X_r_x->push_back(x);
		X_r_y->push_back(y);  
	}

	for (int j = 0; j < X_r_x->size(); j++)
	{
		T_r->push_back(1);
	}
}




bool check_breakage(std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *T_l, std::vector<int> *L_to_R, float DX, float DY, float DTHETA)
{//vÃ©rifie si la diffusion ne brise pas un pont
	bool bk = false;
	int receptor = 0;
    float theta;
    float x;
    float y;

    theta = DTHETA + X_CM_theta->at((X_CM_theta->size()) - 1);
	//std::cout << "overextended bridges : ";
	for (int i = 0; i < X_l_x->size(); ++i)
	{
		if (T_l->at(i) == 2)
		{
			receptor = L_to_R->at(i);
		    
		    x = (X_l_x->at(i))*cos(theta) - (X_l_y->at(i))*sin(theta) + X_CM_x->at((X_CM_x->size()) -1) ;
		    y = (X_l_x->at(i))*sin(theta) + (X_l_y->at(i))*cos(theta) + X_CM_y->at((X_CM_y->size()) -1) ;
		    
		    
			if (pow(x + DX - X_r_x->at(receptor - 1),2) + pow(y + DY - X_r_y->at(receptor - 1),2) > D*D)
			{
				//std::cout << "check_bk: " << L_to_R->at(i) << " ";
				//std::cout << pow(x + DX - X_r_x->at(receptor - 1),2) + pow(y + DY - X_r_y->at(receptor - 1),2) << '\n';
				bk = true;
			}
		}
	}
	
	return bk;
}



int diffusion(std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *T_l, std::vector<int> *L_to_R, long *idum, float delta_t)
{//effectue la diffusion aprÃ¨s avoir vÃ©rifiÃ© si elle est autorisÃ©e
    std::vector<float> random_numbers;
    int acc;

    for(int k = 0; k<3; k++)
    {
        random_numbers.push_back(RGV(idum));
    }
    
    float old_theta = X_CM_theta->at((X_CM_theta->size()) - 1);
    
	float DX_Brwn = (-1)*sqrt(2*D_parallel*delta_t)*sin(old_theta)*random_numbers[0] + sqrt(2*D_orthogonal*delta_t)*cos(old_theta)*random_numbers[1];
	// version de Fletcher: (sqrt(2*D_parallel*delta_t)*pow(cos(old_theta),2) + sqrt(2*D_orthogonal*delta_t)*pow(sin(old_theta),2))*random_numbers[0] + ((sqrt(2*D_parallel*delta_t) - sqrt(2*D_orthogonal*delta_t))*cos(old_theta)*sin(old_theta))*random_numbers[1];
	float DY_Brwn = sqrt(2*D_parallel*delta_t)*cos(old_theta)*random_numbers[0] + sqrt(2*D_orthogonal*delta_t)*sin(old_theta)*random_numbers[1]; 
	// version de Fletcher: ((sqrt(2*D_parallel*delta_t) - sqrt(2*D_orthogonal*delta_t))*cos(old_theta)*sin(old_theta))*random_numbers[0] + (sqrt(2*D_parallel*delta_t)*pow(sin(old_theta),2) + sqrt(2*D_orthogonal*delta_t)*pow(cos(old_theta),2))*random_numbers[1];
	float DTHETA_Brwn = sqrt(2*D_THETA*delta_t)*random_numbers[2];
	
	bool bk;
	

    
    bk = check_breakage(X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x, X_r_y, T_l, L_to_R, DX_Brwn, DY_Brwn, DTHETA_Brwn);
        
        if(bk == false) //on effectue la diffusion
        {
        	X_CM_x->push_back(X_CM_x->at((X_CM_x->size()) -1) + DX_Brwn);
        	X_CM_y->push_back(X_CM_y->at((X_CM_y->size()) -1) + DY_Brwn);
        	X_CM_theta->push_back(X_CM_theta->at((X_CM_theta->size()) -1) + DTHETA_Brwn);
            //std::cout << "check breakage false " << '\n'; 
            acc = 1;
        }   
        
        if(bk == true) //le systÃ¨me ne bouge pas
        {
            X_CM_x->push_back(X_CM_x->at(X_CM_x->size() - 1));
        	X_CM_y->push_back(X_CM_y->at(X_CM_y->size() - 1));
        	X_CM_theta->push_back(X_CM_theta->at(X_CM_theta->size() - 1));
            //std::cout << "check breakage true" << '\n'; 
        	acc = 0;
        }
    return acc;
}

/*void update_position_receptor(float rx, float ry)
{//conditions aux bords
	float size = Taille_Systeme;
	float r_x = float(rx);
	float r_y = float(ry);
	if(r_x > size/2)
		rx = r_x - size;
	if(r_x < -1*size/2)
		rx = r_x + size;

	if(r_y > size/2)
		ry = r_y - size;
	if(r_y < -1*size/2)
		ry = r_y + size;

	std::cout << "inside: " << r_x << ", " << r_y << "\n";
}*/

bool check_breakage_receptors(std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<float> *X_l_x, std::vector<float> *X_l_y, float r_x, float r_y, std::vector<int> *T_l, std::vector<int> *L_to_R, std::vector<int> *list_bridged_receptors, float DX, float DY, int rec)
{//vÃ©rifie si la diffusion  des recepteurs ne brise pas un pont
	//std::cout << "start check bk rec" << rec << '\n';
	bool bk_rec;
	int receptor = 0;
    float theta;
    float x;
    float y;
    int i = -1;
    
    theta = X_CM_theta->at((X_CM_theta->size()) - 1);
	//std::cout << "overextended bridges : ";
    /*std::cout << "L_to_R:";
	for(int j = 0; j < L_to_R->size(); j++)
	{

			std::cout << L_to_R->at(j) << ' ';
		
	}
	std::cout << '\n';*/

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
	//std::cout << "before if: " << pow(x - (r_x + DX),2) + pow(y - (r_y + DY),2) << " " << D*D << " ----- receptor: " << receptor << '\n';
    
    
	if (pow(x - (r_x + DX),2) + pow(y - (r_y + DY),2) > D*D)
	{
		//std::cout << "in if: " << pow(x - (r_x + DX),2) + pow(y - (r_y + DY),2) << " " << D*D << "\n";
		bk_rec = true;
	}

	//std::cout << "end check bk rec: " << '\n';	
	return bk_rec;
}


int diffusion_receptors(std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *T_l, std::vector<int> *L_to_R, std::vector<int> *list_bridged_receptors, long *idum, float delta_t)
{//effectue la diffusion aprÃ¨s avoir vÃ©rifiÃ© si elle est autorisÃ©e
    //std::cout << "start diffusion_receptors" << "\n";
    std::vector<float> random_numbers;
    int acc;
    //std::cout << "\n" << "start DR" << "\n";
    
    for(int k = 0; k<2; k++)
    {
        	random_numbers.push_back(0.0);
    }

    for (int rec = 0; rec < N_R; ++rec)
	{	
		//std::cout << rec << '\n';
    	random_numbers[0] = RGV(idum);
    	random_numbers[1] = RGV(idum);
       	//std::cout << random_numbers[0] << ", " << random_numbers[1] << "\n"; 
    
		float DX_Brwn = sqrt(2*D_receptor*delta_t)*random_numbers[0];
		float DY_Brwn = sqrt(2*D_receptor*delta_t)*random_numbers[1];

		//std::cout << "=> " << DX_Brwn << ", " << DY_Brwn << "\n";
		/*std::cout << "list_bridged: ";
		for(int j = 0; j < list_bridged_receptors->size(); j++)
    	{

    			std::cout <<  list_bridged_receptors->at(j) << ' ';
    		
    	}
    	std::cout << '\n';*/
	
    	bool bridged = false;
    	for(int j = 0; j < list_bridged_receptors->size(); j++)
    	{
    		if(list_bridged_receptors->at(j) - 1 == rec)
    		{
    			//std::cout << "bridged: " << list_bridged_receptors->at(j) << '\n';
    			bridged = true;
    		}
    	}

    	if(bridged == true)
		{
			bool bk_rec = false;
			bk_rec = check_breakage_receptors(X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x->at(rec), X_r_y->at(rec), T_l, L_to_R, list_bridged_receptors, DX_Brwn, DY_Brwn, rec);
			//std::cout << "bk_rec: " << bk_rec << "\n";
			if(bk_rec == true)
			{
				acc = 0;
			}
			if(bk_rec == false)
			{
				X_r_x->at(rec) = X_r_x->at(rec) + DX_Brwn;
	       		X_r_y->at(rec) = X_r_y->at(rec) + DY_Brwn;
	        	
	            	//std::cout << "check breakage false " << '\n'; 
	           	acc = 1;
	           	float size = Taille_Systeme;
				if(X_r_x->at(rec) > size/2)
					X_r_x->at(rec) = X_r_x->at(rec) - size;
				if(X_r_x->at(rec) < -1*size/2)
					X_r_x->at(rec) = X_r_x->at(rec) + size;

				if(X_r_y->at(rec) > size/2)
					X_r_y->at(rec) = X_r_y->at(rec) - size;
				if(X_r_y->at(rec) < -1*size/2)
					X_r_y->at(rec) = X_r_y->at(rec) + size;
	           	//update_position_receptor(X_r_x->at(rec),X_r_y->at(rec));
	        }
		}
		else
		{
			X_r_x->at(rec) = X_r_x->at(rec) + DX_Brwn;
	       	X_r_y->at(rec) = X_r_y->at(rec) + DY_Brwn;
	       	//std::cout << "\n" << "before: " << X_r_x->at(rec) << ", " << X_r_y->at(rec) << "\n";

	       	float size = Taille_Systeme;
			if(X_r_x->at(rec) > size/2)
				X_r_x->at(rec) = X_r_x->at(rec) - size;
			if(X_r_x->at(rec) < -1*size/2)
				X_r_x->at(rec) = X_r_x->at(rec) + size;

			if(X_r_y->at(rec) > size/2)
				X_r_y->at(rec) = X_r_y->at(rec) - size;
			if(X_r_y->at(rec) < -1*size/2)
				X_r_y->at(rec) = X_r_y->at(rec) + size;
	       	//update_position_receptor(&X_r_x->at(rec),&X_r_y->at(rec));
	       	//std::cout << "after: " << X_r_x->at(rec) << ", " << X_r_y->at(rec) << "\n";
		}
		    	
    } 
    //std::cout << "end diffusion_receptors" << "\n";
    return acc;
}



int Cell_to_Index(int N_x, int N_y, int ix, int iy) //reÃ§oit les coordonnÃ©es de la cellule et renvoie l'indice de la cellule
{
	int ixp;
	int iyp;
	int nc;

	ixp= ix;
	iyp = iy;

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



void Initialize_L_L(MATRIX_INT *L_L, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> X_r_x, std::vector<float> X_r_y, std::vector<float> X_CM_x, std::vector<float> X_CM_y, std::vector<float> X_CM_theta, MATRIX_INT Neigh, std::vector<int> hoc_R, std::vector<int> ll_R, std::vector<int> T_l, std::vector<int> L_to_R, int N_x, int N_y, float R_Cell, int N_L, int step)
{//Initialise L_L en utilisant les listes de voisinage
    float theta;
    float x, y;
	for (int alpha = 0; alpha < N_L; ++alpha)
	{ 
		//std::cout << "test LL, ligand " << alpha << '\n';
		theta = X_CM_theta.at((X_CM_theta.size()) - 1);
	    x = (X_l_x->at(alpha))*cos(theta) - (X_l_y->at(alpha))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
		y = (X_l_x->at(alpha))*sin(theta) + (X_l_y->at(alpha))*cos(theta) + X_CM_y[X_CM_y.size() - 1];

		int ix = int((x + Taille_Systeme/2)/R_Cell) + 1;
		int iy = int((y + Taille_Systeme/2)/R_Cell) + 1;
		int index = Cell_to_Index(N_x, N_y, ix, iy);
        
		for (int j = 0; j < 9; ++j)
		{
			//if(alpha == 54)
			//	std::cout << "in LL loop " << j << "\n";
			int i = Neigh[index-1][j];
			
			if(hoc_R[i] != 0)
			{
			//	if(alpha == 54)
			//		std::cout << "in LL loop test 2" << "\n";
	            if(pow(x - X_r_x[hoc_R[i] - 1],2) + pow(y - X_r_y[hoc_R[i] - 1],2) < R_Verlet*R_Verlet)
    				(L_L->at(alpha)).push_back(hoc_R[i]); 
				int g = ll_R[hoc_R[i]];

			//	if(alpha == 54)
			//		std::cout << "in LL loop test 3: " << hoc_R[i] << " and " << ll_R[hoc_R[i]] << "\n";

				while(g != 0)
				{
					//if(alpha == 54)
						//std::cout << "in LL loop test 4: " << g << "\n";
	                if(pow(x - X_r_x[g - 1],2) + pow(y - X_r_y[g - 1],2) < R_Verlet*R_Verlet)
					    (L_L->at(alpha)).push_back(g);
					int k = g;	
					g = ll_R[k];
				}
			}
		}
		//std::cout << "test LL 2" << '\n';
		if(T_l[alpha] == 2 && pow(x - X_r_x[L_to_R[alpha] - 1],2) + pow(y - X_r_y[L_to_R[alpha] - 1],2) > R_Verlet*R_Verlet) //gère l'exception due à l'imprécision dans check_breakage
			(L_L->at(alpha)).push_back(L_to_R[alpha]);

		//std::cout << "test LL 3" << '\n';
		int N = (L_L->at(alpha)).size();
		(L_L->at(alpha)).insert((L_L->at(alpha)).begin(), N);
	}
}



void Initialize_L_R(MATRIX_INT *L_R, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> X_r_x, std::vector<float> X_r_y, std::vector<float> X_CM_x, std::vector<float> X_CM_y, std::vector<float> X_CM_theta, MATRIX_INT Neigh, std::vector<int> hoc_L, std::vector<int> ll_L, std::vector<int> T_r, std::vector<int> T_l, std::vector<int> L_to_R, int N_x, int N_y, float R_Cell, int N_L, int step)
{//initialise L_R en utilisant les listes de voisinage
    float theta;
    float x, y;
    float count = 0;

    //std::cout << "un" << '\n'; 

	for (int alpha = 0; alpha < N_R; ++alpha)
	{ 
		//std::cout << "in LR loop " << alpha << ", LR1" << "\n";
		int ix = int((X_r_x[alpha] + Taille_Systeme/2)/R_Cell) + 1;
		int iy = int((X_r_y[alpha] + Taille_Systeme/2)/R_Cell) + 1;
		int index = Cell_to_Index(N_x, N_y, ix, iy);
        //std::cout << "LR2" << "\n";
		for (int j = 0; j < 9; ++j)
		{
			//if(alpha == 1059)
				//std::cout << "in LR 2nd loop " << j << ", LR3" << "\n";
			int i = Neigh[index - 1][j];
			
			if(hoc_L[i] != 0)
			{
				//if(alpha == 1059)
				//	std::cout << "in LR 2nd loop " << j << ", LR4 - hoc[i] = " << hoc_L[i] << "\n";
				theta = X_CM_theta.at((X_CM_theta.size()) - 1);
	    		x = (X_l_x->at(hoc_L[i] - 1))*cos(theta) - (X_l_y->at(hoc_L[i] - 1))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
				y = (X_l_x->at(hoc_L[i] - 1))*sin(theta) + (X_l_y->at(hoc_L[i] - 1))*cos(theta) + X_CM_y[X_CM_y.size() - 1];

	            if(pow(x - X_r_x[alpha],2) + pow(y - X_r_y[alpha],2) < R_Verlet*R_Verlet)
    				(L_R->at(alpha)).push_back(hoc_L[i]); 
		    	
		    	//std::cout << alpha << '\n'; 			
				
				int g = ll_L[hoc_L[i]];
				
				//if(alpha == 1059)
				//	std::cout << "in LR 2nd loop " << j << ", LR5 - ll_L[hoc_L[i]] = " << ll_L[hoc_L[i]] << "\n";				
				count = 0;
				while(g != 0)
				{
					//if(alpha == 1059)
					//	std::cout << "in LR 2nd loop " << j << ", LR6 - g = " << g << "\n";
					//std::cout << count << " " << ll_L.size() << '\n';
					theta = X_CM_theta.at((X_CM_theta.size()) - 1);
				    x = (X_l_x->at(g-1))*cos(theta) - (X_l_y->at(g-1))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
					y = (X_l_x->at(g-1))*sin(theta) + (X_l_y->at(g-1))*cos(theta) + X_CM_y[X_CM_y.size() - 1];

	                if(pow(x - X_r_x[alpha],2) + pow(y - X_r_y[alpha],2) < R_Verlet*R_Verlet)
					    (L_R->at(alpha)).push_back(g);
					int k = g;	
					g = ll_L[k];
					count += 1;
				}
				//if(alpha == 1059)
					//std::cout << "in LR 2nd loop - out" << "\n";	

			}
			//if(alpha == 1059)
				//std::cout << "in LR 2nd loop - out2" << "\n";
		}
		//if(alpha == 1059)
			//std::cout << "after LR 2nd loop" << "\n";

		
		if(T_r[alpha] == 2) //corrige l'exception due à l'imprécision numérique dans check_breakage
		{
			//std::cout << "are we bridged ? - check for receptor " << alpha << "\n";
			for(int z = 0; z < N_L; z++)
			{
				//std::cout << "LR in if bridged - ligands: " << z << "\n";
				//std::cout << "T_l = " << T_l[z] << " and L_to_R = " << L_to_R[z] << '\n';
				if(T_l[z] == 2 && L_to_R[z] == alpha + 1)
				{
				//	std::cout << "do we go here ?" << "\n";
					theta = X_CM_theta.at((X_CM_theta.size()) - 1);
					x = (X_l_x->at(z))*cos(theta) - (X_l_y->at(z))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
					y = (X_l_x->at(z))*sin(theta) + (X_l_y->at(z))*cos(theta) + X_CM_y[X_CM_y.size() - 1];
					if(pow(x - X_r_x[alpha],2) + pow(y - X_r_y[alpha],2) > R_Verlet*R_Verlet)
						(L_R->at(alpha)).push_back(z+1);
				}
			}
		}

		int N = (L_R->at(alpha)).size();
		(L_R->at(alpha)).insert((L_R->at(alpha)).begin(), N);
	}
}


void Initialize_Affinity(MATRIX_FLOAT *a, MATRIX_INT *L_L, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *L_to_R, std::vector<int> T_l, std::vector<int> T_r, float k, int N_L, int t1, int t2)
{//initialise les matrices a_d, a_on et a_off
	float a_sum = 0;
	float theta, x, y;

	for (int alpha = 0; alpha < N_L; ++alpha)
	{
		theta = X_CM_theta->at((X_CM_theta->size()) - 1);
	    x = (X_l_x->at(alpha))*cos(theta) - (X_l_y->at(alpha))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
		y = (X_l_x->at(alpha))*sin(theta) + (X_l_y->at(alpha))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);
	
		a_sum = 0;
		for (int i = 1; i < (L_L->at(alpha)).size(); ++i)
		{
			//if(alpha == 80 && L_L->at(alpha)[i] == 124029)
				//std::cout << " -- distance in Initialize_Affinity: " << pow(x - X_r_x->at(124029-1),2) + pow(y - X_r_y->at(124029-1),2) << " while max = " << D*D << "\n";
 
			if((pow(x - X_r_x->at(L_L->at(alpha)[i] - 1),2) + pow(y - X_r_y->at(L_L->at(alpha)[i] - 1),2) < D*D))
			{
				//if(alpha == 81 && L_L->at(alpha)[i] == 124029)
					//std::cout << "we go in the if" << "\n";
		        if(T_l[alpha] == t1 && T_r[((L_L->at(alpha)).at(i)) - 1] == t2)
		        {
	    	       	if(t1 == 2 && t2 == 2) //gère le cas où on a un pont
	    	       	{
	    	       	    if(L_to_R->at(alpha) == L_L->at(alpha)[i])
	    	       	    {
	    	       	        a_sum += k;
			                (a->at(alpha)).push_back(k);
	    	       	    }
	    	       	    else
	    	       	    {
		                    (a->at(alpha)).push_back(0);
	    	       	    }
	    	       	}
	    	       	else //cas où on a pas de pont
	    	       	{
			            a_sum += k;
			            (a->at(alpha)).push_back(k);
	    	       	}
		        }
		        else
		        {
		            (a->at(alpha)).push_back(0);
			    }
			}
			else
			{
				(a->at(alpha)).push_back(0);
			}    
		}
		(a->at(alpha)).insert((a->at(alpha)).begin(), a_sum);
	}
}



std::vector<int> Index_to_Cell(int N_x, int N_y, int nc) //reÃ§oit l'indice de la cellule et renvoie ses coordonnÃ©es
{
	std::vector<int> ixy;
	int ix;
	int iy;

	iy = (nc-1)/N_x + 1;
	ix = nc - (iy - 1)*N_x;

	ixy.push_back(ix);
	ixy.push_back(iy);

	return ixy;
}



MATRIX_INT New_List(float R_Cell, int N_Cell, int N_part, int N_x, int N_y, std::vector<float> *X_part_x, std::vector<float> *X_part_y)
{//crÃ©e les listes hoc, ll et bl
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


MATRIX_INT Set_Neigh(int N_Cell, int N_x, int N_y) //crÃ©e les listes de voisinage
{
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
		for (int icx = ix-1; icx < ix + 2; ++icx)
		{
			for (int icy = iy-1 ; icy < iy + 2; ++icy)
			{
				int index = Cell_to_Index(N_x,N_y,icx,icy);
				Neigh[i][cnt] = index;
				cnt=cnt+1;
			}
		}
	}
	return Neigh;
}





void UpGrade(int Old_Cell, int New_Cell, int ipar, std::vector<int> *ll, std::vector<int> *bl, std::vector<int> *hoc)
{//met Ã  jour les cellules si des particules rentrent ou sortent
	//supprimer particule de l'ancienne cellule
	if (ll->at(ipar) != 0)
	{
		if (bl->at(ipar) != 0) //particule au milieu de la cellule
		{
			ll->at(bl->at(ipar)) = ll->at(ipar);
			bl->at(ll->at(ipar)) = bl->at(ipar);
		}
		else //particule est la premiÃ¨re de la cellule
		{
			bl->at(ll->at(ipar)) = 0;
			hoc->at(Old_Cell) = ll->at(ipar);
		}
	}
	else
	{
		if (bl->at(ipar) != 0) //particule est la derniÃ¨re de la celulle
		{
			ll->at(bl->at(ipar)) = 0;
		}
		else //particule est la seule dans la cellule
		{
			hoc->at(Old_Cell) = 0;
		}
	}


	//ajouter particule Ã  la nouvelle cellule
	bl->at(ipar) = 0;
	ll->at(ipar) = hoc->at(New_Cell);
	if (hoc->at(New_Cell) !=0)
	{
		bl->at(hoc->at(New_Cell)) = ipar;
	}
	hoc->at(New_Cell) = ipar;
}




int Choose(std::vector<float> V, long *g) //choisit une rÃ©action, un ligand ou un rÃ©cepteur
{
	float r = ran2(g); 
	float sumpar = 0.0;
	int i = 0;

        if(V[0] == 0)
        return -1;
        
	while(sumpar <= r && i<V.size()-1) 
        { 
		sumpar += (V[i+1]/V[0]);
                i += 1;
	}
	
	return i;
}


int Get_Id(std::vector<int> L_L_line, int i) //donne la position du rÃ©cepteur i dans le vecteur L_L[alpha] (alpha in [0,N_L-1], i in [1,N_R])
{
    int id = 1;
    for(int w = 0; w < (L_L_line.size()); w++)
    {
        if(L_L_line[w + 1] == i)
        {
            return id;
        }
        else
        {
            id += 1;   
        }
    }
    return -1; //si erreur
}


std::vector<float> Get_a_X_tot(MATRIX_FLOAT *a_X, int N_L)  //crÃ©e le vecteur nÃ©cessaire pour que Choose puisse choisir un ligand une fois la rÃ©action choisie
{
	std::vector<float> a_X_tot;
	a_X_tot.clear();
	float a_X_sum = 0.0;

	for (int i = 0; i < N_L; ++i)
	{
		a_X_tot.push_back((a_X->at(i)).at(0));
		a_X_sum += (a_X->at(i)).at(0);
	}
	a_X_tot.insert(a_X_tot.begin(), a_X_sum);
	return a_X_tot;
}


std::vector<float> Get_a_tot(std::vector<float> *a_d_tot, std::vector<float> *a_on_tot, std::vector<float> *a_off_tot, int N_L) 
{//crÃ©e le vecteur nÃ©cessaire pour choisir une rÃ©action
	float a_tot;
	float a_d = 0.0;
	float a_on = 0.0;
	float a_off = 0.0;
	float epsilon = 0.0000001;
	std::vector<float> a;
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
		if(a_on < epsilon) //Ã  remplacer par k_on quand on remet k_on !=0
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


float Get_time_reaction(float a_tot, long *idum)
{
    float t;
    t =  (-1)*(log(1 - ran2(idum)))/a_tot;
    
    return t;
}

std::vector<float> center(std::vector<float> w1,std::vector<float> w2, float R)
{//Return the center of a circle of radius R supporting v1 and v2
	//std::cout << "A" << '\n';
	float dist;
	std::vector<float> m, w12, u, c;

	m.clear();
	w12.clear();
	u.clear();
	c.clear();

	for(int i = 0; i<w1.size();i++)
	{
		m.push_back((w1[i]+w2[i])/2);
		w12.push_back(w2[i] - w1[i]);
	}	

	u.push_back(-w12[1]);
	u.push_back(w12[0]);

	float u0_old = 0;
	u0_old = u[0];
	u[0] = u[0]/(sqrt(u[0]*u[0] + u[1]*u[1]));
	u[1] = u[1]/(sqrt(u0_old*u0_old + u[1]*u[1]));

	if((w12[0]*u[1])-(w12[1]*u[0]) < 0)
	{
		u[0] = -u[0];
		u[1] = -u[1];
	}
	dist = sqrt(R*R - 0.25*w12[0]*w12[0] - 0.25*w12[1]*w12[1]);

	c.push_back(m[0]+dist*u[0]);
	c.push_back(m[1]+dist*u[1]);

	return c;
}

float distance2(std::vector<float> r1,std::vector<float> r2)
{//square distance between two vectors 
	float dist2;


	dist2=pow((r1[0]-r2[0]),2)+pow((r1[1]-r2[1]),2);

	return dist2;
}

std::vector<float> get_position_linked_ligand(std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<int> *L_to_R, int receptor, float theta_mod)
{//trouve la position du ligands lié au récepteur i, dans le repère du laboratoire
	std::vector<float> position;
	position.clear();
	int ligand = -1;
	float theta;

	for(int j = 0; j < L_to_R->size(); j++)
	{
		if(L_to_R->at(j) == receptor)
			ligand = j + 1;
	}

	if(ligand == -1)
	{
		position.push_back(X_CM_x->at(X_CM_x->size() - 1));
		position.push_back(X_CM_y->at(X_CM_y->size() - 1));

		return position;
	}
	
	theta = X_CM_theta->at((X_CM_theta->size()) - 1) + theta_mod;

	position.push_back((X_l_x->at(ligand-1))*cos(theta) - (X_l_y->at(ligand-1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1));
	position.push_back((X_l_x->at(ligand-1))*sin(theta) + (X_l_y->at(ligand-1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1));

	/*if(repere == 0) //dans le repère du laboratoire
	{	theta = X_CM_theta->at((X_CM_theta->size()) - 1);

		position.push_back((X_l_x->at(ligand-1))*cos(theta) - (X_l_y->at(ligand-1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1));
		position.push_back((X_l_x->at(ligand-1))*sin(theta) + (X_l_y->at(ligand-1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1));
	}

	if(repere == 1) //dans le repère du CM
	{
		position.push_back((X_l_x->at(ligand-1)));
		position.push_back((X_l_y->at(ligand-1)));
	}*/

	return position;
}

std::vector<float> get_new_position_linked_ligand(std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> r_trial, std::vector<float> *X_CM_theta, std::vector<int> *L_to_R, int receptor)
{//trouve la position du ligands lié au récepteur i, dans le repère du laboratoire
	std::vector<float> position;
	position.clear();
	int ligand = -1;
	float theta;

	for(int j = 0; j < L_to_R->size(); j++)
	{
		if(L_to_R->at(j) == receptor)
			ligand = j + 1;
	}

	if(ligand == -1)
	{
		std::cout << "WARNING get_new_position_linked_ligand" << "\n";
		position.push_back(r_trial[0]);
		position.push_back(r_trial[1]);

		return position;
	}

	theta = X_CM_theta->at(X_CM_theta->size() - 1);

	position.push_back((X_l_x->at(ligand-1))*cos(theta) - (X_l_y->at(ligand-1))*sin(theta) + r_trial[0]);
	position.push_back((X_l_x->at(ligand-1))*sin(theta) + (X_l_y->at(ligand-1))*cos(theta) + r_trial[1]);

	return position;
}


std::vector<int> update_GCH(std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<int> *L_to_R, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *list_bridged_receptors, std::vector<int> *list_CB, int irc)
{
	//std::cout << "C" << '\n';
	float dwc2;
	int f,already_inserted;

	std::vector<float> wp;
	std::vector<float> w_test;
	std::vector<int> new_list_CB;

	std::vector<float> wp1;
	std::vector<float> wp2;
	std::vector<float> w1;
	std::vector<float> w2;

	std::vector<float> c;
	std::vector<float> X_ligand;

	X_ligand.clear();
	X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_bridged_receptors->at(irc),0);

	//std::cout << "position recepteur : " << X_r_x->at(list_bridged_receptors->at(irc) - 1) << " " << X_r_y->at(list_bridged_receptors->at(irc) - 1) << "\n";
	wp.push_back(X_r_x->at(list_bridged_receptors->at(irc) - 1)  - X_ligand[0]);
	wp.push_back(X_r_y->at(list_bridged_receptors->at(irc) - 1)  - X_ligand[1]);

	w_test = wp;

	already_inserted = 0;

	//std::cout << "C2" << '\n';
	/*std::cout << "list CB before for loop in update_GCH: ";
	for(int k = 0; k<list_CB->size();k++)
	{
		std::cout << list_CB->at(k) << " ";
	}
	std::cout << "\n" << "\n";*/
	//std::cout << '\n' << "list_CB->size(): " << list_CB->size() << "\n";
	for(int i_Hull = 0;i_Hull<list_CB->size();i_Hull++)
	{
		//std::cout << '\n' << "i_Hull: " << i_Hull << "\n";

		wp1.clear();
		wp2.clear();
		w1.clear();
		w2.clear();
		X_ligand.clear();

		//std::cout << "C2 1" << '\n';	

		X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_CB->at(i_Hull),0);
		//std::cout << "recepteur wp1: " << list_CB->at(i_Hull) << ", et CB position: " << i_Hull << "\n";
		//std::cout << "position recepteur wp1: " << X_r_x->at(list_CB->at(i_Hull) - 1) << " " << X_r_y->at(list_CB->at(i_Hull) - 1) << "\n" << "\n";


		wp1.push_back(X_r_x->at(list_CB->at(i_Hull) - 1) - X_ligand[0]);
		wp1.push_back(X_r_y->at(list_CB->at(i_Hull) - 1) - X_ligand[1]);

		X_ligand.clear();
		X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_CB->at((i_Hull+1) % list_CB->size()),0);
		//std::cout << "recepteur wp2: " << list_CB->at((i_Hull + 1) % list_CB->size()) << ", et CB position: " << (i_Hull + 1) % list_CB->size() << "\n";
		//std::cout << "position recepteur wp2: " << X_r_x->at(list_CB->at((i_Hull + 1) % list_CB->size()) - 1) << " " << X_r_y->at(list_CB->at((i_Hull + 1) % list_CB->size()) - 1) << "\n";

		wp2.push_back(X_r_x->at(list_CB->at((i_Hull + 1) % list_CB->size()) - 1) - X_ligand[0]);
		wp2.push_back(X_r_y->at(list_CB->at((i_Hull + 1) % list_CB->size()) - 1) - X_ligand[1]);		

		//std::cout << "C2 2" << '\n';

		w1 = wp1;
		w2 = wp2;

		/*std::cout << "w1: ";
		for(int i=0;i<2;i++)
		{
			std::cout << w1[i] << " ";
		}
		std::cout << "\n" << "w2: ";
		for(int i=0;i<2;i++)
		{
			std::cout << w2[i] << " ";
		}*/


		c = center(w1,w2,D);

		//std::cout << "\n" << "center: ";
		/*for(int i=0; i<2; i++)
		{
			std::cout << c[i] << " ";
		}*/
		//std::cout << "\n";

		dwc2 = distance2(c,w_test);

		//std::cout << "dwc2 = " << dwc2 << '\n';
		//std::cout << "D*D = " << D*D << '\n';


		//déterminer si w_test est inclus ou non
		if(dwc2>=D*D) //w_test not contained
		{
			//std::cout << "not contained" << "\n";
			f = 1;  
		}
		if(dwc2<D*D) //w_test contained
		{
			//std::cout << "contained" << "\n";
			f = 0;
		}

		//std::cout << "C4" << '\n';

		//compléter la liste des CB en fonction de ça
		if(f == 1) //not contained
		{
			//std::cout << "C4 1" << '\n';
			if(already_inserted == 0)
			{
				//std::cout << "inserted: " << list_bridged_receptors->at(irc) << "\n";
				new_list_CB.push_back(list_bridged_receptors->at(irc));
				already_inserted = 1;	
			}
		}
		if(f == 0) //contained
		{
			//std::cout << "C4 2" << '\n';
			if(i_Hull == 0)
			{
				//std::cout << "C4 3" << '\n';	
				//std::cout << list_CB->at(i_Hull) << "\n";
				new_list_CB.push_back(list_CB->at(i_Hull));
			}
			else
			{
				//std::cout << "C4 4" << '\n';	
				if(new_list_CB[new_list_CB.size()-1] != list_CB->at(i_Hull))	
				{	
					//std::cout << "C4 5" << '\n';		
					//std::cout << list_CB->at(i_Hull) << "\n";
					new_list_CB.push_back(list_CB->at(i_Hull));
				}	
			}

			if(new_list_CB[0] != list_CB->at((i_Hull + 1) % list_CB->size()))
			{
				//std::cout << "C4 6" << '\n';
				//std::cout << list_CB->at((i_Hull+1) % list_CB->size()) << "\n";
				new_list_CB.push_back(list_CB->at((i_Hull+1) % list_CB->size()));
			}
		}
	}

	//std::cout << "C5" << '\n';


	//ERRORS
	std::vector<float> wp_iir;
	std::vector<float> w_iir;
	float dist, rc2;

	if(new_list_CB.size() == 1)
	{
		//std::cout << "if 1 CB in new_CB" << '\n';

		if(new_list_CB[0] != list_bridged_receptors->at(irc))
			std::cout << "some problems when new CB contains a single point" << "\n";
		
		dist = 0;
		rc2 = -1;
		//std::cout << list_CB->size() << "\n";
		for(int iir=0; iir<list_CB->size(); iir++)
		{
			//std::cout << "iir: " << iir << '\n';
			wp_iir.clear();

			X_ligand.clear();
			X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_CB->at(iir),0);
			
			wp_iir.push_back(X_r_x->at(list_CB->at(iir) - 1) - X_ligand[0]);
			wp_iir.push_back(X_r_y->at(list_CB->at(iir) - 1) - X_ligand[1]);

			w_iir = wp_iir;

			float dist_test = distance2(w_test,w_iir);
			//std::cout << "dist_test: " << dist_test << "\n";
			//std::cout << "distances to compare: " << dist_test << ", " << dist << '\n';

			if(dist_test > dist)
			{
				dist = dist_test;
				//std::cout << list_CB->at(iir) << "\n";
				rc2 = list_CB->at(iir);
			}		
			//iir = iir+1;		
		}
		//std::cout << "inserted in 1CB if: " << rc2 << "\n";
		new_list_CB.push_back(rc2);
		if(rc2 == -1)
			std::cout << "problem in Update GCH" << "\n";
	}

	return new_list_CB;
}

std::vector<int> gen_GCH(std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<int> *L_to_R, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *list_bridged_receptors, std::vector<int> *list_CB)
{
	//std::cout << "D" << '\n';
	//std::cout << "gen_GCH used" << "\n";
	std::vector<int> new_list_CB;

	list_CB->clear();
	list_CB->push_back(list_bridged_receptors->at(0));
	list_CB->push_back(list_bridged_receptors->at(1));

	//std::cout << "list CB in beginning gen_GCH: " << list_CB->at(0) << " " << list_CB->at(1) << '\n';

	int i;
	i = 2;

	while(i<list_bridged_receptors->size())
	{
		//std::cout << '\n' << "WHILE: i = " << i << '\n';
		new_list_CB.clear();
		new_list_CB = update_GCH(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, X_r_x, X_r_y, list_bridged_receptors, list_CB, i);
		list_CB->clear();
		//std::cout << "new_list_CB->size() after update: " << new_list_CB.size() << '\n' ;
		//std::cout << "new_list_CB after update: ";
		/*for(int k = 0; k<new_list_CB.size();k++)
		{
			std::cout << new_list_CB[k] << ' ';
		}
		std::cout << '\n';*/

		for(int j = 0;j<new_list_CB.size();j++)
		{
			list_CB->push_back(new_list_CB[j]);
		}

		/*std::cout << "list_CB after copy: ";
		for(int k = 0; k<list_CB->size();k++)
		{
			std::cout << list_CB->at(k) << ' ';
		}
		std::cout << '\n';*/
		i += 1;
	}

	/*for(int itest = 0;itest<2;itest++)
	{	new_list_CB.clear();
		new_list_CB = update_GCH(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, X_r_x, X_r_y, list_bridged_receptors, list_CB, itest);
		list_CB->clear();	

		for(int j = 0;j<new_list_CB.size();j++)
		{
			list_CB->push_back(new_list_CB[j]);
		}
	}*/

	if(list_bridged_receptors->size() == 2 && new_list_CB.size() == 0)
	{
		new_list_CB.push_back(list_bridged_receptors->at(0));
		new_list_CB.push_back(list_bridged_receptors->at(1));
	}

	return new_list_CB;
}


void Make_reaction(MATRIX_FLOAT *a_d, MATRIX_FLOAT *a_on, MATRIX_FLOAT *a_off, std::vector<float> a_d_tot, std::vector<float> a_on_tot, std::vector<float> a_off_tot, MATRIX_INT L_L, MATRIX_INT L_R, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *T_l, std::vector<int> *T_r, std::vector<int> *L_to_R, std::vector<int> *list_bridged_receptors, std::vector<int> *list_CB, int r, long *idum, int step, std::ofstream *depleted)
{//choisit le ligand et le rÃ©cepteur et effectue la rÃ©action choisie auparavant 
    float a_sum;
    float theta, x, y;
    int nb;
    std::vector<int> new_list_CB;

	if (r == 1) //destruction SA
	{
		int alpha = Choose(a_d_tot, idum); //choix du ligand qui rÃ©agit

		int id = Choose((a_d->at(alpha-1)), idum); //choix du rÃ©cepteur qui rÃ©agit
		int i = L_L[alpha-1][id];

		*depleted << step << "," << i << "\n"; //print le récepteur qui est détruit, numéroté de 1 à N_R

		//MISE A JOUR
		(a_d->at(alpha-1))[0] -= (a_d->at(alpha-1))[id];
		(a_d->at(alpha-1))[id] = 0;

		for (int j = 1; j < L_R[i-1].size(); j++) 
		{
			theta = X_CM_theta->at((X_CM_theta->size()) - 1);
		    x = (X_l_x->at(L_R[i-1][j] - 1))*cos(theta) - (X_l_y->at(L_R[i-1][j] - 1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
			y = (X_l_x->at(L_R[i-1][j] - 1))*sin(theta) + (X_l_y->at(L_R[i-1][j] - 1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);

			if((pow(x - X_r_x->at(i-1),2) + pow(y - X_r_y->at(i-1),2) < D*D))
			{
				//clock_t startTime = clock();
				int id2 = Get_Id(L_L[(L_R[i-1][j]) - 1], i);
				//std::cout << "test L_R id - r1" << '\n';
				//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';
			    if(T_l->at(L_R[i-1][j] - 1) == 0)
				{
					(a_d->at(L_R[i-1][j] - 1))[0] -= (a_d->at(L_R[i-1][j] - 1))[id2];
				    (a_d->at(L_R[i-1][j] - 1))[id2] = 0.0;
				}
				if(T_l->at(L_R[i-1][j] - 1) == 1)
				{
					(a_on->at(L_R[i-1][j] - 1))[0] -= (a_on->at(L_R[i-1][j] - 1))[id2];
					(a_on->at(L_R[i-1][j] - 1))[id2] = 0.0;
				}
			}
		}
	    T_r->at(i-1) = 0;
	}


	else if(r == 2) //crÃ©ation pont
	{
		//std::cout << "reaction type 2" << '\n';
		int alpha = Choose(a_on_tot, idum); //choix du ligand qui rÃ©agit
		int id = Choose((a_on->at(alpha - 1)), idum); //choix du rÃ©cepteur qui rÃ©agit
		int i = L_L[alpha - 1][id];

		//MISE A JOUR
		(a_off->at(alpha - 1))[0] += k_off;
		(a_off->at(alpha - 1))[id] = k_off; 
	
		for (int k = 1; k < L_R[i-1].size(); k++)
		{
			theta = X_CM_theta->at((X_CM_theta->size()) - 1);
		    x = (X_l_x->at(L_R[i-1][k] - 1))*cos(theta) - (X_l_y->at(L_R[i-1][k] - 1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
			y = (X_l_x->at(L_R[i-1][k] - 1))*sin(theta) + (X_l_y->at(L_R[i-1][k] - 1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);

			//if(i == 124029 && alpha == 81)
			//	std::cout << "step " << step << " -- distance in reaction: " << pow(x - X_r_x->at(i-1),2) + pow(y - X_r_y->at(i-1),2) << " while max = " << D*D << "\n";
			if((pow(x - X_r_x->at(i-1),2) + pow(y - X_r_y->at(i-1),2) < D*D))
			{
				//clock_t startTime = clock();
				int id2 = Get_Id(L_L[(L_R[i-1][k]) - 1], i);
				//std::cout << "test L_R id - r2" << '\n';
				//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';			    
				if (T_l->at(L_R[i-1][k] - 1) == 1)
				{
					(a_on->at(L_R[i-1][k] - 1))[0] -= (a_on->at(L_R[i-1][k] - 1))[id2];
					(a_on->at(L_R[i-1][k] - 1))[id2] = 0.0;				
				}
				else if (T_l->at(L_R[i-1][k] - 1) == 0)
				{
					(a_d->at(L_R[i-1][k] - 1))[0] -= (a_d->at(L_R[i-1][k] - 1))[id2];
					(a_d->at(L_R[i-1][k] - 1))[id2] = 0.0;						
				}
			}
		}
		
		for (int k = 0; k < (a_on->at(alpha-1)).size(); k++)
		{
		    a_on->at(alpha-1).at(k) = 0.0;
		}
		
		T_l->at(alpha - 1) = 2;
		T_r->at(i - 1) = 2;
		L_to_R->at(alpha - 1) = i;	

		//update list_bridged_receptors et list_CB
		list_bridged_receptors->push_back(i);
		nb = list_bridged_receptors->size() - 1;
		new_list_CB.clear();

		new_list_CB = update_GCH(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, X_r_x, X_r_y, list_bridged_receptors, list_CB, nb);
		list_CB->clear();


		for(int k = 0;k<new_list_CB.size();k++)
		{
			list_CB->push_back(new_list_CB[k]);
		}

		/*std::cout << '\n' << "list_bridged_receptors after adding bridge: ";
		for(int k = 0; k<list_bridged_receptors->size();k++)
		{
			std::cout << list_bridged_receptors->at(k) << ' ';
		}
		std::cout << '\n';

		std::cout << "list_CB after adding bridge: ";
		for(int k = 0; k<list_CB->size();k++)
		{
			std::cout << list_CB->at(k) << ' ';
		}
		std::cout << '\n';*/
		//std::cout << "reaction: on lie ligand " << alpha << " et recepteur " << i << '\n';

		//std::cout << "end reaction 2" << '\n';
	}


	else if(r == 3) //destruction pont
	{
		//std::cout << "reaction type 3" << '\n';
		int alpha = Choose(a_off_tot, idum); //choix du ligand qui rÃ©agit

		int id = Choose((a_off->at(alpha-1)), idum); //choix de la position du rÃ©cepteur qui rÃ©agit
		int i = L_L[alpha-1][id]; //indice du rÃ©cepteur qui rÃ©agit


		//MISE A JOUR
		(a_off->at(alpha-1))[0] = 0.0;
		(a_off->at(alpha-1))[id] = 0.0; 

		T_l->at(alpha-1) = 1;
		T_r->at(i - 1) = 1;

		for (int m = 1; m < L_R[i-1].size(); m++)
		{

			theta = X_CM_theta->at((X_CM_theta->size()) - 1);
		    x = (X_l_x->at(L_R[i-1][m] - 1))*cos(theta) - (X_l_y->at(L_R[i-1][m] - 1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
			y = (X_l_x->at(L_R[i-1][m] - 1))*sin(theta) + (X_l_y->at(L_R[i-1][m] - 1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);

			if((pow(x - X_r_x->at(i-1),2) + pow(y - X_r_y->at(i-1),2) < D*D))
			{
				//clock_t startTime = clock();
				int id2 = Get_Id(L_L[(L_R[i-1][m]) - 1], i);
				//std::cout << "test L_R id - r3" << '\n';
				//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';		

				if (T_l->at(L_R[i-1][m] - 1) == 1 && L_R[i-1][m] != alpha)	
				{
					(a_on->at(L_R[i-1][m] - 1))[0] += k_on;
					(a_on->at(L_R[i-1][m] - 1))[id2] = k_on;				
				}
				else if (T_l->at(L_R[i-1][m] - 1) == 0)
				{
					(a_d->at(L_R[i-1][m] - 1))[0] += k_d;
					(a_d->at(L_R[i-1][m] - 1))[id2] = k_d;						
				}
			}
		}
		
		a_sum = 0;
		for (int k = 1; k < L_L[alpha-1].size(); k++)
		{
			theta = X_CM_theta->at((X_CM_theta->size()) - 1);
		    x = (X_l_x->at(alpha-1))*cos(theta) - (X_l_y->at(alpha-1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
			y = (X_l_x->at(alpha-1))*sin(theta) + (X_l_y->at(alpha-1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);
		
		    if ((T_r->at(L_L[alpha-1][k] - 1) == 1)  && (pow(x - X_r_x->at(L_L[alpha-1][k] - 1),2) + pow(y - X_r_y->at(L_L[alpha-1][k] - 1),2) < D*D))
		    {
		        a_on->at(alpha-1)[k] += k_on;
		        a_sum += k_on;
		    }
		    else
		    {
		        a_on->at(alpha-1)[k] = 0.0;
		    }
		}
		a_on->at(alpha-1)[0] = a_sum; 
		L_to_R->at(alpha - 1) = -1;


		//PRINTS
/*		std::cout << '\n' << '\n' << "NEW REACTION: " << '\n';
		std::cout << '\n' << "list_bridged_receptors before removing bridge: ";
		for(int k = 0; k<list_bridged_receptors->size();k++)
		{
			std::cout << list_bridged_receptors->at(k) << ' ';
		}
		std::cout << '\n';

		std::cout << "list_CB before removing bridge: ";
		for(int k = 0; k<list_CB->size();k++)
		{
			std::cout << list_CB->at(k) << ' ';
		}
		std::cout << '\n'<< "removed receptor: " << i << "\n" << "removed ligand: " << alpha << "\n";		
*/
		//update list_bridged_receptors et list_CB
		int idb = Get_Id(*list_bridged_receptors, i);
		//std::cout << "Get_Id in bridged receptors: " << idb << '\n';
		if(idb == -1 && list_bridged_receptors->at(0) == i)
			idb = 0;
		if(idb == list_bridged_receptors->size() && list_bridged_receptors->at(0) == i)
			idb = 0;
		if(i != list_bridged_receptors->at(idb))
			std::cout << "WARNING MKREACT" << "\n";
		list_bridged_receptors->erase(list_bridged_receptors->begin()+idb);


		int in_the_GCH = 0;
		int id_GCH = 0;
		for(int iGCH = 0; iGCH<list_CB->size(); iGCH++) //check si le pont était un CB
		{
			if(list_CB->at(iGCH) == i)
			{	//std::cout << "in the GCH" << "\n";
				in_the_GCH = 1;
				id_GCH = iGCH;
			}
		}

		/*for(int iGCH = 0; iGCH<list_CB->size(); iGCH++)
		{
			if(list_CB->at(iGCH) == nb - 1)
				list_CB->at(iGCH) = i;
		}*/

		//std::cout << "id_GCH: " << id_GCH << '\n';
		//std::cout << alpha << " " << i << '\n';
		if(in_the_GCH == 1)
		{
			if(list_bridged_receptors->size() < 2)
			{
				list_CB->erase(list_CB->begin()+id_GCH);
			}
			else
			{
				//std::cout << "bridge " << i << " in list_CB " << "\n";
				new_list_CB.clear();
				new_list_CB = gen_GCH(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, X_r_x, X_r_y, list_bridged_receptors, list_CB);
				list_CB->clear();

				//std::cout << '\n' << "new_list_CB generated after removing bridge: ";
				/*for(int k = 0; k<new_list_CB.size();k++)
				{
					std::cout << new_list_CB.at(k) << ' ';
				}
				std::cout << '\n';
*/

				for(int k = 0;k<new_list_CB.size();k++)
				{
					list_CB->push_back(new_list_CB[k]);
				}
			}
		}


		/*std::cout << '\n' << "list_bridged_receptors after removing bridge: ";
		for(int k = 0; k<list_bridged_receptors->size();k++)
		{
			std::cout << list_bridged_receptors->at(k) << ' ';
		}
		std::cout << '\n';

		std::cout << "list_CB after removing bridge: ";
		for(int k = 0; k<list_CB->size();k++)
		{
			std::cout << list_CB->at(k) << ' ';
		}
		std::cout << '\n';	*/
		//std::cout << "end reaction 3" << '\n';
	}    
	
	else
	{
	    std::cout << "error in choice of reaction: no reaction has been chosen" << '\n';
	}
}



void Update_Cell_L(std::vector<int> *ll, std::vector<int> *bl, std::vector<int> *hoc, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<int> *check, int N_L, int N_x, int N_y, float R_Cell, float X_CM_x_init, float X_CM_y_init, float X_CM_theta_init)
{//met à jour toutes les cellules après déplacement du virus
	float theta = X_CM_theta->at(X_CM_theta->size() - 1);
	float otheta = X_CM_theta_init;


	for (int i = 1; i < N_L + 1; i++)
	{

				float x = (X_l_x->at(i-1))*cos(theta) - (X_l_y->at(i-1))*sin(theta);
				float y = (X_l_x->at(i-1))*sin(theta) + (X_l_y->at(i-1))*cos(theta);
				    
				float x_p = x + X_CM_x->at(X_CM_x->size() - 1);
				float y_p = y + X_CM_y->at(X_CM_y->size() - 1);

				int ix = int((x_p + Taille_Systeme/2)/R_Cell) + 1;
				int iy = int((y_p + Taille_Systeme/2)/R_Cell) + 1;
				int nc = Cell_to_Index(N_x,N_y,ix,iy);



				float ox = (X_l_x->at(i-1))*cos(otheta) - (X_l_y->at(i-1))*sin(otheta);
				float oy = (X_l_x->at(i-1))*sin(otheta) + (X_l_y->at(i-1))*cos(otheta);
				    
				float ox_p = ox + X_CM_x_init;
				float oy_p = oy + X_CM_y_init;

				int oix = int((ox_p + Taille_Systeme/2)/R_Cell) + 1;
				int oiy = int((oy_p + Taille_Systeme/2)/R_Cell) + 1;
				int oc = Cell_to_Index(N_x,N_y,oix,oiy);

				if(oc != nc)
				{
					//std::cout << "test a" << "\n";
					UpGrade(oc, nc, i, ll, bl, hoc);
					//std::cout << "test b" << '\n';
				}
			
		
	}

}


void Update_Cell_R(std::vector<int> *ll, std::vector<int> *bl, std::vector<int> *hoc, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<float> *x_init_R, std::vector<float> *y_init_R, int N_x, int N_y, float R_Cell)
{//met à jour toutes les cellules après déplacement du virus

	for (int i = 1; i < N_R + 1; i++)
	{

				float x = (X_r_x->at(i-1));
				float y = (X_r_y->at(i-1));
				    
				float x_p = x;
				float y_p = y;

				int ix = int((x_p + Taille_Systeme/2)/R_Cell) + 1;
				int iy = int((y_p + Taille_Systeme/2)/R_Cell) + 1;
				int nc = Cell_to_Index(N_x,N_y,ix,iy);


				float ox_p = x_init_R->at(i-1);
				float oy_p = y_init_R->at(i-1);

				int oix = int((ox_p + Taille_Systeme/2)/R_Cell) + 1;
				int oiy = int((oy_p + Taille_Systeme/2)/R_Cell) + 1;
				int oc = Cell_to_Index(N_x,N_y,oix,oiy);

				if(oc != nc)
				{
					//std::cout << "test a" << "\n";
					UpGrade(oc, nc, i, ll, bl, hoc);
					//std::cout << "test b" << '\n';
				}	
	}

}



void Update_Lists(MATRIX_INT *L_L, MATRIX_INT *L_R, MATRIX_FLOAT *a_d, MATRIX_FLOAT *a_on, MATRIX_FLOAT *a_off, std::vector<float> *a_d_tot, std::vector<float> *a_on_tot, std::vector<float> *a_off_tot, std::vector<float> *affinity, MATRIX_INT Neigh, std::vector<int> hoc_R, std::vector<int> hoc_L, std::vector<int> ll_L, std::vector<int> ll_R, std::vector<float> X_l_x, std::vector<float> X_l_y, std::vector<float> X_r_x, std::vector<float> X_r_y, std::vector<float> X_CM_x, std::vector<float> X_CM_y, std::vector<float> X_CM_theta, int N_x, int N_y, float R_Cell, int N_L, std::vector<int> T_l, std::vector<int> T_r, std::vector<int> *L_to_R, int step)
{//update L_L, L_R, a_d, a_on, a_off en les vidant, puis en les remplissant Ã  nouveau avec les fonctions Initialize
    for(int i=0; i<N_L; i++)
    {
        (L_L->at(i)).clear();
        (a_d->at(i)).clear();
        (a_on->at(i)).clear();
        (a_off->at(i)).clear();
    }
    //std::cout << "A" << '\n';
    for(int i=0; i<N_R; i++)
    {
        (L_R->at(i)).clear();
    }
    //std::cout << "B" << '\n';
    a_d_tot->clear();
    a_on_tot->clear();
    a_off_tot->clear();
    affinity->clear();
    //std::cout << "C" << '\n';   
	Initialize_L_L(L_L, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_R, ll_R, T_l, *L_to_R, N_x, N_y, R_Cell, N_L, step);
	//std::cout << "C bis" << '\n';
	Initialize_L_R(L_R, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_L, ll_L, T_r, T_l, *L_to_R, N_x, N_y, R_Cell, N_L, step);

    //std::cout << "D" << '\n';
	Initialize_Affinity(a_d, L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, L_to_R, T_l, T_r, k_d, N_L, 0, 1);
	Initialize_Affinity(a_on, L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, L_to_R, T_l, T_r, k_on, N_L, 1, 1);
    Initialize_Affinity(a_off, L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, L_to_R, T_l, T_r, k_off, N_L, 2, 2);
    //std::cout << "E" << '\n';
    *a_d_tot = Get_a_X_tot(a_d, N_L);
	*a_on_tot = Get_a_X_tot(a_on, N_L);
	*a_off_tot = Get_a_X_tot(a_off, N_L);
    //std::cout << "F" << '\n';
	*affinity = Get_a_tot(a_d_tot, a_on_tot, a_off_tot, N_L);
}


void Update_L_L(MATRIX_INT *L_L, MATRIX_INT *L_R, MATRIX_INT *L_L_old, int N_L, std::vector<int> *check, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> X_r_x, std::vector<float> X_r_y, std::vector<float> X_CM_x, std::vector<float> X_CM_y, std::vector<float> X_CM_theta, MATRIX_INT Neigh, std::vector<int> hoc_R, std::vector<int> ll_R, std::vector<int> T_l, std::vector<int> L_to_R, int N_x, int N_y, float R_Cell)
{
	for(int k = 0; k < N_L; k++)
	{	
		(L_L_old->at(k)).clear();
		L_L_old->at(k) = L_L->at(k);
	}

	float theta;
    float x, y;
    int alpha;

	//std::cout << "test 1bis" << '\n';

	for (int beta = 0; beta < check->size(); ++beta)
	{ 
		//std::cout << check->size() << '\n';
		//std::cout << beta << '\n';
		alpha = check->at(beta) - 1;
		(L_L->at(alpha)).clear();
		theta = X_CM_theta.at((X_CM_theta.size()) - 1);
	    x = (X_l_x->at(alpha))*cos(theta) - (X_l_y->at(alpha))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
		y = (X_l_x->at(alpha))*sin(theta) + (X_l_y->at(alpha))*cos(theta) + X_CM_y[X_CM_y.size() - 1];

		int ix = int((x + Taille_Systeme/2)/R_Cell) + 1;
		int iy = int((y + Taille_Systeme/2)/R_Cell) + 1;
		int index = Cell_to_Index(N_x, N_y, ix, iy);
        
		for (int j = 0; j < 9; ++j)
		{
			int i = Neigh[index-1][j];
			
			if(hoc_R[i] != 0)
			{
	            if(pow(x - X_r_x[hoc_R[i] - 1],2) + pow(y - X_r_y[hoc_R[i] - 1],2) < R_Verlet*R_Verlet)
    				(L_L->at(alpha)).push_back(hoc_R[i]); 
				int g = ll_R[hoc_R[i]];
				
				while(g != 0)
				{
	                if(pow(x - X_r_x[g - 1],2) + pow(y - X_r_y[g - 1],2) < R_Verlet*R_Verlet)
					    (L_L->at(alpha)).push_back(g);
					int k = g;	
					g = ll_R[k];
				}
			}
		}
		if(T_l[alpha] == 2 && pow(x - X_r_x[L_to_R[alpha] - 1],2) + pow(y - X_r_y[L_to_R[alpha] - 1],2) > R_Verlet*R_Verlet) //gère l'exception due à l'imprécision dans check_breakage
			(L_L->at(alpha)).push_back(L_to_R[alpha]);

		int N = (L_L->at(alpha)).size();
		(L_L->at(alpha)).insert((L_L->at(alpha)).begin(), N);
	}
}

bool check_intersection(int element, std::vector<int> v)
{
	for(int j = 0; j < v.size(); j++)
	{
		if(element == v[j])
		{
			return true;
		}		
	}
	return false;
}


void Update_L_R(MATRIX_INT *L_L_old, MATRIX_INT *L_L, MATRIX_INT *L_R, std::vector<int> *check)
{
	clock_t startTime = clock();
	int alpha;
	bool check_inter;
	std::vector<int> v_old;
	std::vector<int> v_new;

	for (int beta = 0; beta < check->size(); ++beta)
	{ 
		alpha = check->at(beta);

		v_old = L_L_old->at(alpha - 1);
		v_new = L_L->at(alpha - 1);

		v_old.erase(v_old.begin());
		v_new.erase(v_new.begin());

	    std::sort(v_old.begin(), v_old.end());
	    std::sort(v_new.begin(), v_new.end());
	 
	    std::vector<int> v_intersection;
	 
	    std::set_intersection(v_old.begin(), v_old.end(),
	                          v_new.begin(), v_new.end(),
	                          std::back_inserter(v_intersection));


	    //ENLEVER LIGAND POUR LES RECEPTEURS SORTANTS
	    for(int j = 0; j < v_old.size(); j++)
	    {
	    	if(v_intersection.size() != 0)
	    	{
	    		check_inter = check_intersection(v_old[j],v_intersection); 

	    		if(check_inter == false)
	    		{
					    int i = v_old[j];
	    				int id = Get_Id(L_R->at(i-1), alpha);    				

	    				(L_R->at(i-1)).erase(L_R->at(i-1).begin() + id);
	    				(L_R->at(i-1).at(0)) -= 1;    				
	    		}
	    	}
	    	else
	    	{
		    	int i = v_old[j];
				int id = Get_Id(L_R->at(i-1), alpha);

				(L_R->at(i-1)).erase(L_R->at(i-1).begin() + id);	
				(L_R->at(i-1).at(0)) -= 1;			
	    	}
	    }


	    //AJOUTER LIGAND POUR LES RECEPTEURS ENTRANTS
	    for(int j = 0; j < v_new.size(); j++)
	    {
	    	if(v_intersection.size() != 0)
	    	{
	    		check_inter = check_intersection(v_new[j],v_intersection);
	    		if(check_inter == false)
	    		{
	    				int i = v_new[j];
	    				(L_R->at(i-1)).push_back(alpha);
	    				L_R->at(i-1).at(0) += 1;
	    		}
	    	}
	    	else
	    	{
		    	int i = v_new[j];
				(L_R->at(i-1)).push_back(alpha);	
				L_R->at(i-1).at(0) += 1;
	    	}
	    }
	}

	v_old.clear();
	v_new.clear();

}


void Update_Affinity(MATRIX_FLOAT *a_on, MATRIX_FLOAT *a_off, MATRIX_FLOAT *a_d, std::vector<float> *a_d_tot, std::vector<float> *a_on_tot, std::vector<float> *a_off_tot, std::vector<float> *affinity, MATRIX_INT *L_L, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *L_to_R, std::vector<int> T_l, std::vector<int> T_r, float N_L)
{
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

	Initialize_Affinity(a_d, L_L, X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x, X_r_y, L_to_R, T_l, T_r, k_d, N_L, 0, 1);
	Initialize_Affinity(a_on, L_L, X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x, X_r_y, L_to_R, T_l, T_r, k_on, N_L, 1, 1);
    Initialize_Affinity(a_off, L_L, X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x, X_r_y, L_to_R, T_l, T_r, k_off, N_L, 2, 2);

    *a_d_tot = Get_a_X_tot(a_d, N_L);
	*a_on_tot = Get_a_X_tot(a_on, N_L);
	*a_off_tot = Get_a_X_tot(a_off, N_L);

	*affinity = Get_a_tot(a_d_tot, a_on_tot, a_off_tot, N_L);

}



std::vector<float> intersection_r1_r2(std::vector<float> v,std::vector<float> a,std::vector<float> w,std::vector<float> b)
{   
	//std::cout << "E" << '\n'; 
	float den, t;
	std::vector<float> inter;

	inter.clear();

	//std::cout << "E1" << '\n'; 	
	den=b[0]*a[1]-b[1]*a[0];
	//std::cout << "den: " << den << '\n';

    if(abs(den) < pow(10.,-10))
    {	 
      	//std::cout << "E2" << '\n'; 
        //std::cout << "overlapping CBs" << '\n';
        inter.push_back((v[0]+w[0])/2);
        inter.push_back((v[1]+w[1])/2);

        return inter;
    }
   	//std::cout << "E3" << '\n'; 

    t=b[0]*(w[1]-v[1])-b[1]*(w[0]-v[0]);
    //std::cout << "t: " << t << '\n';
    t=t/den;

	//std::cout << "E4" << '\n'; 

	//std::cout << "inter: " << v[0] + a[0]*t << " " << v[1] + a[1]*t << '\n';
    inter.push_back(v[0] + a[0]*t);
    inter.push_back(v[1] + a[1]*t);
	//std::cout << "E5" << '\n';

    return inter;
}

MATRIX_FLOAT Limit_CM_adv(MATRIX_FLOAT X_GCH, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<int> *list_CB, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<int> *L_to_R)
{
	//std::cout << "F" << '\n';
	MATRIX_FLOAT r;
	r.resize(2, std::vector<float>(0));

	int ncb;
	ncb = X_GCH.size();
	//std::cout << "X_GCH.size(): " << X_GCH.size() << '\n'; 

	float dd;
	std::vector<float> r1, r2, r1p, r2p, z1, z2, zm, u, k1, k2, k3, k4, wx, wy;

	std::vector<float> X_ligand;

	r1.clear();
	r2.clear();
	r1p.clear();
	r2p.clear();
	z1.clear();
	z2.clear();
	zm.clear();
	u.clear();
	k1.clear();
	k2.clear();
	k3.clear();
	k4.clear();
	wx.clear();
	wy.clear();

	if(ncb == 2)
	{
		//std::cout << "Test A1" << "\n";

		//std::cout << "list_CB: " << list_CB->at(0) << " " << list_CB->at(1) << '\n';
		for(int j = 0; j<2 ; j++)
		{	
			X_ligand.clear();
			X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_CB->at(0),0);

			r1p.push_back(X_GCH[0][j] - X_ligand[j]);
			
			X_ligand.clear();
			X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_CB->at(1),0);

			r2p.push_back(X_GCH[1][j] - X_ligand[j]);
		}	
		//std::cout << "position: on regarde recepteur " << list_CB->at(0) << '\n';	
		//std::cout << "position: on regarde recepteur " << list_CB->at(1) << '\n';
		

        //std::cout << "r1p: " << r1p[0] << " " << r1p[1] << '\n';
        //std::cout << "r2p: " << r2p[0] << " " << r2p[1] << '\n';
        r1 = r1p;
        r2 = r2p;
        z1 = center(r1, r2, D);
        z2 = center(r2, r1, D);

        //std::cout << "z1: " << z1[0] << " " << z1[1] << '\n';
        //std::cout << "z2: " << z2[0] << " " << z2[1] << '\n';

        zm.push_back((z1[0]+z2[0])/2);
        zm.push_back((z1[1]+z2[1])/2);
        dd = sqrt(pow(zm[0]-r1[0],2)+pow(zm[1]-r1[1],2));
        dd = D - dd;

        if(dd<0)
            std::cout << "problems in Limit_CM_adv, ncb==2" << '\n';
        
        u.push_back(r1[0] - zm[0]);
        u.push_back(r1[1] - zm[1]);

        u[0]=u[0]/(sqrt(pow(r1[0]-zm[0],2)+pow(r1[1]-zm[1],2)));
        u[1]=u[1]/(sqrt(pow(r1[0]-zm[0],2)+pow(r1[1]-zm[1],2)));

        u[0] = u[0]*dd;
        u[1] = u[1]*dd;

        //std::cout << "u: " << u[0] << " " << u[1] << '\n';

        k1.push_back(z1[0]+u[0]);
        k1.push_back(z1[1]+u[1]);

        k2.push_back(z1[0]-u[0]);
        k2.push_back(z1[1]-u[1]);

        k3.push_back(z2[0]+u[0]);
        k3.push_back(z2[1]+u[1]);
 
        k4.push_back(z2[0]-u[0]);
        k4.push_back(z2[1]-u[1]);

        wx.push_back(k1[0]);
        wx.push_back(k2[0]);
        wx.push_back(k3[0]);
        wx.push_back(k4[0]);

        wy.push_back(k1[1]);
        wy.push_back(k2[1]);
        wy.push_back(k3[1]);
        wy.push_back(k4[1]);
	}

	else if(ncb>2)
	{
		//std::cout << "Test A2" << "\n";
		std::vector<float> r3, r3p, v, w, a, b, vint;
		r3.clear();
		r3p.clear();
		v.clear();
		w.clear();
		a.clear();
		b.clear();

		for(int i=0; i<ncb;i++)
		{
			r1.clear();
			r1p.clear();
			r2.clear();
			r2p.clear();
			r3.clear();
			r3p.clear();
			v.clear();
			w.clear();
			a.clear();
			b.clear();
			vint.clear();

			//std::cout << "for i = " << i << "\n";

			for(int j = 0; j<X_GCH[i].size();j++)
			{
				X_ligand.clear();
				X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_CB->at(i),0);
				//std::cout << "A2 1: " << i << " " << X_GCH[i][j] - X_ligand[j] << "\n";
				r1p.push_back(X_GCH[i][j] - X_ligand[j]);


				X_ligand.clear();
				X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_CB->at((i+1) % ncb),0);
				//std::cout << "A2 2: " << (i+1) % ncb << " " << (X_GCH[(i+1) % ncb][j]) - X_ligand[j] << "\n";
				r2p.push_back((X_GCH[(i+1) % ncb][j]) - X_ligand[j]);

				X_ligand.clear();
				X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_CB->at((i+2) % ncb),0);
				//std::cout << "A2 3: " << (i+2) % ncb << " " << (X_GCH[(i+2) % ncb][j]) - X_ligand[j] << "\n";
				r3p.push_back((X_GCH[(i+2) % ncb][j]) - X_ligand[j]);
			}

			/*std::cout << "position recepteur r1p: " << X_GCH[i][0] << " " << X_GCH[i][1] << "\n";
			std::cout << "position recepteur r2p: " << X_GCH[(i+1) % ncb][0] << " " << X_GCH[(i+1) % ncb][1] << "\n";
			std::cout << "position recepteur r3p: " << X_GCH[(i+2) % ncb][0] << " " << X_GCH[(i+2) % ncb][1] << "\n";
			*/

			/*std::cout << '\n' << "r1p: " << r1p[0] << " " << r1p[1] << "\n";
			std::cout << "r2p: " << r2p[0] << " " << r2p[1] << "\n";
			std::cout << "r3p: " << r3p[0] << " " << r3p[1] << "\n";*/

			r1 = r1p;
			r2 = r2p;
			r3 = r3p;

			v = center(r1, r2, D);
			w = center(r2, r3, D);

			a.push_back(r1[0]-r2[0]);
			a.push_back(r1[1]-r2[1]);
			b.push_back(r2[0]-r3[0]);
			b.push_back(r2[1]-r3[1]);		

			vint = intersection_r1_r2(v,a,w,b);
			//std::cout << "vint: " << vint[0] << " " << vint[1] << '\n'; 	

			wx.push_back(vint[0]);
			//std::cout << "post E2" << '\n'; 	
			wy.push_back(vint[1]);
			//std::cout << "post E3" << '\n'; 	
		}
		//std::cout << "post E4" << '\n'; 

		/*if(list_CB->at(0) == 142371 && list_CB->size() == 4)
		{
			std::cout << "\n" << "wxy: ";
			for(int k = 0; k<wx.size(); k++)
			{
				std::cout << "[" << wx[k] << "," << wy[k] << "]";
			}
		}	//std::cout << "\n" << "\n";*/
	}

	else if(ncb == 1) //gère le cas où on a un seul pont
	{
		X_ligand.clear();
		X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_CB->at(0),0);

		r[0].push_back(X_GCH[0][0] - X_ligand[0] - D);
		r[0].push_back(X_GCH[0][1] - X_ligand[1] - D);

		r[1].push_back(X_GCH[0][0] - X_ligand[0] + D);
		r[1].push_back(X_GCH[0][1] - X_ligand[1] + D);

		return r;
	}


	else
	{
		std::cout << "Moins de 1 CB dans Limit_CM_adv" << '\n';
	}

	//std::cout << "post ifs 1" << '\n'; 	

	//r[0].push_back(X_CM_x->at(X_CM_x->size() -1));
	//r[0].push_back(X_CM_y->at(X_CM_y->size() -1));
	//r[0].push_back(X_CM_theta->at(X_CM_theta->size() -1));
	r[0].push_back(0.0);
	r[0].push_back(0.0);

	//std::cout << "post ifs 2" << '\n'; 
	
	//r[1].push_back(X_CM_x->at(X_CM_x->size() -1));
	//r[1].push_back(X_CM_y->at(X_CM_y->size() -1));
	//r[1].push_back(X_CM_theta->at(X_CM_theta->size() -1));

	r[1].push_back(0.0);
	r[1].push_back(0.0);

	//std::cout << "Test B: " << wx.size() << "\n";
	for(int irc=0; irc < wx.size(); irc++)
	{

		//std::cout << "1: " << wx[irc] << '\n';
		if(r[0][0]>wx[irc])
		{	
			//std::cout << "test1" <<"\n";
           	r[0][0]=wx[irc];
        }
  		//std::cout << "2: " << wy[irc] << '\n';
        if(r[0][1]>wy[irc])
        {	
        	//std::cout << "test2" <<"\n";
            r[0][1]=wy[irc];
        }
		//std::cout << "3: " << wx[irc] << '\n';
        if(r[1][0]<wx[irc])
        {
			//std::cout << "test3" <<"\n";
			r[1][0]=wx[irc];
        }
       // std::cout << "4: " << wy[irc] << '\n';
        if(r[1][1]<wy[irc])
		{
			//std::cout << "test4" <<"\n";
            r[1][1]=wy[irc];	
		}			
	}

	//std::cout << "Test C" << "\n";

	return r;
}

std::vector<float> UnSamp(std::vector<float> r_min, std::vector<float> r_max, long *idum)
{ //return a vector inside the rectangle identified by r_min, r_max
	//std::cout << "G" << '\n';  	
  	std::vector<float> rnd, r;
  	r.clear();
  	rnd.clear();

	rnd.push_back(ran2(idum));
	rnd.push_back(ran2(idum));
	//std::cout << "random numbers: " << rnd[0] << " " << rnd[1] << '\n';
    r.push_back(r_min[0]+(r_max[0]-r_min[0])*rnd[0]);
    //std::cout << r_min[0] << " " << r_max[0]-r_min[0] << '\n'; 
    r.push_back(r_min[1]+(r_max[1]-r_min[1])*rnd[1]);
	//std::cout << r_min[1] << " " << r_max[1]-r_min[1] << '\n'; 
    
    return r;
}


int check_breakage_CB(std::vector<std::vector<float>> X_GCH, std::vector<float> r_trial, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y,std::vector<float> *X_CM_theta, std::vector<float> *X_l_x,std::vector<float> *X_l_y, std::vector<int> *L_to_R, std::vector<int> *list_CB, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *list_bridged_receptors)
{
	//std::cout << "H" << '\n';
	int nbk = 1;
	float dist = 0;
	//std::cout << "check_break_CB begin: " << '\n';
	std::vector<float> X_ligand;
	std::vector<float> new_X_ligand;
	std::vector<float> X_receptor;
	std::vector<float> center;
	center.push_back(0.0);
	center.push_back(0.0);

	for(int irc = 0; irc<X_GCH.size(); irc++)
	{
		//std::cout << "irc: " << irc << '\n';
		X_ligand.clear();
		X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_CB->at(irc),0);
		//std::cout << "H1" << '\n';
		new_X_ligand.clear();
		new_X_ligand = get_new_position_linked_ligand(X_l_x, X_l_y, r_trial, X_CM_theta, L_to_R, list_CB->at(irc));

		X_GCH[irc][0] -= new_X_ligand[0];
		X_GCH[irc][1] -= new_X_ligand[1];

		//std::cout << "H2" << '\n';
		dist = distance2(X_GCH[irc],center);
		//std::cout << "H3" << '\n';
		X_GCH[irc][0] += new_X_ligand[0];
		X_GCH[irc][1] += new_X_ligand[1]; 
		//std::cout << "H4" << '\n';
	/*for(int irc = 0; irc < list_bridged_receptors->size();irc++)
	{	
		X_receptor.clear();
		X_ligand.clear();
		X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_bridged_receptors->at(irc));
		X_receptor.push_back(X_r_x->at(list_bridged_receptors->at(irc)-1) - X_ligand[0]);
		X_receptor.push_back(X_r_y->at(list_bridged_receptors->at(irc)-1) - X_ligand[1]);
		dist = distance2(X_receptor,X_ligand);*/
		//std::cout << dist << " " << D*D << '\n' << '\n';
		if(dist > D*D)
		{
			
			nbk = 0;
			return nbk;
		}
		
	}
	return nbk;
}

std::vector<float> Sample_Conf_CGH(std::vector<float> *X_r_x, std::vector<float> *X_r_y,std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<int> *list_bridged_receptors, std::vector<int> *list_CB, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<int> *L_to_R, long *idum)
{
	//std::cout << "I" << '\n';
	float theta;
	theta = X_CM_theta->at(X_CM_theta->size()-1);
	int ir = 0;
	int ngch;
	ngch = list_CB->size();

	std::vector<float> wp;
	wp.push_back(X_r_x->at(list_CB->at(ir) - 1));
	wp.push_back(X_r_y->at(list_CB->at(ir) - 1));

	MATRIX_FLOAT X_GCH;
	X_GCH.resize(ngch, std::vector<float>(0));
	X_GCH[0].push_back(wp[0]);
	X_GCH[0].push_back(wp[1]);

	std::vector<float> r_trial, r_trial_temp;

	ir += 1;

	while(ir<ngch)
	{
		wp.clear();
		wp.push_back(X_r_x->at(list_CB->at(ir) - 1));
		wp.push_back(X_r_y->at(list_CB->at(ir) - 1));

		r_trial = wp;

		X_GCH[ir].push_back(r_trial[0]);
		X_GCH[ir].push_back(r_trial[1]);

		ir += 1;
	}

	MATRIX_FLOAT r;

	r = Limit_CM_adv(X_GCH,X_CM_x,X_CM_y,X_CM_theta,list_CB,X_l_x,X_l_y,L_to_R); //r[0] = r_min, r[1] = r_max
	
	
	int no_bk = 0;
	int test = 0;

	float x,y,temp_theta;
	while(no_bk < 0.5)
	{
		//std::cout << "Test while" << "\n";
		r_trial.clear();
		r_trial_temp.clear();
		//std::cout << "I 1" << '\n';
		//std::cout << "list_CB size: " << list_CB->size() << '\n';
		//std::cout << '\n' << "idum: " << *idum << '\n';
		r_trial_temp = UnSamp(r[0],r[1], idum);

		//On se remet dans le repère du CM
		r_trial.push_back(r_trial_temp[0]*cos(theta) - r_trial_temp[1]*sin(theta) + X_CM_x->at(X_CM_x->size()-1));
		r_trial.push_back(r_trial_temp[0]*sin(theta) + r_trial_temp[1]*cos(theta) + X_CM_y->at(X_CM_y->size()-1));  
		//std::cout << "all r_trial :" << r_trial[0] << " " << r_trial[1] << '\n' << "rectangle: " << r[0][0] << ' ' << r[0][1] << '\n' << r[1][0] << " " << r[1][1] << '\n';
		//std::cout << "r_trial :" << r_trial[0] << " " << r_trial[1] << '\n' << "rectangle: " << r[0][0] << ' ' << r[0][1] << '\n' << r[1][0] << " " << r[1][1] << '\n';

		//std::cout << "r_trial: " << r_trial[0] << " " << r_trial[1] << "\n";
		//std::cout << "I 2" << '\n';
		no_bk = check_breakage_CB(X_GCH,r_trial,X_CM_x, X_CM_y,X_CM_theta,X_l_x,X_l_y,L_to_R, list_CB, X_r_x, X_r_y, list_bridged_receptors);
		//std::cout << "I 3" << '\n';		
		test += 1;

		if(test % 100000 == 0)
		{
			std::cout << "while loop for " << test << " times " << '\n';
			//std::cout << "L_to_R[81]: " << L_to_R->at(80) << " =? " << list_CB->at(2) << "\n";
			//std::cout << "distance: "
			/*temp_theta = X_CM_theta->at((X_CM_theta->size()) - 1);
			x = (X_l_x->at(80))*cos(temp_theta) - (X_l_y->at(80))*sin(temp_theta) + X_CM_x->at(X_CM_x->size() - 1);
			y = (X_l_x->at(80))*sin(temp_theta) + (X_l_y->at(80))*cos(temp_theta) + X_CM_y->at(X_CM_y->size() - 1);
			std::cout << " during Sample_Conf_CGH -> distance entre ligand 81 et receptor 124029: " << pow(x - X_r_x->at(124029-1),2) + pow(y - X_r_y->at(124029-1),2) << " while max = " << D*D << "\n";
*/

			std::cout << "x_CM = [" << X_CM_x->at(X_CM_x->size()-1) << "," << X_CM_y->at(X_CM_y->size()-1) << "]" << '\n';
			std::cout << "r_trial :" << r_trial[0] << " " << r_trial[1] << '\n' << "rectangle: " << r[0][0] << ' ' << r[0][1] << '\n' << r[1][0] << " " << r[1][1] << '\n';
			std::cout << "list_CB size: " << list_CB->size() << '\n' << '\n';
	
		}
		//if(test == 1)
			//std::cout << "r_trial :" << r_trial[0] << " " << r_trial[1] << '\n' << "rectangle: " << r[0][0] << ' ' << r[0][1] << '\n' << r[1][0] << " " << r[1][1] << '\n';

		//std::cout << "I 4" << '\n';
	}
	//std::cout << "I 5" << '\n';

	//temporary theta (no move)
	//r_trial.push_back(0.0);
	//std::cout << "I 6" << '\n';
	return r_trial;
}


std::vector<float> theta_minmax(std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<int> *L_to_R, std::vector<int> *list_bridged_receptors, int step, std::ofstream *thetaminmax)
{
	float theta_min,theta_max, theta,alpha;
	float new_theta_min,new_theta_max,true_new_theta_min,true_new_theta_max;
	float xrl2, xr2, xl2, xr, xl;
	float dtest1, dtest2;
	std::vector<float> X_receptor, X_ligand, cm, Lim_theta;
	std::vector<float> X_l_max, X_l_min;

	theta_max = 0;
	theta_min = 0;
	theta = X_CM_theta->at(X_CM_theta->size()-1);
	cm.push_back(X_CM_x->at(X_CM_x->size()-1));
	cm.push_back(X_CM_y->at(X_CM_y->size()-1));
	//std::cout << X_CM_x->at(X_CM_x->size()-1) << ", " << X_CM_y->at(X_CM_y->size()-1) << ", " << X_CM_theta->at(X_CM_theta->size()-1) << '\n';

	for(int i = 0; i<list_bridged_receptors->size(); i++)
	{
		X_receptor.clear();
		//X_receptor.push_back(X_r_x->at(list_bridged_receptors->at(i) - 1)*cos(theta) + X_r_y->at(list_bridged_receptors->at(i) - 1)*sin(theta) - X_CM_x->at(X_CM_x->size()-1));
		//X_receptor.push_back(X_r_y->at(list_bridged_receptors->at(i) - 1)*cos(theta) - X_r_x->at(list_bridged_receptors->at(i) - 1)*sin(theta) - X_CM_y->at(X_CM_y->size()-1));

		X_receptor.push_back(X_r_x->at(list_bridged_receptors->at(i) - 1));
		X_receptor.push_back(X_r_y->at(list_bridged_receptors->at(i) - 1));

		X_ligand.clear();
		//X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_bridged_receptors->at(i),1);
		X_ligand = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_bridged_receptors->at(i),0);
		//std::cout << X_ligand[0] << ", " << X_ligand[1] << "\n";


		xrl2 = distance2(X_receptor,X_ligand);
		//std::cout << "initial distance for ligand " << i << ": " << xrl2 << "\n";
		xr2 = distance2(X_receptor,cm);
		xl2 = distance2(X_ligand,cm);
		xr = sqrt(xr2);
		xl = sqrt(xl2);

		alpha = acos((xr2 + xl2 - xrl2) / (2*xr*xl));
		//if(step == 2685)
			//std::cout << "ligand " << i << ": " << alpha << "\n";
		/*if(i == 31)
		{
			std::cout << (xr2 + xl2 - (D*D)) / (2*xr*xl) << "\n";
			std::cout << xr << " " << xl << "\n";
		}*/

		new_theta_min = acos((xr2 + xl2 - (D*D)) / (2*xr*xl)) - alpha; 
		new_theta_max = - new_theta_min - 2*alpha;
		//if(step == 2685)
			//std::cout << new_theta_min << "," << new_theta_max << "\n";
		//std::cout << acos((xr2 + xl2 - (D*D)) / (2*xr*xl)) << "," ;

		X_l_max.clear();
		X_l_min.clear();
		X_l_max = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_bridged_receptors->at(i),new_theta_max);
		X_l_min = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_bridged_receptors->at(i),new_theta_min);
		//std::cout << X_l_max[0] << ", " << X_l_max[1] << "\n";
		//std::cout << X_l_min[0] << ", " << X_l_min[1] << "\n";
		//X_l_max.push_back(X_ligand[0]*cos(new_theta_max) - X_ligand[1]*sin(new_theta_max));
		//X_l_max.push_back(X_ligand[0]*sin(new_theta_max) + X_ligand[1]*cos(new_theta_max));
		//X_l_min.push_back(X_ligand[0]*cos(new_theta_min) - X_ligand[1]*sin(new_theta_min));
		//X_l_min.push_back(X_ligand[0]*sin(new_theta_min) + X_ligand[1]*cos(new_theta_min));

		dtest1 = distance2(X_receptor,X_l_max);
		dtest2 = distance2(X_receptor,X_l_min);
		//if(step == 2685)
			//std::cout << i << " -> distances: " << dtest1 << ", " << dtest2 << ", " << D*D << "\n";		
		//REGLE LES PROBLEMES D'ARRONDIS DANS LE ARCCOS
		if(dtest1 > D*D + 0.2 || dtest2 > D*D + 0.2 || dtest1 < D*D - 0.2 || dtest2 < D*D - 0.2)
		{
			//if(step == 2685)
				//std::cout << i << " in if: " << dtest1 << ", " << dtest2 << ", " << D*D << "\n";
			true_new_theta_min = - new_theta_max;
			true_new_theta_max = - new_theta_min;
		}	
		else
		{
			/*if(dtest1 > D*D)
				true_new_theta_min = new_theta_min;
			if(dtest2 > D*D)
				true_new_theta_max = new_theta_max;
			if(dtest1 <= D*D && dtest2 <= D*D)
			{*/
			true_new_theta_min = new_theta_min;
			true_new_theta_max = new_theta_max;
			//}

		}

		X_l_max.clear();
		X_l_min.clear();
		X_l_max = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_bridged_receptors->at(i), true_new_theta_max);
		X_l_min = get_position_linked_ligand(X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_bridged_receptors->at(i), true_new_theta_min);
		//X_l_max.push_back(X_ligand[0]*cos(new_theta_max) - X_ligand[1]*sin(new_theta_max));
		//X_l_max.push_back(X_ligand[0]*sin(new_theta_max) + X_ligand[1]*cos(new_theta_max));
		//X_l_min.push_back(X_ligand[0]*cos(new_theta_min) - X_ligand[1]*sin(new_theta_min));
		//X_l_min.push_back(X_ligand[0]*sin(new_theta_min) + X_ligand[1]*cos(new_theta_min));

		dtest1 = distance2(X_receptor,X_l_max);
		dtest2 = distance2(X_receptor,X_l_min);
		//if(step == 2685)
			//std::cout << i << ": " << dtest1 << ", " << dtest2 << ", " << D*D << "\n";
		//std::cout << "initial distance for ligand " << i << ": " << dtest1 << ", " << dtest2 << ", compared to " << D*D << "\n";
		//if(step == 2685)
			//std::cout << true_new_theta_min << "," << true_new_theta_max << "\n";

		if(i == 0)
		{
			if(alpha != alpha || true_new_theta_min != true_new_theta_min || true_new_theta_max != true_new_theta_max)
			{
				theta_min = 1.0;
				theta_max = -1.0;
			}
			else
			{
				theta_max = true_new_theta_max;
				theta_min = true_new_theta_min;
			}
			//if(true_new_theta_max > 0)
				//std::cout << "WARNING WARNING WARNING ligand " << i << "\n";
			//if(true_new_theta_min < 0)
				//std::cout << "WARNING WARNING WARNING ligand " << i << "\n";
		}
		else
		{
			if(true_new_theta_max*true_new_theta_max<theta_max*theta_max)
				theta_max = true_new_theta_max;

			if(true_new_theta_min*true_new_theta_min<theta_min*theta_min)
				theta_min = true_new_theta_min;

			//if(true_new_theta_max > 0)
				//std::cout << "WARNING WARNING WARNING ligand " << i << "\n";
			//if(true_new_theta_min <0)
				//std::cout << "WARNING WARNING WARNING ligand " << i << "\n";
		}
	}
	if(theta_min < 0)
		theta_min = 0.0;
	if(theta_max > 0)
		theta_max = 0.0;
	//std::cout << '\n' << "Lim_theta: " << theta_min << " " << theta_max << '\n';
	*thetaminmax << theta_min << "," << theta_max << '\n';
	Lim_theta.push_back(theta_min);
	Lim_theta.push_back(theta_max);

	return Lim_theta;
} 

void rotation(std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<int> *L_to_R, std::vector<int> *list_bridged_receptors, std::vector<int> *T_l, long *idum, int step, std::ofstream *thetaminmax)
{
	bool bk;
	float new_theta, rnd;
	std::vector<float> Lim_theta;

	int test;
	test = 0;

	Lim_theta = theta_minmax(X_r_x, X_r_y, X_l_x, X_l_y, X_CM_x, X_CM_y, X_CM_theta, L_to_R, list_bridged_receptors, step, thetaminmax);
	bk = true;
	while(bk == true)
	{	rnd = ran2(idum);
		new_theta = rnd*(Lim_theta[0]-Lim_theta[1]) + Lim_theta[1];
		bk = check_breakage(X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x, X_r_y, T_l, L_to_R, 0.0, 0.0, new_theta);
		test += 1;
		//std::cout << test << ", ";
	}
	//std::cout << "new_theta: " << new_theta << '\n';
	//std::cout << test << "\n" << "\n";
	X_CM_theta->push_back((X_CM_theta->at(X_CM_theta->size() - 1) + new_theta));
}

std::vector<int> Get_bridges(std::vector<int> *T_l, std::vector<int> *L_to_R)
{
	//std::cout << "J" << '\n';
	std::vector<int> list_bridged_receptors;

	for(int i = 0; i<T_l->size(); i++)
	{
		if(T_l->at(i) == 2)
			list_bridged_receptors.push_back(L_to_R->at(i));
	}

	return list_bridged_receptors;
}

int main()
{
	clock_t startTime = clock();
    std::ofstream trajectory;
    trajectory.open ("trajectory.txt");

    trajectory << "X_CM,Y_CM,THETA_CM,N_P" << '\n';

    std::ofstream trajectory_2;
    trajectory_2.open ("trajectory_2.txt");

    std::ofstream receptor;
    receptor.open ("recetpor.txt");

    std::ofstream ligand;
    ligand.open ("ligand.txt");

    std::ofstream run_info;
    run_info.open ("run_info.txt");

    std::ofstream Liste_Ligand;
    Liste_Ligand.open ("L_L.txt");

    std::ofstream r_1;
    r_1.open ("r_1.txt");
    std::ofstream r_2;
    r_2.open ("r_2.txt");
    std::ofstream r_3;
    r_3.open ("r_3.txt");

    std::ofstream depleted;
    depleted.open("depleted_SA.txt");
    depleted << "STEP,SA" << '\n';

    std::ofstream opti;
    opti.open("opti.txt");

    std::ofstream time_sim;
    time_sim.open("time.txt");
    time_sim << "TIME,NBrow" << '\n';

    opti << "REACTION,DIFFUSION,UPDATE_VERLET,UPDATE_AFFINITY,UPDATE_TOTAL" << '\n';

    r_1 << "TIME" << '\n';
    r_2 << "TIME" << '\n';
    r_3 << "TIME" << '\n';

    /*std::ofstream Liste_Receptor;
    Liste_Receptor.open ("L_R.txt");*/

    std::ofstream Affinity;
    Affinity.open ("Affinity.txt");
  
    std::ofstream CB;
    CB.open("CB.txt");

    std::ofstream ltr;
    ltr.open("L_to_R.txt");

    std::ofstream bl;
    bl.open("bridges.txt");

    std::ofstream thetaminmax;
    thetaminmax.open("theta_minmax.txt");
    thetaminmax << "THETA_MIN,THETA_MAX" << "\n";
	
    srand (time(NULL));
    long idum =  -1*(rand() % 888888888 + 111111111); //-111142450;-111121278;-111121278;-111140767; -111121051; 

    run_info << "idum: " << idum << '\n';
    run_info << "k_d: " << k_d << "  k_off: " << k_off << "  k_on: " << k_on << "  DT: " << DT << "  N_BRWN: " << N_BRWN << "  Taille: " << Taille_Systeme << "  N_R: " << N_R << "  R_Verlet: " << R_Verlet << '\n';
    run_info.flush();
    std::cout << "seed du GNPA: " << idum << '\n';
    
  	int pct = 0;
  	int accepted_iteration = 0;
  	int accepted_iteration_temp = 0;
  	int accepted_iteration_step = 0;
  	int accept = 0;

	std::vector<float> a_d_tot;
	std::vector<float> a_on_tot;
	std::vector<float> a_off_tot;

	std::vector<float> affinity;

	MATRIX_INT v_L;
	std::vector<int> ll_L;
	std::vector<int> bl_L;

	MATRIX_INT v_R;
	std::vector<int> ll_R;
	std::vector<int> bl_R;

	std::vector<int> hoc_L;
	std::vector<int> hoc_R;

	float R_Cell;
	int N_x;
	int N_y;
	int N_Cell;
	int nc;
	int ix;
	int iy;

	MATRIX_INT Neigh;

	std::vector<float> X_l_x;  // dans le repÃ¨re du CM
	std::vector<float> X_l_y;
	
	std::vector<float> X_CM_x; // dans le repÃ¨re du laboratoire
	std::vector<float> X_CM_y;
	std::vector<float> X_CM_theta; 

	std::vector<float> X_r_x; // dans le repÃ¨re du laboratoire
	std::vector<float> X_r_y; 

	std::vector<int> T_l;
	std::vector<int> T_r;
	std::vector<int> L_to_R;

	std::vector<int> list_bridged_receptors;
	std::vector<int> list_CB;
	
	float time_react;
	float time_tot;
	std::vector<float> time_vect;
	time_vect.push_back(0);
	
	int N_ponts;
    std::vector<int> ponts_vect;
    ponts_vect.push_back(0.0);
    int count;

	initialize_IAV(&X_l_x, &X_l_y, &T_l, &L_to_R, &idum);


	ligand << "X_l_x = [";
	for(int i = 0; i<X_l_x.size();i++)
	{
		ligand << X_l_x[i] << ',';
	}
	ligand << "]" << '\n' << "X_l_y =[";
	for(int i = 0; i<X_l_y.size();i++)
	{
		ligand << X_l_y[i] << ',';
	}	
	ligand << "]";

	int N_L = X_l_x.size();
	
	MATRIX_INT L_L; 
	L_L.resize(N_L, std::vector<int>(0));
	MATRIX_INT L_R;
	L_R.resize(N_R, std::vector<int>(0));
	MATRIX_INT L_L_old; 
	L_L_old.resize(N_L, std::vector<int>(0));

	MATRIX_FLOAT a_d;
	a_d.resize(N_L, std::vector<float>(0));
	MATRIX_FLOAT a_on;
	a_on.resize(N_L, std::vector<float>(0));
	MATRIX_FLOAT a_off;
	a_off.resize(N_L, std::vector<float>(0));

    X_CM_x.push_back(0.0);
    X_CM_y.push_back(0.0);
    X_CM_theta.push_back(0.0);

    initialize_SA(&X_r_x, &X_r_y, &T_r, &idum);

	receptor << "X_r_x = [";
	for(int i = 0; i<X_r_x.size();i++)
	{
		receptor << X_r_x[i] << ',';
	}
	receptor << "]" << '\n' << "X_r_y =[";
	for(int i = 0; i<X_r_y.size();i++)
	{
		receptor << X_r_y[i] << ',';
	}	
	receptor << "]";


	N_x = int(Taille_Systeme/R_Verlet); //D->R_Verlet
	N_y = N_x;
	N_Cell = N_x*N_y;

	R_Cell = Taille_Systeme/N_x;

	Neigh = Set_Neigh(N_Cell, N_x, N_y);
    
	v_L = New_List(R_Cell, N_Cell, N_L, N_x, N_y, &X_l_x, &X_l_y);
	ll_L = v_L[0];
	bl_L = v_L[1];
	hoc_L = v_L[2];

	v_R = New_List(R_Cell, N_Cell, N_R, N_x, N_y, &X_r_x, &X_r_y);
	ll_R = v_R[0];
	bl_R = v_R[1];
	hoc_R = v_R[2];


	Initialize_L_L(&L_L, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_R, ll_R, T_l, L_to_R, N_x, N_y, R_Cell, N_L, 0);
	Initialize_L_R(&L_R, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_L, ll_L, T_r, T_l, L_to_R, N_x, N_y, R_Cell, N_L, 0);


	Initialize_Affinity(&a_d, &L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &L_to_R, T_l, T_r, k_d, N_L, 0, 1);
	Initialize_Affinity(&a_on, &L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &L_to_R, T_l, T_r, k_on, N_L, 1, 1);
    Initialize_Affinity(&a_off, &L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &L_to_R, T_l, T_r, k_off, N_L, 2, 2);


    a_d_tot = Get_a_X_tot(&a_d, N_L);
	a_on_tot = Get_a_X_tot(&a_on, N_L);
	a_off_tot = Get_a_X_tot(&a_off, N_L);
	
	//std::cout << '\n';

	affinity = Get_a_tot(&a_d_tot, &a_on_tot, &a_off_tot, N_L);

  	//std::vector<float> X_init;
  	std::vector<float> displacement_L;
  	std::vector<float> x_L;
  	std::vector<float> y_L;
  	std::vector<float> x_init_L;
  	std::vector<float> y_init_L;
  	float max_L = 0.0;
  	for(int w=0; w<N_L;w++)
  	{
  		x_init_L.push_back(0.0);
  		y_init_L.push_back(0.0);
  		x_L.push_back(0.0);
  		y_L.push_back(0.0);
  		displacement_L.push_back(0.0);
  	}

  	float X_CM_x_init = 0.0;
  	float X_CM_y_init = 0.0;
  	float X_CM_theta_init = 0.0;
  	std::vector<int> check;
  	float theta;

  	std::vector<float> displacement_R;
  	std::vector<float> displacement_close_R;
  	int N_close_R = 0;
  	std::vector<float> x_init_R;
  	std::vector<float> y_init_R;
  	float max_R = 0.0;
  	for(int w=0; w<N_R;w++)
  	{
  		x_init_R.push_back(0.0);
  		y_init_R.push_back(0.0);
  		displacement_R.push_back(0.0);
  	}

  	for(int w=0; w < N_L; w++)
	{
		theta = X_CM_theta_init;
		x_init_L[w] = X_l_x[w]*cos(theta) - X_l_y[w]*sin(theta) + X_CM_x_init;
	    y_init_L[w] = X_l_x[w]*sin(theta) + X_l_y[w]*cos(theta) + X_CM_y_init;

		//X_init[w] = sqrt(pow(x,2) + pow(y,2));
	}
	for(int w=0;w < N_R; w++)
	{
		x_init_R[w] = X_r_x[w];
		y_init_R[w] = X_r_y[w];
	}

	std::vector<float> new_rCM;
	std::vector<int> bridges_therm;
	int NBrow;

/*
	Liste_Ligand << '\n' << "L_L: " << '\n';
    for(int i = 0; i < N_L; i++)
    {
        Liste_Ligand << "Ligne de Ligand " << i+1 << ": ";
        for(int k = 0; k < L_L[i].size(); k++)
        {
            Liste_Ligand << L_L[i][k] << ' ';
        }
        Liste_Ligand << '\n';
    }
    
    Liste_Ligand << '\n';*/

	int gen = 0;

    for (int i = 0; i < N_BRWN; ++i)
	{	
		
		//std::cout << "\n" << "STEP " << i << '\n';
		time_tot = 0.0;
		count = 0;
		N_ponts = 0;		

		float DT1;

		DT1=DT;

		 if(i==0)
		{
		    DT1=0.1;
		}
		//std::cout << "test 1" << '\n';
		//REACTIONS
		while(time_tot < DT1)
    	{

    		int r = Choose(affinity, &idum); //choix du type de rÃ©action
    	    time_react = Get_time_reaction(affinity[0], &idum);
    	    time_tot += time_react;
    	    if(time_tot < DT1)
    	    {
    	    	if(i == 0 && gen == 0)
    	    	{
    	    		list_bridged_receptors.clear();
    	    		list_bridged_receptors = Get_bridges(&T_l, &L_to_R);

					std::cout << "list_bridged_receptors main loop: ";
					for(int k = 0; k<list_bridged_receptors.size();k++)
					{
						std::cout << list_bridged_receptors.at(k) << ' ';
					}
					std::cout << '\n';

    	    		if(list_bridged_receptors.size() > 1 && gen == 0)
	    			{	
	    				list_CB.clear();
	    				list_CB = gen_GCH(&X_l_x, &X_l_y, &X_CM_x, &X_CM_y, &X_CM_theta, &L_to_R,&X_r_x, &X_r_y, &list_bridged_receptors, &list_CB);
    	    			gen = 1;
    	    			std::cout << '\n' << '\n' << "list CB generated" << '\n';
    	    		}
    	    	}
    	    	

    	    	if(list_CB.size() == 2 && list_CB[0] == list_CB[1])
    	    	{
    	    		std::cout << "WARNING WARNING WARNING";
    	    		return 0;
    	    	}
    	    
        	    Make_reaction(&a_d, &a_on, &a_off, a_d_tot, a_on_tot, a_off_tot, L_L, L_R, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &T_l, &T_r, &L_to_R, &list_bridged_receptors, &list_CB, r, &idum, i, &depleted);
 				

                a_d_tot.clear();
                a_on_tot.clear();
                a_off_tot.clear();
                affinity.clear();
        	    
        	    a_d_tot = Get_a_X_tot(&a_d, N_L);
                a_on_tot = Get_a_X_tot(&a_on, N_L);
            	a_off_tot = Get_a_X_tot(&a_off, N_L);
            	affinity = Get_a_tot(&a_d_tot, &a_on_tot, &a_off_tot, N_L);

        	    count += 1;


        	    if(i==0){
		    		//Affinity_th << affinity[2] << "," << affinity[3] << "\n";
        	    	N_ponts = 0;
        	    	for(int k = 0; k<T_l.size(); k++) //calcul du nombre de ponts
					{
					    if(T_l[k] == 2)
					    {
					        N_ponts += 1;
					    }
					}
					//std::cout << N_ponts << "\n";
					bridges_therm.push_back(N_ponts);
					//bridges_th << bridges_therm[bridges_therm.size()-1] << "\n";	
        	    }
    	    }
    	}
		//std::cout << "test 2" << '\n';
    	depleted.flush();

    	N_ponts = 0;
		for(int k = 0; k<T_l.size(); k++) //calcul du nombre de ponts
		{
		    if(T_l[k] == 2)
		    {
		        N_ponts += 1;
		    }
		}
		//std::cout << N_ponts << '\n';
		ponts_vect.push_back(N_ponts);	

		if(i==0)
    	{//calcul dynamique de NBrow lors de la thermalisation
    		float mean_bridges=0;
    		float sum = std::accumulate(std::begin(bridges_therm), std::end(bridges_therm), 0.0);
			mean_bridges =  sum / (bridges_therm.size());

			run_info << " N_bridges: " << mean_bridges; 
			std::cout << " N_bridges: " << mean_bridges << '\n';

			float dXCM = sqrt(M_PI*4.9)*D/mean_bridges; 
			float DT_max = dXCM*dXCM/(D_parallel*10); 
			std::cout << DT_max << " " << DT << '\n';

			NBrow = ((int)(DT/DT_max) + 1)*5; //1000;
 
			std::cout << NBrow << '\n';

			run_info << " NBrow: " << NBrow << '\n';
   			run_info << '\n' << "thermalisation terminée en ";
			run_info << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';
			run_info.flush();
    	}

    	//DIFFUSION
		new_rCM = Sample_Conf_CGH(&X_r_x, &X_r_y, &X_l_x,&X_l_y, &list_bridged_receptors, &list_CB, &X_CM_x, &X_CM_y, &X_CM_theta,&L_to_R, &idum);
		X_CM_x.push_back(new_rCM[0]);
		X_CM_y.push_back(new_rCM[1]);

		bool check_bk = false;
		check_bk = check_breakage(&X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &T_l, &L_to_R, 0.0, 0.0, 0.0);
		if(check_bk == true)
		{
			std::cout << "\n" << "overextended bridge in step: " << i << "\n";
			std::cout << "x_CM = [" << X_CM_x[X_CM_x.size()-1] << "," << X_CM_y[X_CM_y.size()-1] << "]" << "\n";
			return 0;
		}
		//std::cout << "test 4" << '\n';
		//std::cout << "x_CM = " << X_CM_x[X_CM_x.size()-1] << ", " << X_CM_y[X_CM_y.size()-1] << ", " << X_CM_theta[X_CM_theta.size()-1] << "\n";
		


		rotation(&X_r_x, &X_r_y, &X_l_x, &X_l_y, &X_CM_x, &X_CM_y, &X_CM_theta, &L_to_R, &list_bridged_receptors, &T_l, &idum, i, &thetaminmax);	
		



		//std::cout << "x_CM = " << X_CM_x[X_CM_x.size()-1] << ", " << X_CM_y[X_CM_y.size()-1] << ", " << X_CM_theta[X_CM_theta.size()-1] << "\n";

		if(X_CM_x[X_CM_x.size()-1] != X_CM_x[X_CM_x.size()-1])
			return 0;
		if(X_CM_y[X_CM_y.size()-1] != X_CM_y[X_CM_y.size()-1])
			return 0;
		if(X_CM_theta[X_CM_theta.size()-1] != X_CM_theta[X_CM_theta.size()-1])
			return 0;

		int acc_rec = 0;
		acc_rec = diffusion_receptors(&X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &T_l, &L_to_R, &list_bridged_receptors, &idum, DT);


		//UPDATES 
		list_CB.clear();
	   	list_CB = gen_GCH(&X_l_x, &X_l_y, &X_CM_x, &X_CM_y, &X_CM_theta, &L_to_R,&X_r_x, &X_r_y, &list_bridged_receptors, &list_CB); //LIST CB
		
		new_rCM.clear();

		check.clear();
		//std::cout << "test 5" << '\n';
		int warning, bridged;
		warning = 0;
		for(int k =0;k<L_to_R.size();k++) //ERROR CHECK
		{
			if(L_to_R[k] == -1 && T_l[k] == 2)
			{	
				std::cout << "WARNING wrong ligand type in step " << i << "\n";
				warning +=1;
			}
			if(L_to_R[k] != -1 && T_r[L_to_R[k]-1] != 2)
			{	
				std::cout << "WARNING wrong receptor type in step " << i << "\n";
				warning +=1;
			}
		}
		for(int k = 0;k<list_bridged_receptors.size();k++) //LIST BRIDGED RECEPTORS
		{	
			bridged = 0;
			for(int l = 0;l<L_to_R.size();l++)
			{
				if(L_to_R[l] == list_bridged_receptors[k])
				{
					//std::cout << list_bridged_receptors[k] << " confirmed in L_to_R" << '\n';
					bridged += 1;
				}
			}
			if(bridged != 1)
				std::cout << "warning#2: receptor " << list_bridged_receptors[k] << "\n";
		}

		if(bridged != 1)
		{
			warning +=1;
			std::cout << "WARNING #2 wrong receptor type in step " << i << "\n";
			std::cout << T_r[24843];
		}

		if(warning>0)
			return 0;


		check_bk = false;
		check_bk = check_breakage(&X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &T_l, &L_to_R, 0.0, 0.0, 0.0);
		if(check_bk == true)
		{
			std::cout << "\n" << "overextended bridge 2nd test in step: " << i << "\n";
			std::cout << "x_CM = [" << X_CM_x[X_CM_x.size()-1] << "," << X_CM_y[X_CM_y.size()-1] << "]" << "\n";
			return 0;
		}

		//std::cout << "test E" << "\n";
		for(int w=0; w < N_L; w++)
		{
			theta = X_CM_theta[X_CM_theta.size() - 1];
			x_L[w] = X_l_x[w]*cos(theta) - X_l_y[w]*sin(theta) + X_CM_x[X_CM_x.size() -1];
		    y_L[w] = X_l_x[w]*sin(theta) + X_l_y[w]*cos(theta) + X_CM_y[X_CM_y.size() -1];
			displacement_L[w] = pow(x_L[w]-x_init_L[w],2) + pow(y_L[w]-y_init_L[w],2);
		}

		//std::cout << "test F" << "\n";
		for (int w=0; w < N_R; w++)
		{
			if(pow(X_r_x[w]-X_CM_x[X_CM_x.size()-1],2)+pow(X_r_y[w]-X_CM_y[X_CM_y.size()-1],2) < 200*200 )
			{
				displacement_R[w] = pow(X_r_x[w]-x_init_R[w],2) + pow(X_r_y[w]-y_init_R[w],2);
				displacement_close_R.push_back(displacement_R[w]);
			}
		}

		//std::cout << "test G" << '\n';
		N_close_R = displacement_close_R.size();

		max_L = sqrt(*max_element(displacement_L.begin(), displacement_L.end()));
		max_R = sqrt(*max_element(displacement_close_R.begin(), displacement_close_R.end()));

		//std::cout << "test H" << "\n";
		if(max_R + max_L > R_Verlet - D) //pow((X_init - X_current),2) >= pow(((D - R_Cell)*0.99),2) )
	    {

	    	//std::cout << "test A" << "\n";
	    	Update_Cell_L(&ll_L, &bl_L, &hoc_L, &X_l_x, &X_l_y, &X_CM_x, &X_CM_y, &X_CM_theta, &check, N_L, N_x, N_y, R_Cell, X_CM_x_init, X_CM_y_init, X_CM_theta_init);
	    	//std::cout << "test A 2" << "\n";
	    	Update_Cell_R(&ll_R, &bl_R, &hoc_R, &X_r_x, &X_r_y, &x_init_R, &y_init_R, N_x, N_y, R_Cell);
			//std::cout << "test B" << "\n";
			Update_Lists(&L_L, &L_R, &a_d, &a_on, &a_off, &a_d_tot, &a_on_tot, &a_off_tot, &affinity, Neigh, hoc_R, hoc_L, ll_L, ll_R, X_l_x, X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, N_x, N_y, R_Cell, N_L, T_l, T_r, &L_to_R, i);
			//std::cout << "test C" << "\n";
    		//Update_L_L(&L_L, &L_R, &L_L_old, N_L, &check, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_R, ll_R, T_l, L_to_R, N_x, N_y, R_Cell);

    		//Update_L_R(&L_L_old, &L_L, &L_R, &check);


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
			
			//std::cout << "test D" << "\n";
	    
	    }


    	else
    	{		
    		//std::cout << "test A bis" << "\n";
			Update_Affinity(&a_on, &a_off, &a_d, &a_d_tot, &a_on_tot, &a_off_tot, &affinity, &L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &L_to_R, T_l, T_r, N_L);
			//Update_Lists(&L_L, &L_R, &a_d, &a_on, &a_off, &a_d_tot, &a_on_tot, &a_off_tot, &affinity, Neigh, hoc_R, hoc_L, ll_L, ll_R, X_l_x, X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, N_x, N_y, R_Cell, N_L, T_l, T_r, &L_to_R, i);
			//std::cout << "test B bis" << "\n";
    	}


    //PRINTS
	    if(i % ((N_BRWN - 1)/200) == 0)
	    {
	    	trajectory << X_CM_x.at(X_CM_x.size() - 1) << "," << X_CM_y.at(X_CM_y.size() - 1) << "," << X_CM_theta.at(X_CM_theta.size() - 1) << "," << ponts_vect.at(ponts_vect.size() - 1) << '\n';
	      	trajectory.flush();		

		    for(int k = 0; k<list_bridged_receptors.size(); k++)
		    {
		    	bl << list_bridged_receptors[k] << ",";
		    }
		    bl << '\n';


		    for(int k = 0; k<list_CB.size(); k++)
		    {
		    	CB << list_CB[k] << ",";
		    }
		    CB << '\n';

		    for(int k = 0; k<L_to_R.size(); k++)
		    {
		    	ltr << L_to_R[k] << ",";
		    }
		    ltr << '\n';

			for(int i = 0; i<X_r_x.size();i++)
			{
				receptor << X_r_x[i] << ',';
			}
			receptor << "]" << '\n' << "X_r_y =[";
			for(int i = 0; i<X_r_y.size();i++)
			{
				receptor << X_r_y[i] << ',';
			}	
			receptor << "],";

		    receptor.flush();
		    bl.flush();
		    ltr.flush();
			CB.flush();

	    }	


	    //NETTOYAGE DES LISTES X_CM
	    /*X_CM_x.erase(X_CM_x.begin(),X_CM_x.begin() + (NBrow - 10));
   	    X_CM_y.erase(X_CM_y.begin(),X_CM_y.begin() + (NBrow - 10));
	    X_CM_theta.erase(X_CM_theta.begin(),X_CM_theta.begin()+ (NBrow - 10));*/

	    if(i % ((N_BRWN-1)/100) == 0)
	    {
	    	std::cout << pct << " %" << '\n';
	    	pct += 1;
	    }
		
	}	

	run_info << "run terminé en " << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds.";
	time_sim << double( clock() - startTime ) / (double)CLOCKS_PER_SEC << "," << NBrow;
	run_info.flush();

	depleted.close();
	trajectory.close();
	trajectory_2.close();
	receptor.close();
	ligand.close();
	run_info.close();
	Liste_Ligand.close();
	Affinity.close();
     
    return 0;
}
