/* ----------------------------------------------------------------------
    Ray Tracing for Many Spheres
    
    This is a Monte Carlo code to find the Radiation Distribution Factors 
    or View Factors in groups of uniform-sized spheres. 
    
    This code and more information can be found at:    
    https://github.com/ef-johnson/Ray-Tracing-Many-Spheres
    
    Written by Evan Johnson
    
    This code is provided as-is, as a tool to help other researchers, but
    without any guarantees. It is provided under the GNU General Public
    License, v3.0. 
     
------------------------------------------------------------------------- */


#include <math.h>
#include "mpi.h"
#include <iostream>
#include <cmath>
#include <random>
#include <cstdio>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <sstream>
#include <cstring>


using namespace std;





//------------------------------- MAIN -----------------------------------------------------

int main(int argc, char *argv[]) {

	// Start MPI
	MPI_Status status;
	MPI::Init(argc, argv);
	int n_procs = MPI::COMM_WORLD.Get_size();
	int rank = MPI::COMM_WORLD.Get_rank();

	// Frequently Changed Parameters

	unsigned long int n_photons;
	char char_n_photons[4]; // char strings need one more than the number of digits! More is ok though
	double vol_frac; 		
	int n_particles; 
	double abs;
	double abs_wall;
	int LTS_ts;
	char char_print_dist_vs_RDF[2];
	char char_print_all_RDF[2];
	double radius;
	double x_hi, x_lo, y_hi, y_lo;
	char is_wall_char;
	bool is_wall; 

if (rank == 0){


// Asks all of the questions to initiate the run:
cout << endl << "Welcome to the Monte-Carlo code for finding View Factors and Radiative Transfer within a particle bed." << endl << endl;
cout << "Enter information - n_particles and vol_frac must match the pos file!" << endl << endl ;

cout << "Enter number of photons in scientific notation (ex: 1e5): ";
cin >> char_n_photons;
cout << "Enter nominal volume fraction (must be in format 0.XX): ";
cin >> vol_frac;
cout << "Enter number of particles: ";
cin >> n_particles;
cout << "Enter absorptivity of the particles (format 0.XX): ";
cin >> abs;
cout << "Enter LIGGGHTS timestep number (format XXXXXX): ";
cin >> LTS_ts;
cout << "Enter radius of the particles (meters): ";
cin >> radius;
cout << "Do you want to export the matrix of distance vs. RDF? (Y/N): ";
cin >> char_print_dist_vs_RDF;
cout << "Do you want to export square matrix of all RDF's? (Y/N): ";
cin >> char_print_all_RDF;
cout << "Would you like to model the wall at Z=0? (Y/N): ";
cin >> is_wall_char;


// If this simulation is for PW:
if (is_wall_char == 'Y'){


cout << "Enter absorptivity of the wall (format 0.XX): ";
cin >> abs_wall;
cout << "Enter min x value of wall. (Photons hitting Z=0 wall outside of this will be ignored): ";
cin >> x_lo;
cout << "Enter max x value of wall. (Photons hitting Z=0 wall outside of this will be ignored): ";
cin >> x_hi;
cout << "Enter min y value of wall. (Photons hitting Z=0 wall outside of this will be ignored): ";
cin >> y_lo;
cout << "Enter max y value of wall. (Photons hitting Z=0 wall outside of this will be ignored): ";
cin >> y_hi;

}
else
{

abs_wall = 0;
x_lo = 0;
x_hi = 0;
y_lo = 0;
y_hi = 0;

}


} // end if rank = 0


// Broadcast from proc0 to all other processors: 
// (Address of the array, number of elements to pass, data type, root processor, communication group)
MPI::COMM_WORLD.Bcast(&char_n_photons, 4, MPI_CHAR, 0);
MPI::COMM_WORLD.Bcast(&vol_frac, 1, MPI_DOUBLE, 0);
MPI::COMM_WORLD.Bcast(&n_particles, 1, MPI_INT, 0);
MPI::COMM_WORLD.Bcast(&abs, 1, MPI_DOUBLE, 0);
MPI::COMM_WORLD.Bcast(&abs_wall, 1, MPI_DOUBLE, 0);
MPI::COMM_WORLD.Bcast(&LTS_ts, 1, MPI_INT, 0);
MPI::COMM_WORLD.Bcast(&char_print_dist_vs_RDF, 2, MPI_CHAR, 0);
MPI::COMM_WORLD.Bcast(&char_print_all_RDF, 2, MPI_CHAR, 0);
MPI::COMM_WORLD.Bcast(&is_wall_char, 2, MPI_CHAR, 0);
MPI::COMM_WORLD.Bcast(&radius, 1, MPI_DOUBLE, 0);
MPI::COMM_WORLD.Bcast(&x_lo, 1, MPI_DOUBLE, 0);
MPI::COMM_WORLD.Bcast(&x_hi, 1, MPI_DOUBLE, 0);
MPI::COMM_WORLD.Bcast(&y_lo, 1, MPI_DOUBLE, 0);
MPI::COMM_WORLD.Bcast(&y_hi, 1, MPI_DOUBLE, 0);

string str_n_photons(char_n_photons); // convert to string for later use in the file name

int base = char_n_photons[0] - 48; // converting from char to int requires subtracting 48 since it's from ASCII table!!
int exp  = char_n_photons[2] - 48;

n_photons = base * pow(10,exp);


// Set the wall to true or false (on all processors)
if (is_wall_char == 'Y'){
	is_wall = true;
}
else{
	is_wall = false;
}


if (rank == 0){
printf("\nCheck that the values were read in correctly and file names are right: \n\n");
}

printf(" --------------- Processor: %d -------------------\n", rank);
printf("n_photons: %s, vol_frac: %1.2f, n_particles: %d, abs: %f \n",char_n_photons, vol_frac, n_particles, abs);




	// Which outputs to print? all_RDF and/or dist_vs_RDF?
	bool print_dist_vs_RDF = true;
	bool print_all_RDF = true;

	//printf(" char dist RDF: %s",char_print_dist_vs_RDF);

	if (char_print_dist_vs_RDF[0] == 'Y'){ 
	print_dist_vs_RDF = true ;
	//printf("printing dist vs. RDF matrix\n");
	}
	else{ print_dist_vs_RDF = false;
	//printf("NOT printing dist vs. RDF matrix\n");
	}

	if (char_print_all_RDF[0] == 'Y'){ 
	print_all_RDF = true ;
	//printf("printing all RDF matrix\n");
	}
	else{ 
	print_all_RDF = false ;
	//printf("NOT printing all RDF matrix\n");
	}


	// Assign a cutoff value. Leave very high, but a user could change it to be lower if desired.
	double cutoff_radii = 1000;
	double cutoff_dist = cutoff_radii * radius;
		


	// Initialize strings
	ostringstream stream_abs;
	stream_abs << fixed << setprecision(2) << abs;
	string string_abs = stream_abs.str();
	
	ostringstream stream_vol_frac;
	stream_vol_frac << fixed << setprecision(2) << vol_frac;
	string string_vol_frac = stream_vol_frac.str();


	ostringstream stream_abs_wall;
	stream_abs_wall << fixed << setprecision(2) << abs_wall;
	string string_abs_wall = stream_abs_wall.str();

	ostringstream stream_LTS_ts;
	stream_LTS_ts << LTS_ts;
	string string_LTS_ts = stream_LTS_ts.str();

	ostringstream stream_n_particles;
	stream_n_particles << n_particles;
	string string_n_particles = stream_n_particles.str();

	ostringstream stream_rank;
	stream_rank << rank;
	string string_rank = stream_rank.str();



	// Position file (input) name
	string pos_filename = "pos/pos_" + string_n_particles + "_" + string_vol_frac + "_" + string_LTS_ts + ".txt";
	ifstream file(pos_filename);

	// home_id file (input) name
	string home_id_filename = "pos/home_id_" + string_n_particles + "_" + string_vol_frac + "_" + string_LTS_ts + ".txt";
	ifstream home_id_file(home_id_filename);

	// For the temporary files which write out results from each processor
	string cc_tmp = "tmp_data/dist_vs_RDF_" + string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_" + string_abs_wall  + "_" + string_LTS_ts + "_proc" +  string_rank + ".txt";
	string all_RDF_tmp = "tmp_data/all_RDF_" + string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_" + string_abs_wall + "_" + string_LTS_ts + "_proc" + string_rank + ".txt";
		string PW_RDF_tmp = "tmp_data/PW_RDF_" + string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_" + string_abs_wall + "_" + string_LTS_ts + "_proc" + string_rank + ".txt";

	// For final file after combining data from all processors
	string cc_filename =    "RDF_files/dist_vs_RDF_" + string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_" + string_abs_wall + "_" + string_LTS_ts + ".txt";
	string all_RDF_filename = "RDF_files/all_RDF_" + string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_" + string_abs_wall + "_" + string_LTS_ts + ".txt";
	string PW_RDF_filename = "RDF_files/PW_RDF_" + string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_" + string_abs_wall + "_" + string_LTS_ts + ".txt";
	
	// Print them out so user can inspect them, in case file isn't found
	cout << "File names:"  << endl;
	cout << pos_filename << endl;
	cout << cc_tmp << endl;
	cout << all_RDF_tmp << endl;
	cout << cc_filename << endl;
	cout << all_RDF_filename << endl << endl;

	
	// Output file streams

	// CC vs View Factor output
	ofstream out_file2;	
	if (print_dist_vs_RDF == true){ 	
	out_file2.open (cc_tmp);}
	
	// All View Factors output	
	ofstream out_file3;
	if (print_all_RDF == true){ 
	out_file3.open (all_RDF_tmp);}
	
	// PW RDF output	
	ofstream out_file6;
	if (is_wall == true){ 
	out_file6.open (PW_RDF_tmp);}






	// Random number generator
	random_device rd;   // non-deterministic generator
	mt19937 gen(rd());  // rd() is a random seed. Replace rd() with a number to get a consistent seed
	uniform_real_distribution<> dist(0,1); // distribute results between 0 and 1, inclusive

	// Initialize variables
	constexpr int n_read_col = 3;
	double pos[n_particles][n_read_col]; // x,y,z position and temps of each particle
	int home_id[n_particles]; // the ids of the particles which photons will be emitted from
	int n_home;
		
	double cc_dist[n_particles]; // center-center distance from home particle to all other particles
	double pi = atan(1)*4; // find pi

	int particles_simulated = 0;
	int particles_simulated_all[n_procs];

	clock_t t_before;
	clock_t t_after;


	// Declare variables for later
	double emiss_loc_rel[3];
	double emiss_loc[3];
	int current_emit_id;
	double alpha;
	double beta;


	// Geometry info:
	// V1: vector from emission location to center of particle being checked
	// V1_mag: distance from emission location to center of particle being checked
	// V2: Vector from emission location to the destination of the photon, in XYZ (global) coords 
	// V2_mag: Magnitude of V2 doesn't matter - but it's necessary to have a value for angle calculations
	// zeta: angle formed between V1 and V2
	// zeta_max: maximum angle between V1 and V2 to hit the sphere




if (rank == 0){

	// ----------------------- Read pos file with columns: x, y, z -------------------------------

	int counter = 0;
	if(file.is_open()){
			for (int row = 0; row<n_particles; row++){
				for(int col = 0; col < n_read_col; col++){
					file >> pos[row][col];
					counter++;
				}
			}
	}
	else {printf("\n\nCould not open pos/pos_xxxxxx file. Exiting.\n\n");
		exit (EXIT_FAILURE);
	}

	
		
		
	// --------- Read home_id file, a 1d array, where 1 = home (emitting) particle, 0 = not --------------------
		
	counter = 0;
	if(home_id_file.is_open()){
			for (int row = 0; row<n_particles; row++){
					home_id_file >> home_id[row];
					counter++;
			}
	}
	else {printf("\n\nCould not open pos/home_id_xxxxx.txt file. Exiting.\n\n");
		exit (EXIT_FAILURE);
	}
	
	
printf("\nResults from processor 0:\n");
	

} // end if rank = 0


// Broadcast: Address of the array, number of elements to pass, data type, root processor, communication group
MPI::COMM_WORLD.Bcast(&pos, n_particles*3, MPI_DOUBLE, 0);
MPI::COMM_WORLD.Bcast(&home_id, n_particles, MPI_INT, 0);

// *********************START HOME PARTICLE LOOP ************************  //

// Find the number of home particles to simulate
n_home = 0;
for (int i = 0; i < n_particles; i++)	{
	n_home = n_home + home_id[i];
}


if (rank == 0){
cout << endl << "Number of home particles to simulate: " << n_home << endl << endl;
}


// Distribute the home_id's among the processors, taking into account the case where n_home_id's is not divisible by n_procs
double n_home_proc_exact = (double)n_home / (double)n_procs ; // exact value of home particles per processor (convert from into to double).
int n_home_proc_nom = floor(n_home_proc_exact);
int n_remainder = n_home % n_procs;
int n_home_proc; // the number of home particles on this processor (will either be n_home_proc_nom or n_home_proc_nom+1)

// Assign n_home_proc to be either the nominal (lower value) or the nominal + 1 for the processors that are 0 up to whatever the n_remainder is
if (rank < n_remainder){
n_home_proc = n_home_proc_nom + 1;
}
else{
n_home_proc = n_home_proc_nom;
}
 
int home_id_proc[n_home_proc]; // home ids that this processor will handle


// Find the ids of all the home particles
// NOTE: because of the arrays starting at 0 in c++, the "home_num" is one less than the row number in the home_id file
// In other words, home_num of 17 shown on the terminal output from this program means that #18 is a home particle as 
// designated by the row number in Matlab/Octave (or whatever the home_id file is made in).
// This doesn't change anything, but it if you notice a mismatch between numbers output on the terminal and the original file, this is why. 
int home_id_all[n_home];
int j = 0;
for (int i = 0; i < n_particles; i++)	{
	if (home_id[i] == 1) {
		home_id_all[j] = i;
		j++;
	}
}

// Determine the home IDs that each processor will handle
// row_num is the row in the home_id_all vector where this processor will start, and in the for loop below it will be used along with n_home_proc to get the id values from home_id_all
int row_num;
// if the rank is in the first group of processors that have the extra, the row number is the rank times the (nominal plus one)
if (rank < n_remainder){
	row_num = rank * (n_home_proc_nom + 1); 
}
// if the rank is in the second group of processors (no extra home particle), the row number is the nominal times the rank, plus the number of extras (n_remainder)
else{
	row_num = rank * n_home_proc_nom + n_remainder; 
}

for (int i = 0; i < n_home_proc; i++)	{
	home_id_proc[i] = home_id_all[row_num];
	row_num++;
}



int home_num;
int home_counter = 0;

//  For each processor, check if this home particle id is on this processor's home_id_proc list
for (home_num = 0; home_num < n_particles; home_num++)	{

	// Check to see if this particle with home_num is on this processor's list (home_id_proc)
	bool home_bool = false;
	for (int i = 0; i < n_home_proc; i++)	{
		if (home_id_proc[i] == home_num) {
			home_bool = true; // The home_num particle IS in the list of home id's for this processor
		 }
	}
	
	// Only proceed with this home particle on this processor if home_bool is true (that is, if this home particle is on the home id list for this processor)
	if (home_bool == true){
	particles_simulated++;
	
	// If this is the first home particle simulated, get the clock tick to estimate the time the simulation will take.
	if (particles_simulated == 1){
		t_before = clock();
	}


		// -------------------- Create neighbor list for this home particle  ----------------------------------------
		// Originally the neighbor list was built to only include nearby particles, but to make it clearer/simpler 
		// (though with higher computation time), the neighbor list now has ALL the particles in it.

		double home_coords[3]= {pos[home_num][0], pos[home_num][1], pos[home_num][2] }; // coordinates of the current home particle
		double neighbor_list[n_particles+1][7]; // neighbor_list has one row for every particle, plus one row for the z=0 wall at the end
		for (int nei_id=0; nei_id<n_particles; nei_id++){

			double p_check[3]; // the particle being checked
			p_check[0] = pos[nei_id][0];
			p_check[1] = pos[nei_id][1];
			p_check[2] = pos[nei_id][2];


			// Find center-center distance between the "home" and the "check" particle
			cc_dist[nei_id] = {sqrt (pow(home_coords[0]-p_check[0],2) + pow(home_coords[1]-p_check[1],2) + pow(home_coords[2]-p_check[2],2) )};

			// Build the neighbor list with the id, cc_dist, XYZ coordinates, and space for the RDF calcs which come later
			neighbor_list[nei_id][0] = nei_id; // id of the particle
			neighbor_list[nei_id][1] = cc_dist[nei_id]; // center-center distance
			neighbor_list[nei_id][2] = p_check[0]; // x-pos
			neighbor_list[nei_id][3] = p_check[1]; // y-pos
			neighbor_list[nei_id][4] = p_check[2]; // z-pos
			neighbor_list[nei_id][5] = 0.0; // Total actual hits ("hit" and the closest V1_mag value for the photon)
			neighbor_list[nei_id][6] = 0.0; // PP RDF (Total absorptions / total photons sent for this home particle)


		} // end of building neighbor list

		// For the PW RDF: The last entry in the neighbor list is for the wall, so manually set this here
		// The entry for the wall is actually at the row number "n_particles" (not n_particles+1, since arrays start at 0)
		neighbor_list[n_particles][0] = -1; // id, where -1 indicates the wall
		neighbor_list[n_particles][1] = home_coords[2]; // Particle center to wall distance - this assumes wall at Z=0
		neighbor_list[n_particles][2] = 0; // x-pos (these don't matter for the wall)
		neighbor_list[n_particles][3] = 0; // y-pos (these don't matter for the wall)
		neighbor_list[n_particles][4] = 0; // z-pos (these don't matter for the wall)
		neighbor_list[n_particles][5] = 0; // Total photons absorbing into the wall
		neighbor_list[n_particles][6] = 0; // PW RDF (Total wall absorptions / total photons sent for this particle)




		// ---------------------------------------------------------------------------------------	

		
		// INFORMATION ABOUT GEOMETRY
		
		// XYZ is the global coord sys. It's easiest to visualize where Z is up, X is towards you, and Y is to the right.
		// x'y'z' is the local coord syst, which has the center on the sphere surface, 
		// at the location where the reflection happens, the incident location (inc_loc)
		// z' always points outward from sphere surface. 
		// x' is in the plane formed by the vectors: A) from the particle center to the incident location, and 
		// B) the projection of that vector onto the XY axis. 
		// y' follows based on the right hand rule
		
		// Beta is the angle from the Z axis down to the incident location.
		// Alpha is the azimuth angle from the X axis to the incident location, measured on the XY plane.
		// Eta is the angle from z' to the outgoing ray.
		// Gamma is the angle from the x' axis to the outgoing ray, measured on the x'y' plane.
		
		// For the initial photon emission, the emission location is chosen randomly anywhere on the sphere surface, 
		// by randomly generating the alpha and beta angles.
		// For diffuse emissions and reflections, eta and gamma are found randomly, in the outward facing hemisphere 
		// from the emission location. But the angles are generated in the x'y'z' coord sys, so they must be converted 
		// back to the XYZ coord system to proceed. This requires an "Euler angle" transformation.
		// Euler angles come in many different orientations, but implemented here are the "ZYZ" convention Euler angles.
		// In general this means a rotation around the Z axis, then a rotation around the (new) Y axis, then another 
		// rotation around the (new) Z axis. 
		// But here, with the angles defined as described above, only two rotations are needed, so the first rotation 
		// is just 0. So in this case, the coord sys is rotated by negative beta around the y' axis, and then by
		// negative alpha around the new z' axis (could call it the z" axis). The x'y'z' has been rotated to be in line with the XYZ axis. 
		
		
		
		
		//----------------- Monte Carlo for this home particle -------------------------------------
		
		// For this home particle, send photons, then check with each neighbor. Use Unsigned long int since i and n_photons can be very large.
		for (unsigned long int i = 0; i<n_photons; ++i){

			// This loops through as many reflections as needed until the photon absorbs.
			// The first time through, refl is set to true so that the original emission happens, 
			// and the emission location is set by alpha and beta
			// On successive runs through (until refl is set to false, representing an absorption or a lost 
			// photon) the emission location is calculated based on the incident location on the sphere surface.
			
			bool refl = true;
			bool first_run = true;
			bool overlapping = false; // set to true if the photon is released inside a neighboring sphere

			int refl_count = 0;

			while (refl == true){

				// On the first run through the reflection loop, the emission location is a random place on the surface, 
				// so alpha and beta must be chosen. 
				// if it's not the first run (i.e. if this is a reflection, not the initial emission), the emission location 
				// (and alpha and beta) are already known from the last iteration through the loop, so no need to find them here. 
				if (first_run == true) {

					// Choose angles for calculation of emission location and direction
					beta = acos(2*dist(gen) - 1); // Beta is angle down from vertical axis - choose a random angle between 0 and pi
					alpha = 2*pi*dist(gen); // Alpha is azimuth angle from x axis - choose a random angle between 0 and 2*pi


					// Find emission location from alpha and beta
					emiss_loc_rel[0] = radius*sin(beta)*cos(alpha);
					emiss_loc_rel[1] = radius*sin(beta)*sin(alpha);
					emiss_loc_rel[2] = radius*cos(beta);

					emiss_loc[0] = home_coords[0] + emiss_loc_rel[0];
					emiss_loc[1] = home_coords[1] + emiss_loc_rel[1];
					emiss_loc[2] = home_coords[2] + emiss_loc_rel[2];

					// Track the id of the particle where photon is emitted from. It starts here with the home particle.
					current_emit_id = home_num;

					// The first emission location has been set. Make first_run false since alpha and beta aren't 
					// chosen at random for the reflected photon.
					first_run = false;
					
				}
				
				

				// Choose angles for direction of emission

				// Initially, choose angles in temporary x'y'z' axes, 
				double eta = asin( sqrt(dist(gen)) ); // Eta is angle down from vertical axis - choose a random angle between 0 and pi/2
				double gamma = 2*pi*dist(gen); // Gamma is azimuth angle, measured from x-axis - choose a random angle between 0 and 2*pi


				// Calculate vector V2 (the ray, from the emission location) in the temporary x'y'z' axes (denoted by "prime")
				double V2_mag = 10; // magnitude doesn't matter
				double V2_prime[]= {V2_mag*sin(eta)*cos(gamma), V2_mag*sin(eta)*sin(gamma), V2_mag*cos(eta)};


				//Euler transform matrix ZYZ format, with the first Z rotation being zero.
				beta = -beta;
				alpha = -alpha;
				double a11 =cos(alpha)*cos(beta);
				double a12= sin(alpha);
				double a13=-cos(alpha)*sin(beta);
				double a21=-sin(alpha)*cos(beta);
				double a22= cos(alpha);
				double a23= sin(alpha)*sin(beta);
				double a31= sin(beta);
				double a32= 0;
				double a33= cos(beta);
				beta = -beta;
				alpha = -alpha;

				//Multiply the ray vector in x'y'z' (V2_prime) by the rotation matrix to the get vector in the new coord sys aligned with XYZ
				double V2[3];
				V2[0] = a11*V2_prime[0] + a12*V2_prime[1] + a13*V2_prime[2];
				V2[1] = a21*V2_prime[0] + a22*V2_prime[1] + a23*V2_prime[2];
				V2[2] = a31*V2_prime[0] + a32*V2_prime[1] + a33*V2_prime[2];




				//----------------- Check neighbors for hits and absorption --------------------------

				// Re-Initialize the ID and distance for the closest neighbor to emission location, for this photon.
				
				// NOTE: this code finds all the neighboring spheres which the photon would hit, and then it takes 
				// the closest one based on the distance between emission location and the center of the hit spheres - not the
				//  distance to the actual location where the photon strikes the sphere surface. 
				// This was done initially becasue it was thought that it was a simpler way which would yield the same result. 
				// But there is actually one rare case where it's not the same - when you have two spheres which are overlapping 
				// and the photon would hit both spheres in the overlap location, it is possible that the photon can strike the 
				// further sphere first, and then the closer sphere. This is a very rare case, unless the spheres are extremely 
				// soft and overlapping a high degree, which is not physically realistic. So it is left in the current state, 
				// where the actual hit sphere is judged based on the distance from emission location to the hit-spheres center locations. 
				
				
				// Keep track of the smallest distance from emission (or reflection) location to hit sphere's center (start excessively high)  
				double smallest_V1_mag = 1e9; 
				
				// keep track of the id of the smallest distance from emission location to hit sphere's center 
				// (start with a fake id, will get changed on first iteration)
				int closest_id = -99; 

				// Check to see which neighbors the photon would hit
				// Iterate from zero until the number of neighbors has been reached
				for (int nei_id = 0; nei_id < n_particles; nei_id++){

					// Don't check the particle that is currently emitting the photon
					if (nei_id != current_emit_id){
					
						// Get the xyz position of the particle being checked
						double p_check[] = {neighbor_list[nei_id][2],neighbor_list[nei_id][3],neighbor_list[nei_id][4]};

						// Calculate vector V1, from emission location to p_check center
						double V1[3] = {p_check[0]-emiss_loc[0], p_check[1]-emiss_loc[1], p_check[2]-emiss_loc[2]};
						double V1_mag = sqrt(pow(V1[0],2) + pow(V1[1],2) + pow(V1[2],2));
						
						// Special case: if the photon emits from inside the overlapping spheres, consider it an
						// automatic absorption. In reality, it would hit, but where would the incident location be
						// calculated? It's not physically realistic. It would likely reflect a couple times at
						// most until absorbing.
						
						// if V1_mag (the dist from emiss_loc to center of particle we're checking) is less than the radius
						// then the photon must be emitting from inside the neighboring sphere.
						if (V1_mag < radius) { 
							overlapping = true;
							closest_id = nei_id; // directly set this as the closest neighbor that hits, and setting overlapping=true prevents 
							// having to do further calculations for all the other neighbors due to the (overlapping==false) condition in the 
							// IF statement below. 
						}

						// Only check for hits if current neighbor is closer than the current closest hit,
						// AND only check if neighbor is closer than the cutoff distance.
						// AND only check if the photon is NOT being released inside of a neighboring sphere.
						if ((V1_mag < smallest_V1_mag) && (V1_mag < cutoff_dist) && overlapping == false){

							// Calculate angle between V1 and V2 (vector to center of particle being checked, and the photon path)
							double zeta = acos( (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2])/(V1_mag * V2_mag)); // find angle with dot product

							// Calculate maximum angle to qualify as a hit
							double zeta_max = asin(radius/V1_mag);

		
							// The photon will hit the neighbor particle if the angle is less than the max angle
							if (zeta < zeta_max){
							
								// If the V1_mag (emission location to neighbor_center) for this neighbor is less than the previous closest value,
								// save the value and id for the new closest hit

								smallest_V1_mag = V1_mag;
								closest_id = nei_id;
								
								
							} // end if zeta < zeta_max
						
						} // end if current neighbor is closer than the current closest hit
						
					} // end exclude the current emitting particle from the neighbor list
					
				} // end neighbor list for loop, which checks to see which neighbors the photon would hit




				//----------------- Check WALL for hits and absorption --------------------------
				
				
				// For the current home particle, if the current photon didn't hit any neighbors, check if it hits the wall or not
				if (closest_id == -99){
					
					// Only check wall for hit/absorption if a wall is specified by user
					if (is_wall == true) {
			
					// Wall is assumed to be under all of the particles (at Z=0)
					// Check for wall hit: any photon not hitting another particle AND going in the negative z-direction must hit the wall
						if (V2[2] < 0) { // photon hits wall at Z=0
					
					
							// We only want to count it as a wall-hit if it's hits within a square, defined by x_lo, x_hi, y_lo, y_hi.
							
							// Calculate the incident location where the photon crosses the Z=0 plane.
							// With the parametric equation x2i + y2j + z2k = (x1+s*x_vec)i + (y1+s*y_vec)j + (z1+s*z_vec)k
							// where x2 is the point we are looking for where the vector (x_vec i, y_vec j, z_vec k) passes through the z=0 plane
							// Here, x_vec, y_vec, z_vec is the V2[] vector, which has magnitude set with V2_mag,
							// and (x1, y1, z1) is the emission location where the vector leaves from, and s is the parametric variable.
							// For the special case of the plane Z=0, s can be solved with s=z1/z_vec	
								
							double inc_loc[3];						
							double s = -emiss_loc[2]/V2[2];
							inc_loc[0] = emiss_loc[0] + s*V2[0];
							inc_loc[1] = emiss_loc[1] + s*V2[1];
							inc_loc[2] = 0;
						
										
						
							// If it hits the Z=0 surface inside the bounds of x_lo, x_hi, etc, then the photon is checked 
							// for absorption. If it's outside, then the photon is lost.  
												
							if ( (inc_loc[0] < x_hi) && (inc_loc[0] > x_lo) && (inc_loc[1] < y_hi) && (inc_loc[1] > y_lo) ) {

								//  Check to see if photon is absorbed into the wall
								double rand_num = dist(gen);
								
								if (rand_num<abs_wall){ // photon is absorbed into the wall
									refl = false; // photon terminated
									neighbor_list[n_particles][5]++; // add one to the count of photons absorbed by the wall, 
								}

								else{ // Photon is reflected from the wall

									refl = true;
									current_emit_id = -1; // Set to -1, so that all the particles will be checked for hits. 
							
									// Set up the next run: incident location for this run is the emiss location for next run
									// Set the alpha and beta values. For the wall at Z=0, alpha and beta are both 0, 
									// since the 'prime' coord sys aligned with the global XYZ coord sys. In fact, this is why Z=0 was
									// chosen for the wall.  
									emiss_loc[0] = inc_loc[0];
									emiss_loc[1] = inc_loc[1];
									emiss_loc[2] = inc_loc[2];
						
									// Set alpha and beta to zero
									alpha = 0;
									beta = 0;
									
								} // end else, photon reflected from the wall
						
							} // end if photon hits the Z=0 plane inside of the domain we care about
						
							
							// Photon crosses the Z=0 plane, but outside the rectangle we care about, designated with x_lo, x_hi, etc.
							else{ 
								refl = false; // photon terminated (lost)
							}
						} // end if photon hits wall

						
						// If photon is not crossing the wall at Z=0, then the photon is done and didn't absorb into any particles or the
						// wall - it must have gone out one of the other sides of the domain.
						else{
						refl = false; // photon terminated (lost)
						}				
						
						} // end if there is a wall (is_wall = true)
						
						// If the photon didn't hit any particles and there is no wall (is_wall = false), terminate the photon 
						else{
						refl = false; // photon terminated (lost)
						}
						
					} // end if the photon hit NO particles
				



				//----------------- PARTICLE absorption or reflection --------------------------
				
				// If the photon did hit some neighbor(s)
				else  {

					// Special case: photon is released from inside the neighboring sphere -> an automatic absorption
					if (overlapping == true){
					
						refl = false; // photon terminated, absorbed into overlapping neighbor
						neighbor_list[closest_id][5]++; // add one to the count of photons absorbed by that neighbor				
					}
					
					// Photon not released inside an overlapping particle (almost always the case)
					else{ 
					
						//  Check to see if photon is absorbed
						double rand_num = dist(gen);
						if (rand_num<abs){
						
							refl = false; // photon terminated, absorbed
							neighbor_list[closest_id][5]++; // add one to the count of photons absorbed by that neighbor
						
						}

						else{ // Photon is reflected

							refl = true;
							
							// Identify the ID of the currently emitting particle. Set this for the next run through the reflection loop.
							current_emit_id = closest_id; 

							// Calculate the incident location based on the intersection of V2 with the closest hit-sphere's surface.

							// Get the center coordinates of the sphere that it actually hits
							double dest_ctr_coords[] = {pos[closest_id][0], pos[closest_id][1], pos[closest_id][2] };


							// Prep for quadratic:
							double a = V2[0]*V2[0] + V2[1]*V2[1] + V2[2]*V2[2];

							double dx = emiss_loc[0]-dest_ctr_coords[0];
							double dy = emiss_loc[1]-dest_ctr_coords[1];
							double dz = emiss_loc[2]-dest_ctr_coords[2];

							double b = 2*(dx*V2[0] + dy*V2[1] + dz*V2[2]);
							double c = -(radius*radius) + dx*dx + dy*dy + dz*dz;

							// Quadratic to find t:
							double t1 = (-b + sqrt( b*b - 4*a*c))/(2*a);
							double t2 = (-b - sqrt( b*b - 4*a*c))/(2*a);

							// Find the xyz locations where it hits
							double x1 = emiss_loc[0] + t1*V2[0];
							double y1 = emiss_loc[1] + t1*V2[1];
							double z1 = emiss_loc[2] + t1*V2[2];

							double x2 = emiss_loc[0] + t2*V2[0];
							double y2 = emiss_loc[1] + t2*V2[1];
							double z2 = emiss_loc[2] + t2*V2[2];

							double mag1 = sqrt ( (x1-emiss_loc[0])*(x1-emiss_loc[0]) + (y1-emiss_loc[1])*(y1-emiss_loc[1]) + (z1-emiss_loc[2])*(z1-emiss_loc[2]) );
							double mag2 = sqrt ( (x2-emiss_loc[0])*(x2-emiss_loc[0]) + (y2-emiss_loc[1])*(y2-emiss_loc[1]) + (z2-emiss_loc[2])*(z2-emiss_loc[2]) );

							// Vector will pass thru sphere twice, one going in, and one coming out
							// Take whichever point is closest to the emission location to be the incident location
							double inc_loc[3];
							if ( mag1 < mag2){
								inc_loc[0] = x1;
								inc_loc[1] = y1; // first point is closer,
								inc_loc[2] = z1;
							}
							else{
								inc_loc[0] = x2;
								inc_loc[1] = y2; // second point is closer,
								inc_loc[2] = z2;
							}

							// Find Alpha and Beta angles for next run through.
							// Because the photon is reflecting from the point where it hit the incident sphere's surface,
							// we need to find the alpha and beta angles of that point for the next run through the reflection loop

							double x_tmp = inc_loc[0] - dest_ctr_coords[0];
							double y_tmp = inc_loc[1] - dest_ctr_coords[1];
							double z_tmp = inc_loc[2] - dest_ctr_coords[2];

							// Calculate alpha. atan gives an angle in the range -pi/2 to +pi/2, so use atan2, which adjusts as necessary
							alpha = atan2( y_tmp, x_tmp );
							beta  = acos( z_tmp / radius );

							// Set the emission location for the next run as the incident location in this run
							emiss_loc[0] = inc_loc[0];
							emiss_loc[1] = inc_loc[1];
							emiss_loc[2] = inc_loc[2];

						} // end check if reflected
						
					} // end if not overlapping 

				} // else if this photon hit at least one neighbor

				refl_count++;

			} // end while refl = true loop (aka the reflection loop)

		} // end photon for loop



		// Find total photons absorbed, total missing photons, and total view factor for this particular home particle

		double total_abs = 0;
		for (int i=0;i<n_particles+1;i++){
			 neighbor_list[i][6] = neighbor_list[i][5] / n_photons; // Calculate the View Factor for each neighbor
			 total_abs = total_abs + neighbor_list[i][5];
		}


// Can use this warning if you should have no escapted photons (otherwise it's an annoying and unnecessary warning):

/*		 Warn if there is an escaped photon
		if (n_photons != total_abs){
			printf("Lost a photon!! %lu photons sent, but only %f absorbed!\n", n_photons,total_abs);
		}
		double total_vf = total_hits / n_photons;
*/




		// Print these center-center distances and view factors to the output text file
		if (print_dist_vs_RDF == true){ 
			for (int row=0;row<n_particles+1;row++){
			
				if (neighbor_list[row][1] < cutoff_dist){ // only print if less than the cutoff distance, saves a lot of file space.

					out_file2 << std::fixed << std::setprecision(0) <<  home_num; // "from" particle id
					out_file2 <<  " ";
					out_file2 << std::fixed << std::setprecision(0) <<  neighbor_list[row][0]; // "to" particle id
					out_file2 <<  " ";
					out_file2 << std::fixed << std::setprecision(8) <<  neighbor_list[row][1]; // center-center dist
					out_file2 <<  " "; 
					out_file2 << std::fixed << std::setprecision(8) <<  neighbor_list[row][6]; // RDF of PP or PW
					out_file2 <<  endl;
					}

			}
		}
		
		// Print all RDFs to the output file in a CSV format which has a square matrix of dimensions n_particles
		if (print_all_RDF == true){
			for (int to_particle=0;to_particle<n_particles+1;to_particle++){ // +1 is for the wall
				if (neighbor_list[to_particle][6]>0){
					out_file3 << std::fixed << std::setprecision(8) <<  neighbor_list[to_particle][6]; // high precision if there is an RDF
				}
				else{
					out_file3 << std::fixed << std::setprecision(0) <<  neighbor_list[to_particle][6]; // low precision if there is not an RDF (to save file space)
				}

				if (to_particle < n_particles){ 
					out_file3 <<  " ";
				}
				
			} // end for loop through all of the "to" particles for the current home particle
			
			out_file3 <<  endl;
		}
		
		// Print PW RDF vs. dist to the output file in a CSV format. Columns are: particle id, PW dist, PW RDF
		if (is_wall == true){ 
			for (int row=0;row<n_particles+1;row++){
			
				// Only print the RDF if it is a PW RDF (if the neighbor id is -1) 
				if (neighbor_list[row][0] == -1){

					out_file6 << std::fixed << std::setprecision(0) <<  home_num; // "from" particle id
					out_file6 <<  " ";

					out_file6 << std::fixed << std::setprecision(8) <<  neighbor_list[row][1]; // center-center dist
					out_file6 <<  " "; 

					out_file6 << std::fixed << std::setprecision(8) <<  neighbor_list[row][6]; // RDF of PP or PW
					out_file6 <<  endl;
					} // end if id is -1 (wall)
			}
		}



		// Outputs to terminal

		if (rank == 0){
			printf("Finished simulating home_num %d, home particle %d out of %d per processor. Sent %lu photons, %1.0f absorbed.\n", home_num, particles_simulated, n_home_proc, n_photons, total_abs);
		
			// Output the time estimate, only if it's the first home particle to be simulated
			double t_elapsed;
			if (particles_simulated == 1){

				t_after = clock();
				t_elapsed = (double)(t_after - t_before) / CLOCKS_PER_SEC; // the double is needed since the clock returns number of ticks, not seconds. "Double" turns the difference into a double.
				printf("Time estimate for finishing simulation: %f minutes, or %f hours \n\n", t_elapsed*n_home_proc/60, t_elapsed*n_home_proc/3600);
			
			}
			
		} // end if rank = 0

	} // end if home particle is in the home_id_proc list

	home_counter++;
	
}// end home particle loop

MPI_Barrier(MPI::COMM_WORLD);


if (rank == 0){

	// Open the files to print the entire output
	cout << endl << "Iterating through the temporary files: " << endl;

	// CC vs View Factor output
	ofstream out_file4;
	if (print_dist_vs_RDF == true){ 
		out_file4.open (cc_filename);
	}
	// All View Factors output
	ofstream out_file5;
	if (print_all_RDF == true){ 
		out_file5.open (all_RDF_filename);
	}
	// PW RDF output
	ofstream out_file6;
	if (is_wall == true){ 
		out_file6.open (PW_RDF_filename);
	}
	

// Iterate through the temporary output files, and combine them into one file
	for (int rnk = 0; rnk <n_procs;rnk++){

		ostringstream stream_rnk;
		stream_rnk << rnk;
		string string_rnk = stream_rnk.str();

		string inFile_cc = "tmp_data/dist_vs_RDF_" + string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_" + string_abs_wall + "_" + string_LTS_ts + "_proc" + string_rnk + ".txt";
		string inFile_all_RDF = "tmp_data/all_RDF_"+ string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_" + string_abs_wall + "_" + string_LTS_ts + "_proc" + string_rnk + ".txt";
		string inFile_PW_RDF = "tmp_data/PW_RDF_"+ string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_" + string_abs_wall + "_" + string_LTS_ts + "_proc" + string_rnk + ".txt";
		

		
		string line;

		// Print the dist vs RDF matrix if desired (not needed for Keff calc)
		if (print_dist_vs_RDF == true){
			ifstream cc_filestream(inFile_cc);
			cout << inFile_cc << endl;
			// Check to make sure the files to read are open
			    if (!cc_filestream.is_open() )
			    {
				cout << "Path Wrong (in dist vs. RDF)!!!! Check the file name and location." << endl;
				exit(EXIT_FAILURE);
			    }

			// Iterates through each line of the file, putting all data into "line" string. Then insert that line
			// into the overall output file.			
			while (getline(cc_filestream,line))
			{
			    istringstream iss(line);
			    string lineStream;
			    out_file4 << line << endl;
			}
		} // end if print_dist_vs_RDF = true

		// Print out the view all_RDF matrix
		if (print_all_RDF == true){
			ifstream all_RDF_filestream(inFile_all_RDF);
			cout << inFile_all_RDF << endl;
			
			// Check to make sure the files to read are open
		    if (!all_RDF_filestream.is_open() )
		    {
		        cout << "Path Wrong (in all_RDF)!!!! Check the file name and location." << endl;
		        exit(EXIT_FAILURE);
		    }
			while (getline(all_RDF_filestream,line))
			{
						 istringstream iss(line);
						 string lineStream;
						 out_file5 << line << endl;
			}


		} // end iterate through the output file written for each processor
				
				
				
		// Print the PW RDF matrix 
		if (is_wall == true){
			ifstream PW_RDF_filestream(inFile_PW_RDF);
			cout << inFile_PW_RDF << endl;
			// Check to make sure the files to read are open
			    if (!PW_RDF_filestream.is_open() )
			    {
				cout << "Path Wrong (in PW RDF)! Check the file name and location." << endl;
				exit(EXIT_FAILURE);
			    }

			// Iterate through each line of the file, putting all data into "line" string. Then insert that line
			// into the overall output file.
			
			while (getline(PW_RDF_filestream,line))
			{
			    istringstream iss(line);
			    string lineStream;
			    out_file6 << line << endl;
			}
		} // end if is_Wall = true
		
	} // end if rank = 0


out_file4.close();
out_file5.close();


printf("\nSuccessfully finished code execution.\n");
printf("Absorptivity: %1.2f \n", abs);
printf("Photons per particle: %lu \n",n_photons);
printf("Volume Fraction: %1.2f \n\n", vol_frac); 



} // end if rank = 0



MPI_Barrier(MPI::COMM_WORLD);
printf("Particles on proc%d: %d \n",rank,particles_simulated);


MPI::Finalize();

return 0;


}
