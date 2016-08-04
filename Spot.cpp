/***************************************************************************************/
/*                                   Spot.cpp

    This code produces a pulse profile once a set of parameter describing the star, 
    spectrum, and hot spot have been inputed.
    
    Based on code written by Coire Cadeau and modified by Sharon Morsink and 
    Abigail Stevens.
    
    PGxx refers to equation xx in Poutanen and Gierlinski 2003, arxiv: 0303084v1
    MLCBxx refers to equation xx in Morsink, Leahy, Cadeau & Braga 2007, arxiv: 0703123v2
    
    (C) Coire Cadeau, 2007; Source (C) Coire Cadeau 2007, all rights reserved.
*/
/***************************************************************************************/

/* Things to FIX 
1. Add Routine to compute look-up table of bending angles
2. Call Bending routine once in every latitude (instead of same for whole spot)
3. Change code so R_eq is input instead of R_spot (almost done!)
4. Change oblate code so it is based on R_eq
5. Add new shape function
6. Add look-up tables for dpsi, delta(t)
7. Dynamic memory allocation
*/


// TEST out GITHUB!!!


// INCLUDE ALL THE THINGS! 
// If you do not get this reference, see 
//       http://hyperboleandahalf.blogspot.ca/2010/06/this-is-why-ill-never-be-adult.html
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <exception>
#include <vector>
#include <string>
#include "OblDeflectionTOA.h"
#include "Chi.h"
#include "PolyOblModelNHQS.h"
#include "PolyOblModelCFLQS.h"
#include "SphericalOblModel.h"
#include "OblModelBase.h"
#include "Units.h"
#include "Exception.h"
#include "Struct.h"
#include "time.h"
#include <string.h>

// MAIN
int main ( int argc, char** argv ) try {  // argc, number of cmd line args; 
                                          // argv, actual character strings of cmd line args

  /*********************************************/
  /* VARIABLE DECLARATIONS AND INITIALIZATIONS */
  /*********************************************/
    
  std::ofstream out;      // output stream; printing information to the output file
  std::ofstream testout;  // testing output stream;
  std::ofstream param_out;// piping out the parameters and chisquared
    
  double incl_1(90.0),               // Inclination angle of the observer, in degrees
    incl_2(90.0),               // PI - incl_1; needed for computing flux from second hot spot, since cannot have a theta greater than 
    theta_1(90.0),              // Emission angle (latitude) of the first upper spot, in degrees, down from spin pole
    theta_2(90.0),              // Emission angle (latitude) of the second lower spot, in degrees, up from spin pole (180 across from first spot)
    mass,                       // Mass of the star, in M_sun
    rspot,                      // Radius of the star at the spot, in km
    mass_over_r, // Dimensionless mass divided by radius ratio
    omega,                      // Frequency of the spin of the star, in Hz
    req,                        // Radius of the star at the equator, in km
    bbrat(1.0),                 // Ratio of blackbody to Compton scattering effects, unitless
    ts(0.0),                    // Phase shift or time off-set from data; Used in chi^2 calculation
    spot_temperature(0.0),      // Inner temperature of the spot, in the star's frame, in keV
    phi_0_1,                    // Equatorial azimuth at the center of the piece of the spot that we're looking at
    phi_0_2,                    // Equatorial azimuth at the center of the piece of the second hot spot that we're looking at
    theta_0_1,                  // Latitude at the center of the piece of the spot we're looking at
    theta_0_2,                  // Latitude at the center of the piece of the second hot spot that we're looking at
    rho(0.0),                   // Angular radius of the inner bullseye part of the spot, in degrees (converted to radians)
    dphi(1.0),                  // Each chunk of azimuthal angle projected onto equator, when broken up into the bins (see numphi)
    dtheta(1.0),                // Each chunk of latitudinal angle, when broken up into the bins (see numtheta)
    phi_edge_2(0.0),            // Equatorial azimuth at the edge of the second spot at some latitude theta_0_2
    theta_edge_2(0.0),          // Latitudinal edge of lower (second) hot spot at -rho
    phishift,
    aniso(0.586),               // Anisotropy parameter
    Gamma1(2.0),                // Spectral index
    Gamma2(2.0),                // Spectral index
    Gamma3(2.0),                // Spectral index
    mu_1(1.0),                  // = cos(theta_1), unitless
    mu_2(1.0),                  // = cos(theta_2), unitless
    cosgamma,                   // Cos of the angle between the radial vector and the vector normal to the surface; defined in equation 13, MLCB
    Flux[NCURVES][MAX_NUMBINS], // Array of fluxes; Each curve gets its own vector of fluxes based on photons in bins.
    Temp[NCURVES][MAX_NUMBINS],
    trueSurfArea(0.0),          // The true geometric surface area of the spot
    E_band_lower_1(2.0),        // Lower bound of first energy band to calculate flux over, in keV.
    E_band_upper_1(3.0),        // Upper bound of first energy band to calculate flux over, in keV.
    E_band_lower_2(5.0),        // Lower bound of second energy band to calculate flux over, in keV.
    E_band_upper_2(6.0),        // Upper bound of second energy band to calculate flux over, in keV.
    background[NCURVES],
    T_mesh[30][30],             // Temperature mesh over the spot; same mesh as theta and phi bins, assuming square mesh
    chisquared(1.0),             // The chi^2 of the data; only used if a data file of fluxes is inputed
    distance(3.0857e22),        // Distance from earth to the NS, in meters; default is 10kpc
    B;                          // from param_degen/equations.pdf 2
   

  double E0, E1, E2, DeltaE;
    
  unsigned int NS_model(1),       // Specifies oblateness (option 3 is spherical)
    spectral_model(0),    // Spectral model choice (initialized to blackbody)
    beaming_model(0),     // Beaming model choice (initialized to isotropic)
    numbins(MAX_NUMBINS), // Number of time or phase bins for one spin period; Also the number of flux data points
    numphi(1),            // Number of azimuthal (projected) angular bins per spot
    numtheta(1),          // Number of latitudinal angular bins per spot
    numbands(NCURVES); // Number of energy bands;

  char out_file[256] = "flux.txt",    // Name of file we send the output to; unused here, done in the shell script
         out_dir[80],                   // Directory we could send to; unused here, done in the shell script
         T_mesh_file[100],              // Input file name for a temperature mesh, to make a spot of any shape
         data_file[256],                // Name of input file for reading in data
	  //testout_file[256] = "test_output.txt", // Name of test output file; currently does time and two energy bands with error bars
         filenameheader[256]="Run";

         
  // flags!
  bool incl_is_set(false),         // True if inclination is set at the command line (inclination is a necessary variable)
    	 theta_is_set(false),        // True if theta is set at the command line (theta is a necessary variable)
    	 mass_is_set(false),         // True if mass is set at the command line (mass is a necessary variable)
    	 rspot_is_set(false),        // True if rspot is set at the command line (rspot is a necessary variable)
    	 omega_is_set(false),        // True if omega is set at the command line (omega is a necessary variable)
    	 model_is_set(false),        // True if NS model is set at the command line (NS model is a necessary variable)
    	 datafile_is_set(false),     // True if a data file for inputting is set at the command line
    	 ignore_time_delays(false),  // True if we are ignoring time delays
    	 T_mesh_in(false),           // True if we are varying spot shape by inputting a temperature mesh for the hot spot's temperature
    	 normalize_flux(false),      // True if we are normalizing the output flux to 1
    	 E_band_lower_2_set(false),  // True if the lower bound of the second energy band is set
    	 E_band_upper_2_set(false),  // True if the upper bound of the second energy band is set
    	 two_spots(false),           // True if we are modelling a NS with two antipodal hot spots
    	 only_second_spot(false),    // True if we only want to see the flux from the second hot spot (does best with normalize_flux = false)
    	 pd_neg_soln(false);
		
  // Create LightCurve data structure
  class LightCurve curve, normcurve;  // variables curve and normalized curve, of type LightCurve
  class DataStruct obsdata;           // observational data as read in from a file


  /*********************************************************/
  /* READING IN PARAMETERS FROM THE COMMAND LINE ARGUMENTS */
  /*********************************************************/
    
    for ( int i(1); i < argc; i++ ) {
        if ( argv[i][0] == '-' ) {  // the '-' flag lets the computer know that we're giving it information from the cmd line
            switch ( argv[i][1] ) {
	    case 'a':  // Anisotropy parameter
	                sscanf(argv[i+1], "%lf", &aniso);
	                break;
	            
	    case 'b': // Blackbody ratio
	            	sscanf(argv[i+1], "%lf", &bbrat);
	            	break;

	    case 'd': //toggle ignore_time_delays (only affects output)
	                ignore_time_delays = true;
	                break;
	            
	    case 'D':  // Distance to NS
	            	sscanf(argv[i+1], "%lf", &distance);
	            	break;
	                
	    case 'e':  // Emission angle of the spot (degrees), latitude
	                sscanf(argv[i+1], "%lf", &theta_1);
	                theta_is_set = true;
	                break;

	    case 'f':  // Spin frequency (Hz)
	                sscanf(argv[i+1], "%lf", &omega);
	                omega_is_set = true;
	                break;

	    case 'g':  // Spectral Model, beaming (graybody factor)
	                sscanf(argv[i+1],"%u", &beaming_model);
	                break;
	                
	    case 'i':  // Inclination angle of the observer (degrees)
	                sscanf(argv[i+1], "%lf", &incl_1);
	                incl_is_set = true;
	                break;
	            
	    case 'I': // Name of input file
	            	sscanf(argv[i+1], "%s", data_file);
	            	datafile_is_set = true;
	            	break;
	                
	    case 'j':  // Flag for calculating only the second (antipodal) hot spot
	            	only_second_spot = true;
	            	break;
	          	          
	    case 'l':  // Time shift, phase shift.
	                sscanf(argv[i+1], "%lf", &ts);
	                break;

	    case 'k': // Background in low energy band (between 0 and 1)
	      sscanf(argv[i+1], "%lf", &background[2]);
	      break;

	    case 'K': // Background in high energy band (between 0 and 1)
	      sscanf(argv[i+1], "%lf", &background[3]);
	      break;
	          	          
	    case 'm':  // Mass of the star (solar mass units)
	                sscanf(argv[i+1], "%lf", &mass);
	                mass_is_set = true;
	                break;
	          
	    case 'n':  // Number of phase or time bins
	                sscanf(argv[i+1], "%u", &numbins);
	                break;
	                
	    case 'N': // Flag for not normalizing the flux output
	            	normalize_flux = true;
	            	break;

	    case 'o':  // Name of output file
	                sscanf(argv[i+1], "%s", out_file);
	                break;

	    case 'O':  case 'A': // This option controls the name of the output directory
	                sscanf(argv[i+1], "%s", out_dir);
	                break;
	          
	    case 'p':  // Angular Radius of spot (radians)
	                sscanf(argv[i+1], "%lf", &rho);
	                break;

	    case 'P':  // Number of phi bins 
	                sscanf(argv[i+1], "%u", &numphi);
	                break;	          

	    case 'q':  // Oblateness model (default is 1)
	                sscanf(argv[i+1], "%u", &NS_model);
	                model_is_set = true;
	                break;
	      	          
	    case 'r':  // Radius of the star at the equator(km)
	                sscanf(argv[i+1], "%lf", &req);
	                rspot_is_set = true;
	                break;

	    case 's':  // Spectral Model
	                sscanf(argv[i+1], "%u", &spectral_model);
			// 0 = blackbody
			// 1 = thin line for nicer
			break;
			
	    case 'S': // Number of Energy bands
	      sscanf(argv[i+1], "%u", &numbands);
	      break;

	    case 't':  // Number of theta bins 
	                sscanf(argv[i+1], "%u", &numtheta);
	                break;
	          
	    case 'T':  // Temperature of the spot, in the star's frame, in keV
	                sscanf(argv[i+1], "%lf", &spot_temperature);
	                break;
	            
	    case 'u': // Lower limit of first energy band, in keV
	            	sscanf(argv[i+1], "%lf", &E1);
	            	//E_band_lower_1_set = true;
	            	break;
	            
	    case 'U': // Upper limit of first energy band, in keV
	            	sscanf(argv[i+1], "%lf", &E2);
	            	//E_band_upper_1_set = true;
	            	break;
	                
	    case 'v': // Lower limit of second energy band, in keV
	            	sscanf(argv[i+1], "%lf", &E_band_lower_2);
	            	E_band_lower_2_set = true;
	            	break;
	            	
	    case 'V': // Upper limit of second energy band, in keV
	            	sscanf(argv[i+1], "%lf", &E_band_upper_2);
	            	E_band_upper_2_set = true;
	            	break;
	            
	    case 'x': // Scattering radius, in kpc
	            	sscanf(argv[i+1], "%lf", &E0);
	            	break;
	            
	    case 'X': // Scattering intensity, units not specified
	            	sscanf(argv[i+1], "%lf", &DeltaE);
	            	break;
	            	
	    case 'z': // Input file for temperature mesh
	            	sscanf(argv[i+1], "%s", T_mesh_file);
	            	T_mesh_in = true;
	            	break;
	            	
	    case '2': // If the user want two spots
	            	two_spots = true;
	            	break;
	            	
	            case '3': // Header for file name
	            	sscanf(argv[i+1],"%s", filenameheader);
	            	break;
	            
	            case '8': // Param_degen gave a negative solution
	            	pd_neg_soln = true;
	            	break;
	            
                case 'h': default: // Prints help
       	            std::cout << "\n\nSpot help:  -flag description [default value]\n" << std::endl
                              << "-a Anisotropy parameter. [0.586]" << std::endl
                              << "-b Ratio of blackbody flux to comptonized flux. [1.0]" << std::endl
                              << "-d Ignores time delays in output (see source). [0]" << std::endl
                              << "-D Distance from earth to star, in meters. [~10kpc]" << std::endl
                              << "-e * Latitudinal location of emission region, in degrees, between 0 and 90." << std::endl
                              << "-f * Spin frequency of star, in Hz." << std::endl
                              << "-g Graybody factor of beaming model: 0 = isotropic, 1 = Gray Atmosphere. [0]" << std::endl
		                      << "-i * Inclination of observer, in degrees, between 0 and 90." << std::endl
                              << "-I Input filename." << std::endl
		                      << "-j Flag for computing only the second (antipodal) hot spot. [false]" << std::endl
		                      << "-l Time shift (or phase shift), in seconds." << std::endl
		                      << "-m * Mass of star in Msun." << std::endl          
      	  	                  << "-n Number of phase or time bins. [128]" << std::endl
      	  	                  << "-N Flag for normalizing the flux. Using this sets it to true. [false]" << std::endl
		                      << "-o Output filename." << std::endl
		                      << "-O Name of the output directory." << std::endl
		                      << "-p Angular radius of spot, rho, in radians. [0.0]" << std::endl
		                      << "-q * Model of star: [3]" << std::endl
		                      << "      1 for Neutron/Hybrid quark star poly model" << std::endl
		                      << "      2 for CFL quark star poly model" << std::endl
		                      << "      3 for spherical model" << std::endl
		                      << "-r * Radius of star (at the equator), in km." << std::endl
		                      << "-s Spectral model of radiation: [0]" << std::endl
		                      << "      0 for bolometric light curve." << std::endl
		                      << "      1 for blackbody in monochromatic energy bands (must include T option)." << std::endl
		                      << "-t Number of theta bins for large spots. Must be < 30. [1]" << std::endl
		                      << "-T Temperature of the spot, in keV. [2]" << std::endl
		                      << "-u Low energy band, lower limit, in keV. [2]" << std::endl
		                      << "-U Low energy band, upper limit, in keV. [3]" << std::endl
		                      << "-v High energy band, lower limit, in keV. [5]" << std::endl
		                      << "-V High energy band, upper limit, in keV. [6]" << std::endl
		                      << "-x Scattering radius, in kpc." << std::endl
		                      << "-X Scattering intensity, units unspecified." << std::endl
		                      << "-z Input file name for temperature mesh." << std::endl
		                      << "-2 Flag for calculating two hot spots, on both magnetic poles. Using this sets it to true. [false]" << std::endl
		                      << "-3 File name header, for use with Ferret and param_degen." << std::endl
		                      << "-8 Sets a flag to indicate that param_degen gave a negative solution." << std::endl
		                      << " Note: '*' next to description means required input parameter." << std::endl
		                      << std::endl;
	                return 0;
            } // end switch	
        } // end if
    } // end for

    /***********************************************/
    /* CHECKING THAT THE NECESSARY VALUES WERE SET */
    /***********************************************/
    
    if( !( incl_is_set && theta_is_set
	    && mass_is_set && rspot_is_set
	    && omega_is_set && model_is_set ) ) {
        throw( Exception(" Not all required parameters were specified. Exiting.\n") );
        return -1;
    }
    
   
    /***************************************************************/
    /* READING TEMPERATURE VALUES FROM AN INPUT FILE, IF SPECIFIED */
    /***************************************************************/
	
    if ( T_mesh_in ) {
    	std::ifstream inStream(T_mesh_file);
    	std::cout << "** Tmeshfile name = " << T_mesh_file << std::endl;
    	unsigned int n(0);

    	if ( !inStream ) {
        	throw( Exception( "Temperature mesh input file could not be opened. \nCheck that the actual input file name is the same as the -z flag parameter value. \nExiting.\n") );
        	return -1;
        }
        std::string inLine;
        while ( std::getline(inStream, inLine) )  // getting numtheta from the file
        	n++;
        if (n != numtheta) {
        	throw( Exception( "Numtheta from mesh file does not match numtheta from command line input. \nExiting.\n") );
        	return -1;
        }
        numtheta = numphi = n;
        // go back to the beginning of the file
        inStream.clear();
        inStream.seekg (0, std::ios::beg);
        for ( unsigned int k(0); k < numtheta; k++ ) {
    		for ( unsigned int j(0); j < numphi; j++ ) {
    			inStream >> T_mesh[k][j];
    		}
    	}
    	inStream.close();
    	
    	// Printing the temperature mesh, for checking that it was read in correctly
    	std::cout << "Temperature Mesh of Spot:" << std::endl;
    	for ( unsigned int i(0); i < numtheta; i++ ) {
    		for ( unsigned int j(0); j < numphi; j++ ) {
    	    	std::cout << T_mesh[i][j] << "  ";
    		}
    		std::cout << std::endl;
    	}
    }
    
    /******************************************/
    /* SENSIBILITY CHECKS ON INPUT PARAMETERS */
    /******************************************/
    
    if ( numtheta < 1 ) {
        throw( Exception(" Illegal number of theta bins. Must be positive. Exiting.\n") );
        return -1;
    }
    if ( spot_temperature < 0.0 || rho < 0.0 || mass < 0.0 || rspot < 0.0 || omega < 0.0 ) {
        throw( Exception(" Cannot have a negative temperature, spot size, NS mass, NS radius, spin frequency. Exiting.\n") );
        return -1;
    }
    if ( rho == 0.0 && numtheta != 1 ) {
        std::cout << "\nWarning: Setting theta bin to 1 for a trivially sized spot.\n" << std::endl;
        numtheta = 1;
    }
   
    if ( numbins > MAX_NUMBINS || numbins <= 0 ) {
    	throw( Exception(" Illegal number of phase bins. Must be between 1 and MAX_NUMBINS, inclusive. Exiting.\n") );
    	return -1;
    }
    
 
 	
    /*****************************************************/
    /* UNIT CONVERSIONS -- MAKE EVERYTHING DIMENSIONLESS */
    /*****************************************************/

    mass_over_r = mass/(req) * Units::GMC2;
    incl_1 *= (Units::PI / 180.0);  // radians
    if ( only_second_spot ) incl_1 = Units::PI - incl_1; // for doing just the 2nd hot spot
    theta_1 *= (Units::PI / 180.0); // radians
    theta_2 = theta_1; // radians
    //rho *= (Units::PI / 180.0);  // rho is input in radians
    mu_1 = cos( theta_1 );
    mu_2 = mu_1; 
    mass = Units::cgs_to_nounits( mass*Units::MSUN, Units::MASS );
    req = Units::cgs_to_nounits( req*1.0e5, Units::LENGTH );
   
    omega = Units::cgs_to_nounits( 2.0*Units::PI*omega, Units::INVTIME );
    distance = Units::cgs_to_nounits( distance*100, Units::LENGTH );
	
     std::cout << "Dimensionless: Mass = " << mass << " Radius = " << req << " M/R = " << mass/req << std::endl; 

    /**********************************/
    /* PASS VALUES INTO THE STRUCTURE */
    /**********************************/    	
    
    curve.para.mass = mass;
    curve.para.mass_over_r = mass_over_r;
    curve.para.omega = omega;
    curve.para.radius = req;
    curve.para.req = req;
    curve.para.theta = theta_1;
    curve.para.incl = incl_1;
    curve.para.aniso = aniso;
    curve.para.Gamma1 = Gamma1;
    curve.para.Gamma2 = Gamma2;
    curve.para.Gamma3 = Gamma3;
    curve.para.temperature = spot_temperature;
    curve.para.ts = ts;
    curve.para.E_band_lower_1 = E_band_lower_1;
    curve.para.E_band_upper_1 = E_band_upper_1;
    curve.para.E_band_lower_2 = E_band_lower_2;
    curve.para.E_band_upper_2 = E_band_upper_2;
    curve.para.distance = distance;
    curve.numbins = numbins;
    //curve.para.rsc = r_sc;
    //curve.para.Isc = I_sc;

    //numphi = numtheta; // code currently only handles a square mesh over the hotspot
  
    curve.flags.ignore_time_delays = ignore_time_delays;
    curve.flags.spectral_model = spectral_model;
    curve.flags.beaming_model = beaming_model;
    curve.numbands = numbands;

   // Define the Spectral Model

    if (curve.flags.spectral_model == 0){ // NICER: Monochromatic Obs at E0=1keV
      curve.para.E0 = 1.0;
      curve.numbands = 1;
    }
    if (curve.flags.spectral_model == 1){ // NICER Line
      curve.para.E0 = E0; // Observed Energy in keV
      curve.para.E1 = E1; // Lowest Energy in keV in Star's frame
      curve.para.E2 = E2; // Highest Energy in keV
      curve.para.DeltaE = DeltaE; // Delta(E) in keV
    }

    // initialize background
    for ( unsigned int p(0); p < numbands; p++ ) background[p] = 0.0;

    /*************************/
    /* OPENING THE DATA FILE */
    /*************************/
   
    if ( datafile_is_set ) {	
      std::ifstream data; //(data_file);      // the data input stream
      data.open( data_file );  // opening the file with observational data
      char line[265]; // line of the data file being read in
      unsigned int numLines(0), i(0);
      if ( data.fail() || data.bad() || !data ) {
	throw( Exception("Couldn't open data file."));
	return -1;
      }
      // need to do the following to use arrays of pointers (as defined in Struct)
      for (unsigned int y(0); y < numbands; y++) {
	obsdata.t = new double[numbins];
	obsdata.f[y] = new double[numbins];
	obsdata.err[y] = new double[numbins];
      }
 
      /****************************************/
      /* READING IN FLUXES FROM THE DATA FILE */
      /****************************************/
    	
      while ( data.getline(line,265) ) {
	i = numLines;
	double get_f1;
	double get_err1;
	double get_f2;
	double get_err2;
			
	sscanf( line, "%lf %lf %lf %lf %lf", &obsdata.t[i], &get_f1, &get_err1, &get_f2, &get_err2 );
	obsdata.f[1][i] = get_f1;
	obsdata.err[1][i] = get_err1;
	obsdata.f[2][i] = get_f2;
	obsdata.err[2][i] = get_err2;
	// PG1808 data has 2 sigma error bars -- divide by 2!
	//obsdata.err[1][i] *= 2.0; // to account for systematic errors
	//obsdata.err[2][i] *= 2.0; // to account for systematic errors
	numLines++;
      } 

      if ( numLines != numbins ) {
	//throw (Exception( "Numbins indicated in command-line not equal to numbins in data file."));
	std::cout << "Warning! Numbins from command-line not equal to numbins in data file." << std::endl;
	std::cout << "Command-line numbins = " << numbins <<", data file numbins = " << numLines << std::endl;
	std::cout << "\t! Setting numbins = numbins from data file." << std::endl;
	numbins = numLines;
	curve.numbins = numLines;
	//return -1;
      }
       
      // Read in data file to structure "obsdata" (observed data)
      // f[1][i] flux in low energy band
      // f[2][i] flux in high energy band
      obsdata.numbins = numbins;

      data.close();
      //ts = obsdata.t[0]; // Don't do this if you want manually setting ts to do anything!!
      //obsdata.shift = obsdata.t[0];
      obsdata.shift = ts;

      //std::cout << "Finished reading data from " << data_file << ". " << std::endl;
    } // Finished reading in the data file
		
    /***************************/
    /* START SETTING THINGS UP */
    /***************************/ 

    // Change to computation of rspot!!!!

    // Calculate the Equatorial Radius of the star.
    //req = calcreq( omega, mass, theta_1, rspot );  // implementation of MLCB11

    //rspot = req;

    /*********************************************************************************/
    /* Set up model describing the shape of the NS; oblate, funky quark, & spherical */
    /*********************************************************************************/
	
    OblModelBase* model;
    if ( NS_model == 1 ) { // Oblate Neutron Hybrid Quark Star model
        // Default model for oblate neutron star
        model = new PolyOblModelNHQS( rspot, req,
		   		    PolyOblModelBase::zetaparam(mass,req),
				    PolyOblModelBase::epsparam(omega, mass, req) );
    }
    else if ( NS_model == 2 ) { // Oblate Colour-Flavour Locked Quark Star model
        // Alternative model for quark stars (not very different)
        model = new PolyOblModelCFLQS( rspot, req,
				     PolyOblModelBase::zetaparam(mass,rspot),
				     PolyOblModelBase::epsparam(omega, mass, rspot) );
    }
    else if ( NS_model == 3 ) { // Standard spherical model
        // Spherical neutron star
      //req = rspot;
      rspot = req;
        model = new SphericalOblModel( rspot );
    }
    else {
        throw(Exception("\nInvalid NS_model parameter. Exiting.\n"));
        return -1;
    }

    printf("R_Spot = %g; R_eq = %g \n", Units::nounits_to_cgs( rspot, Units::LENGTH ), Units::nounits_to_cgs( req, Units::LENGTH ));




    // defltoa is a structure that "points" to routines in the file "OblDeflectionTOA.cpp"
    // used to compute deflection angles and times of arrivals 
    // defltoa is the deflection time of arrival
    OblDeflectionTOA* defltoa = new OblDeflectionTOA(model, mass, mass_over_r, rspot); 
    // defltoa is a pointer (of type OblDeclectionTOA) to a group of functions

    // Values we need in some of the formulas.
    cosgamma = model->cos_gamma(mu_1);   // model is pointing to the function cos_gamma
    curve.para.cosgamma = cosgamma;


    /**********************************************************/
    /* Compute maximum deflection for purely outgoing photons */
    /**********************************************************/
	
    double  b_mid;  // the value of b, the impact parameter, at 90% of b_max
    curve.defl.b_max =  defltoa->bmax_outgoing(rspot); // telling us the largest value of b

    std::cout << "New M/R = " << mass_over_r << std::endl;
    std::cout << "New b_max/R = " << 1.0 / sqrt( 1.0 - 2.0 * mass_over_r ) << std::endl;
    std::cout << "z = " << 1.0 / sqrt( 1.0 - 2.0 * mass_over_r ) - 1.0 << std::endl;
    std::cout << "b_max/R = " << curve.defl.b_max/rspot << std::endl;

    curve.defl.psi_max = defltoa->psi_max_outgoing_u(curve.defl.b_max,rspot,&curve.problem); // telling us the largest value of psi
    std::cout << "psi_max = " << curve.defl.psi_max << std::endl;

    /********************************************************************/
    /* COMPUTE b VS psi LOOKUP TABLE, GOOD FOR THE SPECIFIED M/R AND mu */
    /********************************************************************/
	
    b_mid = curve.defl.b_max * 0.9; 
    // since we want to split up the table between coarse and fine spacing. 
    // 0 - 90% is large spacing, 90% - 100% is small spacing. b_mid is this value at 90%.
    curve.defl.b_psi[0] = 0.0; // definitions of the actual look-up table for b when psi = 0
    curve.defl.psi_b[0] = 0.0; // definitions of the actual look-up table for psi when b = 0
    
    for ( unsigned int i(1); i < NN+1; i++ ) { /* compute table of b vs psi points */
        curve.defl.b_psi[i] = b_mid * i / (NN * 1.0);
        curve.defl.psi_b[i] = defltoa->psi_outgoing_u(curve.defl.b_psi[i], rspot, curve.defl.b_max, curve.defl.psi_max, &curve.problem); 
	// calculates the integral
    }

    // For arcane reasons, the table is not evenly spaced.
    for ( unsigned int i(NN+1); i < 3*NN; i++ ) { /* compute table of b vs psi points */
        curve.defl.b_psi[i] = b_mid + (curve.defl.b_max - b_mid) / 2.0 * (i - NN) / (NN * 1.0); // spacing for the part where the points are closer together
        curve.defl.psi_b[i] = defltoa->psi_outgoing_u(curve.defl.b_psi[i], rspot, curve.defl.b_max, curve.defl.psi_max, &curve.problem); // referenced same as above
    }
    
    curve.defl.b_psi[3*NN] = curve.defl.b_max;   // maximums
    curve.defl.psi_b[3*NN] = curve.defl.psi_max;
    // Finished computing lookup table

    // Print out bending angles
    //std::cout << "Now print out the bending angles " << std::endl;
    //Bend(&curve,defltoa);

    /****************************/
    /* Initialize time and flux */
    /****************************/
	
    for ( unsigned int i(0); i < numbins; i++ ) {
        curve.t[i] = i / (1.0 * numbins);// + ts;  // defining the time used in the lightcurves
        for ( unsigned int p(0); p < numbands; p++ ) {
            Flux[p][i] = 0.0;                   // initializing flux to 0
            curve.f[p][i] = 0.0;
        }
    } 
    
   
    /************************************/
    /* LOCATION OF THE SPOT ON THE STAR */
    /* FOR ONE OR FIRST HOT SPOT        */
    /************************************/
    bool negative_theta(false); // possibly used for a hotspot that is asymmetric over the geometric pole

    /*********************************************/
    /* SPOT IS TRIVIALLY SIZED ON GEOMETRIC POLE */
    /*********************************************/
		
    if ( theta_1 == 0 && rho == 0 ) {
      for ( unsigned int i(0); i < numbins; i++ ) {
	for ( unsigned int p(0); p < numbands; p++ )
	  curve.f[p][i] = 0.0;
      }
      
      // Add curves, load into Flux array
      for (unsigned int i(0); i < numbins; i++) {
	for (unsigned int p(0); p < numbands; p++) {
	  Flux[p][i] += curve.f[p][i];
	}
      } // ending Add curves
    } // ending trivial spot on pole
	
    // This is the STANDARD CASE, and should be the one that is normally executed.
    
    //    std::cout << "Standard Case" << std::endl;
		
    if ( T_mesh_in ) {
      std::cout << "WARNING: code can't handle a spot asymmetric over the pole with a temperature mesh." << std::endl;
      spot_temperature = 2;
    }
    curve.para.temperature = spot_temperature;

    int pieces;
    //Does the spot go over the pole?
    if ( rho > theta_1){ // yes
      pieces=2;
    }
    else //no
      pieces=1;

    //std::cout << "pieces = " << pieces << std::endl;

    for (unsigned int p(0);p<pieces;p++){

      double deltatheta = 2.0*rho/numtheta;

      if (pieces==2){
	if (p==0)
	  deltatheta = (rho-theta_1)/numtheta;
	else
	  deltatheta = (2.0*theta_1)/numtheta;
      }
      
      std::cout << "numtheta=" << numtheta
		<< " theta = " << theta_1
		<< " delta(theta) = " << deltatheta << std::endl;
     
      // Looping through the mesh of the spot
      for (unsigned int k(0); k < numtheta; k++) { // Loop through the circles
	//for (unsigned int k(0); k < 6; k++) { // Loop through the circles
	std::cout << "k=" << k << std::endl;

	double thetak = theta_1 - rho + (k+0.5)*deltatheta; 
	double phi_edge, phij;

	if (pieces==2){
	  if (p==0){
	    thetak = (k+0.5)*deltatheta;
	    phi_edge = Units::PI;
	  }
	  else {
	    thetak = rho - theta_1 + (k+0.5)*deltatheta;
	  }
	}

	std::cout << "k=" << k << " thetak = " << thetak << std::endl;
	dphi = 2.0*Units::PI/(numbins*1.0);

	if ( (pieces==2 && p==1) || (pieces==1)){
	  double cos_phi_edge = (cos(rho) - cos(theta_1)*cos(thetak))/(sin(theta_1)*sin(thetak));
	  /*std::cout << "cos(rho) = " << cos(rho)
		    << " cos(th_1) = " << cos(theta_1)
		    << " cos(th_k) = " << cos(thetak)
		    << " Numerator = " << cos(rho) - cos(theta_1)*cos(thetak)
		    << " sin(th_1) = " << sin(theta_1)
		    << " sin(th_k) = " << sin(thetak)
		    << " Denominator = " << sin(theta_1)*sin(thetak)
		    << " cos(phi_e) = " << cos_phi_edge << std::endl;
	  */
	  if (  cos_phi_edge > 1.0 || cos_phi_edge < -1.0 ) 
	    cos_phi_edge = 1.0;

	  if ( fabs( sin(theta_1) * sin(thetak) ) > 0.0) { // checking for a divide by 0
	    phi_edge = acos( cos_phi_edge );   // value of phi (a.k.a. azimuth projected onto equatorial plane) at the edge of the circular spot at some latitude theta_0_2
	  // dphi = 2.0 * phi_edge / ( numphi * 1.0 );
	    /*
	    std::cout << "phi_edge = " << phi_edge 
		      << "numphi = " << numphi
		      << " delta(phi)= " << dphi << std::endl;
	    */
	  }
	  else {  // trying to divide by zero
	    throw( Exception(" Tried to divide by zero in calculation of phi_edge_2. Likely, theta_0_2 = 0. Exiting.") );
	    return -1;
	  }
	}
      
	numphi = 2.0*phi_edge/dphi;
	/*
	std::cout << "2.0*phi_edge = " << 2.0*phi_edge << " numphi=" << numphi << "  numphi*dphi=" << numphi*dphi << std::endl;
	std::cout << " Diff = " << 2.0*phi_edge - numphi*dphi << std::endl;
	*/
	phishift = 2.0*phi_edge - numphi*dphi;

	curve.para.dS = pow(rspot,2) * sin(thetak) * deltatheta * dphi;
	curve.para.theta = thetak;

	if (numtheta==1){
	  numphi=1;
	  phi_edge=0.0;
	  dphi=0.0;
	  phishift = 0.0;
	  curve.para.dS = 2.0*Units::PI * pow(rspot,2) * (1.0 - cos(rho));
	}
       

      for ( unsigned int j(0); j < numphi; j++ ) {   // looping through the phi divisions
	//  for ( unsigned int j(numphi-1); j < numphi; j++ ) {   // looping through the phi divisions
		  
	phij = -phi_edge + (j+0.5)*dphi;
	curve.para.phi_0 = phij;
	/*
	std::cout << "k=" << k 
		  << " theta-k = " << thetak
		  << " j=" << j 
		  << " phi-j = " << phij
		  << " dS = " << curve.para.dS << std::endl;
	*/
	if ( NS_model == 1 || NS_model == 2 )
	  curve.para.dS /= curve.para.cosgamma;


	if ( j==0){ // Only do computation if its the first phi bin - otherwise just shift
	  curve = ComputeAngles(&curve, defltoa); 
	  curve = ComputeCurve(&curve);
	}
	
	if ( curve.para.temperature == 0.0 ) { 
	  for ( unsigned int i(0); i < numbins; i++ ) {
	    for ( unsigned int p(0); p < numbands; p++ )
	      curve.f[p][i] = 0.0;
	  }
	}
	// Add curves, load into Flux array
	for ( unsigned int i(0); i < numbins; i++ ) {
	  int q(i+j);
	  if (q>=numbins) q+=-numbins;
	  for ( unsigned int p(0); p < numbands; p++ ) {
	    Flux[p][i] += curve.f[p][q];
	  }
	} // ending Add curves
      } // end for-j-loop
      // Add in the missing bit.

      
      if (phishift != 0.0){ // Add light from last bin, which requires shifting
      	for ( unsigned int i(0); i < numbins; i++ ) {
	  int q(i+numphi-1);
	  if (q>=numbins) q+=-numbins;
	  for ( unsigned int p(0); p < numbands; p++ ) {
	    Temp[p][i] = curve.f[p][q];
	  }
	}
	for (unsigned int p(0); p < numbands; p++ )
	  for ( unsigned int i(0); i < numbins; i++ ) 
	    curve.f[p][i] = Temp[p][i];	  
   	     
	curve = ShiftCurve(&curve,phishift);
       
      for( unsigned int p(0); p < numbands; p++ )
	for ( unsigned int i(0); i < numbins; i++ ) 
	  Flux[p][i] += curve.f[p][i]*phishift/dphi      ;
      }
      

      } // closing for loop through theta divisions
    } // end loop through pieces		 
	
     // End Standard Case

    /********************************************************************************/
    /* SECOND HOT SPOT -- Can handle going over geometric pole, but not well-tested */
    /********************************************************************************/

    if ( two_spots ) {
    	incl_2 = Units::PI - incl_1; // keeping theta the same, but changing inclination
    	curve.para.incl = incl_2;
    	curve.para.theta = theta_2;  // keeping theta the same, but changing inclination
    	cosgamma = model->cos_gamma(mu_2);
    	curve.para.cosgamma = cosgamma;
    		
    	if ( rho == 0.0 ) { // Default is an infinitesmal spot, defined to be 0 degrees in radius.
     		phi_0_2 = Units::PI;
     		theta_0_2 = theta_2;
    		curve.para.dS = trueSurfArea = 0.0; // dS is the area of the particular bin of spot we're looking at right now
    	}  
    	else {   // for a nontrivially-sized spot
    		dtheta = 2.0 * rho / (numtheta * 1.0);    // center of spot at theta, rho is angular radius of spot; starts at theta_2-rho to theta_2+rho
    		theta_edge_2 = theta_2 - rho;
    		trueSurfArea = 2 * Units::PI * pow(rspot,2) * (1 - cos(rho));
    	}
    		
    	/****************************************************************/
	/* SECOND HOT SPOT -- SPOT IS TRIVIALLY SIZED ON GEOMETRIC POLE */
	/****************************************************************/
		
	if ( theta_2 == 0 && rho == 0 ) {
	  for ( unsigned int i(0); i < numbins; i++ ) {
	    for ( unsigned int p(0); p < numbands; p++ )
	      curve.f[p][i] = 0.0;
	  }
	    
	  // Add curves, load into Flux array
	  for (unsigned int i(0); i < numbins; i++) {
	    for (unsigned int p(0); p < numbands; p++) {
	      Flux[p][i] += curve.f[p][i];
	    }
	  } // ending Add curves
	      	
	} // ending trivial spot on pole
    
    	/************************************************************/
	/* SECOND HOT SPOT -- SPOT IS SYMMETRIC OVER GEOMETRIC POLE */
	/************************************************************/
		
	else if ( theta_2 == 0 && rho != 0 ) {
	  // Looping through the mesh of the spot
	  for ( unsigned int k(0); k < numtheta; k++ ) {
	    theta_0_2 = theta_2 - rho + k * dtheta + 0.5 * dtheta;   // note: theta_0_2 changes depending on how we've cut up the spot
	    curve.para.theta = theta_0_2;
	    // don't need a phi_edge the way i'm doing the phi_0_2 calculation
	    dphi = Units::PI / numphi;
	    for ( unsigned int j(0); j < numphi; j++ ) {
	      phi_0_2 = -Units::PI/2 + dphi * j + 0.5 * dphi; // don't need to change phase because its second hot spot, because it's symmetric over the spin axis
	      if ( !T_mesh_in ) {
		T_mesh[k][j] = spot_temperature;
	      }
					
	      curve.para.dS = pow(rspot,2) * sin(fabs(theta_0_2)) * dtheta * dphi;
	      // Need to multiply by R^2 here because of my if numtheta=1 statement, 
	      // which sets dS = true surface area
					
	      if ( numtheta == 1 ) 
		curve.para.dS = trueSurfArea;
	      if ( NS_model == 1 || NS_model == 2)
		curve.para.dS /= curve.para.cosgamma;
	      curve.para.temperature = T_mesh[k][j];
	      curve.para.phi_0 = phi_0_2;
	      curve = ComputeAngles(&curve, defltoa);  // Computing the parameters it needs to compute the light curve; defltoa has the routines for what we want to do
	      curve = ComputeCurve(&curve);
	      
	      if (curve.para.temperature == 0.0) { 
		for (unsigned int i(0); i < numbins; i++) {
		  for (unsigned int p(0); p < numbands; p++)
		    curve.f[p][i] = 0.0;  // if temperature is 0, then there is no flux!
		}
	      }
	      // Add curves, load into Flux array
	      for (unsigned int i(0); i < numbins; i++) {
		for (unsigned int p(0); p < numbands; p++) {
		  Flux[p][i] += curve.f[p][i];
		}
	      } // ending Add curves
				
	    }
	  }
	} // ending symmetric spot over pole
    
    	/****************************************************************/
	/* SECOND HOT SPOT -- SPOT IS ASYMMETRIC OVER GEOMETRIC POLE */
	/****************************************************************/
		
	//THIS NEEDS TO BE FIXED IN THE SAME WAY THAT THE OTHER ASYMMETRIC ONE WAS!
		
	else if ( (theta_2 - rho) <= 0 ) {
	  // Looping through the mesh of the spot
	  for (unsigned int k(0); k < numtheta; k++) { // looping through the theta divisions
	    theta_0_2 = theta_2 - rho + 0.5 * dtheta + k * dtheta;   // note: theta_0_2 changes depending on how we've cut up the spot
	    if (theta_0_2 == 0.0) 
	      theta_0_2 = 0.0 + DBL_EPSILON;
	    curve.para.theta = theta_0_2;   // passing this value into curve theta, so that we can use it in another routine; somewhat confusing names
        
	    double cos_phi_edge = (cos(rho) - cos(theta_2)*cos(theta_0_2))/(sin(theta_2)*sin(theta_0_2));
	    if (  cos_phi_edge > 1.0 || cos_phi_edge < -1.0 ) 
	      cos_phi_edge = 1.0;
        
	    if ( fabs( sin(theta_2) * sin(theta_0_2) ) > 0.0) { // checking for a divide by 0
	      phi_edge_2 = acos( cos_phi_edge );   // value of phi (a.k.a. azimuth projected onto equatorial plane) at the edge of the circular spot at some latitude theta_0_2
	      dphi = 2.0 * phi_edge_2 / ( numphi * 1.0 );
	    }
	    else {  // trying to divide by zero
	      throw( Exception(" Tried to divide by zero in calculation of phi_edge_2. Likely, theta_0_2 = 0. Exiting.") );
	      return -1;
	    }

	    for ( unsigned int j(0); j < numphi; j++ ) {   // looping through the phi divisions
	      phi_0_2 = Units::PI -phi_edge_2 + 0.5 * dphi + j * dphi;   // phi_0_2 is at the center of the piece of the spot that we're looking at; defined as distance from the y-axis as projected onto the equator
	      if ( !T_mesh_in ) {
		T_mesh[k][j] = spot_temperature;
	      }
			 	
	      curve.para.dS = pow(rspot,2) * sin(fabs(theta_0_2)) * dtheta * dphi;
	      // Need to multiply by R^2 here because of my if numtheta=1 statement, 
	      // which sets dS = true surface area
            	
	      if( numtheta == 1 ) 
		curve.para.dS = trueSurfArea;
	      if ( NS_model == 1 || NS_model == 2 )
		curve.para.dS /= curve.para.cosgamma;
	      curve.para.temperature = T_mesh[k][j];
	      curve.para.phi_0 = phi_0_2;
	      
	      curve = ComputeAngles(&curve, defltoa);  // Computing the parameters it needs to compute the light curve; defltoa has the routines for what we want to do
	      curve = ComputeCurve(&curve);
	        	
	      if ( curve.para.temperature == 0.0 ) { 
		for ( unsigned int i(0); i < numbins; i++ ) {
		  for ( unsigned int p(0); p < numbands; p++ )
		    curve.f[p][i] = 0.0;
		}
	      }
	      // Add curves, load into Flux array
	      for ( unsigned int i(0); i < numbins; i++ ) {
		for ( unsigned int p(0); p < numbands; p++ ) {
		  Flux[p][i] += curve.f[p][i];
		}
	      } // ending Add curves
       	
	    } // closing for loop through phi divisions
	  } // closing for loop through theta divisions
	} // ending antisymmetric spot over pole

	/***********************************************************/
	/* SECOND HOT SPOT -- SPOT DOES NOT GO OVER GEOMETRIC POLE */
	/***********************************************************/
    
	else {
	  // Looping through the mesh of the second spot
	  printf("Working on the 2nd spot!!!\n");

	  for ( unsigned int k(0); k < numtheta; k++ ) { // looping through the theta divisions
	    theta_0_2 = theta_2 - rho + 0.5 * dtheta + k * dtheta;   // note: theta_0_2 changes depending on how we've cut up the spot
	    curve.para.theta = theta_0_2;   // passing this value into curve theta, so that we can use it in another routine; somewhat confusing names
	    //std::cout //<< "\n " k = " << k << ", theta_0_2 = " << theta_0_2 * 180.0 / Units::PI << std::endl;
       
	    double cos_phi_edge = (cos(rho) - cos(theta_2)*cos(theta_0_2))/(sin(theta_2)*sin(theta_0_2));
	    if ( cos_phi_edge > 1.0 || cos_phi_edge < -1.0 ) 
	      cos_phi_edge = 1.0;
       
	    if ( fabs( sin(theta_2) * sin(theta_0_2) ) > 0.0 ) { // checking for a divide by 0
	      phi_edge_2 = acos ( cos_phi_edge );   // value of phi (a.k.a. azimuth projected onto equatorial plane) at the edge of the circular spot at some latitude theta_0_2
	      dphi = 2.0 * phi_edge_2 / ( numphi * 1.0 ); // want to keep the same dphi as before
	    }
	    else {  // trying to divide by zero
	      throw( Exception(" Tried to divide by zero in calculation of phi_edge_2. \n Check values of sin(theta_2) and sin(theta_0_2). Exiting.") );
	      return -1;
	    }     
        
	    for ( unsigned int j(0); j < numphi; j++ ) {   // looping through the phi divisions
	      phi_0_2 = Units::PI - phi_edge_2 + 0.5 * dphi + j * dphi;
	      if ( !T_mesh_in ) {
		T_mesh[k][j] = spot_temperature;
	      }
           
	      curve.para.temperature = T_mesh[k][j];
	      curve.para.phi_0 = phi_0_2;
	      curve.para.dS = pow(rspot,2) * sin(fabs(theta_0_2)) * dtheta * dphi; // assigning partial dS here
	      if( numtheta == 1 ) 
		curve.para.dS = trueSurfArea;
	      if ( NS_model == 1 || NS_model == 2 )
		curve.para.dS /= curve.para.cosgamma;
           			
	      curve = ComputeAngles(&curve, defltoa);  // Computing the parameters it needs to compute the light curve; defltoa has the routines for what we want to do
	      curve = ComputeCurve(&curve);            // Compute Light Curve, for each separate mesh bit
        	
	      if ( curve.para.temperature == 0.0 ) { 
		for ( unsigned int i(0); i < numbins; i++ ) {
		  for ( unsigned int p(0); p < numbands; p++ )
		    curve.f[p][i] += 0.0;
		}
	      }
	      // Add curves, load into Flux array
	      for ( unsigned int i(0); i < numbins; i++ ) {
		for ( unsigned int p(0); p < numbands; p++ ) {

		  printf("i=%d Flux1=%g Flux2=%g ",i, Flux[p][i], curve.f[p][i]);
		  Flux[p][i] += curve.f[p][i];
		  printf(" Flux(1+2)=%g \n",Flux[p][i]);
		}
	      } // closing Add curves
		      
	    } // closing for loop through phi divisions
	  } // closing for loop through theta divisions
	} // closing spot doesn't go over geometric pole
    		
    } // closing if two spots
	
            
    // You need to be so super sure that ignore_time_delays is set equal to false.
    // It took almost a month to figure out that that was the reason it was messing up.
     
    
    /*******************************/
    /* NORMALIZING THE FLUXES TO 1 */
    /*******************************/
	    
    // Normalizing the flux to 1 in low energy band. 
    if ( normalize_flux ) {
      normcurve = Normalize1( Flux, numbins );

      // Add background to normalized flux
      for ( unsigned int i(0); i < numbins; i++ ) {
	for ( unsigned int p(0); p < numbands; p++ ) {
	  Flux[p][i] = normcurve.f[p][i] + background[p];
	}
      } 
      
      // Renormalize to 1.0
      normcurve = Normalize1( Flux, numbins );

      // fun little way around declaring the Chi.cpp method Normalize as returning a matrix! so that we can still print out Flux
      for ( unsigned int i(0); i < numbins; i++ ) {
	for ( unsigned int p(0); p < numbands; p++ ) {
	  curve.f[p][i] = normcurve.f[p][i];
	}
      }  
    } // Finished Normalizing
    else{ // Curves are not normalized
      for ( unsigned int i(0); i < numbins; i++ ) {
	for ( unsigned int p(0); p < numbands; p++ ) {
	  curve.f[p][i] = Flux[p][i];
	}
      }  
    }
     
    /************************************************************/
    /* If data file is set, calculate chi^2 fit with simulation */
    /************************************************************/
	
    if ( datafile_is_set ) {
    	chisquared = ChiSquare ( &obsdata, &curve );
    }
    
    std::cout << "Spot: m = " << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN 
	      << " Msun, r = " << Units::nounits_to_cgs(rspot, Units::LENGTH )*1.0e-5 
	      << " km, f = " << Units::nounits_to_cgs(omega, Units::INVTIME)/(2.0*Units::PI) 
	      << " Hz, i = " << incl_1 * 180.0 / Units::PI 
	      << ", e = " << theta_1 * 180.0 / Units::PI 
	      << ", X^2 = " << chisquared 
	      << std::endl;    


    /*******************************/
    /* CALCULATING PULSE FRACTIONS */
    /*******************************/
    
    double avgPulseFraction(0.0); // average of pulse fractions across all energy bands
    double sumMaxFlux(0.0), sumMinFlux(0.0);
    double overallPulseFraction(0.0);
    for ( unsigned int j(0); j < numbands; j++ ) {
      sumMaxFlux += curve.maxFlux[j];
      sumMinFlux += curve.minFlux[j];
      avgPulseFraction += curve.pulseFraction[j];
    }
    overallPulseFraction = (sumMaxFlux - sumMinFlux) / (sumMaxFlux + sumMinFlux);
    avgPulseFraction /= (NCURVES);
	
    double U(0.0), Q(0.0), A(0.0);   // as defined in PG19 and PG20
    U = (1 - 2 * mass / rspot) * sin(incl_1) * sin(theta_1); // PG20
    Q = 2 * mass / rspot + (1 - 2 * mass / rspot) * cos(incl_1) * cos(theta_1); // PG20
    A = U / Q; // PG19
    B = rspot * sin(incl_1) * sin(theta_1) * omega / ( sqrt( 1 - ((2*mass) / rspot) ) ); // param_degen/equations.pdf 2    

    /********************************************/
    /* WRITING THE SIMULATION TO AN OUTPUT FILE */
    /********************************************/ 
    	
    out.open(out_file, std::ios_base::trunc);
    if ( out.bad() || out.fail() ) {
        std::cerr << "Couldn't open output file. Exiting." << std::endl;
        return -1;
    }

    if ( rho == 0.0 ) rho = Units::PI/180.0;

    out << "# Photon Flux for hotspot on a NS. \n"
        << "# R_sp = " << Units::nounits_to_cgs(rspot, Units::LENGTH )*1.0e-5 << " km; "
        << "# R_eq = " << Units::nounits_to_cgs(req, Units::LENGTH )*1.0e-5 << " km; "
        << "# M = " << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN << " Msun; "
        << "# Spin = " << Units::nounits_to_cgs(omega, Units::INVTIME)/(2.0*Units::PI) << " Hz \n"
        << "# Gravitational redshift 1+z = " << 1.0/sqrt(1.0 - 2*mass/rspot) << "\n"
        << "# Inclination angle = " << incl_1 * 180.0/Units::PI << " degrees \n"
        << "# Emission angle = " << theta_1 * 180.0/Units::PI << " degrees \n"
        << "# Angular radius of spot = " << rho << " radians \n"
        << "# Number of spot bins = " << numtheta << " \n"
        << "# Distance to NS = " << Units::nounits_to_cgs(distance, Units::LENGTH)*.01 << " m \n" 
        << "# Phase shift or time lag = " << ts << " \n"
      //        << "# Pulse fractions: bolo = " << curve.pulseFraction[0] <<", mono = " << curve.pulseFraction[1] << ", \n"
      //<< "#                  low E band = " << curve.pulseFraction[2] << ", high E band = " << curve.pulseFraction[3] << " \n"
      //<< "# Rise/Fall Asymm: bolo = " << curve.asym[0] <<", mono = " << curve.asym[1] << ", \n"
      //<< "#                  low E band = " << curve.asym[2] << ", high E band = " << curve.asym[3] << " \n"
      //<< "# B = " << B << " \n"
      //<< "# B/(1+z) = " << B * sqrt(1.0 - 2*mass/rspot) << "\n"
      //<< "# U = " << U << " \n"
      //<< "# Q = " << Q << " \n"
      //<< "# A = U/Q = " << A << " \n"
      //<< "# %(Amp-Bolo) = " << (A - curve.pulseFraction[0])/A * 100.0 << " \n"
      //<< "# %(Bolo-low) = " << (curve.pulseFraction[0]-curve.pulseFraction[2])/curve.pulseFraction[0] * 100.0 << " \n"
        << std::endl;
        
    if (datafile_is_set)
    	out << "# Data file " << data_file << ", chisquared = " << chisquared << std::endl;

    if ( T_mesh_in )
    	out << "# Temperature mesh input: "<< T_mesh_file << std::endl;
    else
    	out << "# Spot temperature, (star's frame) kT = " << spot_temperature << " keV " << std::endl;
	if ( NS_model == 1)
    	out << "# Oblate NS model " << std::endl;
    else if (NS_model == 3)
    	out << "# Spherical NS model " << std::endl;
    if ( beaming_model == 0 )
        out << "# Isotropic Emission " << std::endl;
    else
        out << "# Limb Darkening for Gray Atmosphere (Hopf Function) " << std::endl;
    if ( normalize_flux )
    	out << "# Flux normalized to 1 " << std::endl;
    else
    	out << "# Flux not normalized " << std::endl;
    
    /***************************************************/
    /* WRITING COLUMN HEADINGS AND DATA TO OUTPUT FILE */
    /***************************************************/
  
    out << "#\n"
    	<< "# Column 1: phase bins (0 to 1)\n";

    if (curve.flags.spectral_model==0){
      out << "# Column 2: Monochromatic Number flux (photons/(cm^2 s keV) measured at energy (at infinity) of " << curve.para.E0 << " keV\n";
      out << "#"
	  << std::endl;
      for ( unsigned int i(0); i < numbins; i++ ) {
        out << curve.t[i]<< "\t" ;
	for ( unsigned int p(0); p < numbands; p++ ) { 
            out << curve.f[p][i] << "\t" ;
        }
	out << i;
        out << std::endl;
      }
    }

    if (curve.flags.spectral_model==1){
      out << "# Column 2: Photon Energy (keV) in Observer's frame" << std::endl;
      out << "# Column 3: Number flux (photons/(cm^2 s) " << std::endl;
      out << "# " << std::endl;

      for ( unsigned int p(0); p < numbands; p++ ) { 
	for ( unsigned int i(0); i < numbins; i++ ) {
	  out << curve.t[i]<< "\t" ;
	  out << curve.para.E0 + p*curve.para.DeltaE << "\t";
	  out << curve.f[p][i] << "\t" << std::endl;
	  //out << i;
	  //out << std::endl;
	}

      }
    }



      //<< "# Column 4: Number flux (photons/(cm^2 s)) in the energy band " << E_band_lower_1 << " keV to " << E_band_upper_1 << " keV \n"
      //<< "# Column 5: Number flux (photons/(cm^2 s)) in the energy band " << E_band_lower_2 << " keV to " << E_band_upper_2 << " keV \n"

    
  


    out.close();
    
    delete defltoa;
    delete model;
    return 0;
} 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
