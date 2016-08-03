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
           phi_edge_1(0.0),            // Equatorial azimuth at the edge of the spot at some latitude theta_0_1
           phi_edge_2(0.0),            // Equatorial azimuth at the edge of the second spot at some latitude theta_0_2
           theta_edge_1(0.0),          // Latitudinal edge of upper hot spot at -rho
           theta_edge_2(0.0),          // Latitudinal edge of lower (second) hot spot at -rho
           aniso(0.586),               // Anisotropy parameter
           Gamma1(2.0),                // Spectral index
           Gamma2(2.0),                // Spectral index
           Gamma3(2.0),                // Spectral index
           mu_1(1.0),                  // = cos(theta_1), unitless
           mu_2(1.0),                  // = cos(theta_2), unitless
           cosgamma,                   // Cos of the angle between the radial vector and the vector normal to the surface; defined in equation 13, MLCB
           Flux[NCURVES][MAX_NUMBINS], // Array of fluxes; Each curve gets its own vector of fluxes based on photons in bins.
           trueSurfArea(0.0),          // The true geometric surface area of the spot
           E_band_lower_1(2.0),        // Lower bound of first energy band to calculate flux over, in keV.
           E_band_upper_1(3.0),        // Upper bound of first energy band to calculate flux over, in keV.
           E_band_lower_2(5.0),        // Lower bound of second energy band to calculate flux over, in keV.
           E_band_upper_2(6.0),        // Upper bound of second energy band to calculate flux over, in keV.
           T_mesh[30][30],             // Temperature mesh over the spot; same mesh as theta and phi bins, assuming square mesh
           chisquared(1.0),             // The chi^2 of the data; only used if a data file of fluxes is inputed
           distance(3.0857e22),        // Distance from earth to the NS, in meters; default is 10kpc
           B,                          // from param_degen/equations.pdf 2
           r_sc,                       // Scattering radius, distance of a cloud scattering the flux
           I_sc;                       // Scattering intensity, how much the cloud scatters the flux

    unsigned int NS_model(1),       // Specifies oblateness (option 3 is spherical)
                 spectral_model(0),    // Spectral model choice (initialized to blackbody)
                 beaming_model(0),     // Beaming model choice (initialized to isotropic)
                 numbins(MAX_NUMBINS), // Number of time or phase bins for one spin period; Also the number of flux data points
                 numphi(1),            // Number of azimuthal (projected) angular bins per spot
                 numtheta(1);          // Number of latitudinal angular bins per spot

	char out_file[256] = "flux.txt",    // Name of file we send the output to; unused here, done in the shell script
         out_dir[80],                   // Directory we could send to; unused here, done in the shell script
         T_mesh_file[100],              // Input file name for a temperature mesh, to make a spot of any shape
         data_file[256],                // Name of input file for reading in data
         testout_file[256] = "test_output.txt", // Name of test output file; currently does time and two energy bands with error bars
         param_out_file[256], // Writes the parameters ran for all runs that day
         spot_run_out_file_pos[256],
         spot_run_out_file_neg[256],
         spot_run_out_file[256],
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
    	 E_band_lower_1_set(false),  // True if the lower bound of the first energy band is set
    	 E_band_upper_1_set(false),  // True if the upper bound of the first energy band is set
    	 E_band_lower_2_set(false),  // True if the lower bound of the second energy band is set
    	 E_band_upper_2_set(false),  // True if the upper bound of the second energy band is set
    	 two_spots(false),           // True if we are modelling a NS with two antipodal hot spots
    	 only_second_spot(false),    // True if we only want to see the flux from the second hot spot (does best with normalize_flux = false)
    	 pd_neg_soln(false);
	
	// the following is used to make 'todaysdate' be today's date in string format, for filename writing purposes
/*	time_t tim;
	struct tm *ptr;
	char todaysdate[80];
	time( &tim );
	ptr = localtime( &tim );
	strftime(todaysdate,80,"%d-%b-%Y",ptr);*/
	
	//std::cout << "Output file of parameters for all runs today: " << param_out_file << std::endl;
	
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
	          
                case 'p':  // Angular Radius of spot (degrees)
	                sscanf(argv[i+1], "%lf", &rho);
	                break;
	          
	            case 'q':  // Oblateness model (default is 1)
	                sscanf(argv[i+1], "%u", &NS_model);
	                model_is_set = true;
	                break;
	      	          
                case 'r':  // Radius of the star at the spot(km)
	                sscanf(argv[i+1], "%lf", &rspot);
	                rspot_is_set = true;
	                break;

                case 's':  // Spectral Model
	                sscanf(argv[i+1], "%u", &spectral_model);
	                break;
	            
	            case 't':  // Number of theta bins 
	                sscanf(argv[i+1], "%u", &numtheta);
	                break;
	          
                case 'T':  // Temperature of the spot, in the star's frame, in keV
	                sscanf(argv[i+1], "%lf", &spot_temperature);
	                break;
	            
	            case 'u': // Lower limit of first energy band, in keV
	            	sscanf(argv[i+1], "%lf", &E_band_lower_1);
	            	E_band_lower_1_set = true;
	            	break;
	            
	            case 'U': // Upper limit of first energy band, in keV
	            	sscanf(argv[i+1], "%lf", &E_band_upper_1);
	            	E_band_upper_1_set = true;
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
	            	sscanf(argv[i+1], "%lf", &r_sc);
	            	break;
	            
	            case 'X': // Scattering intensity, units not specified
	            	sscanf(argv[i+1], "%lf", &I_sc);
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
		                      << "-p Angular radius of spot, rho, in degrees. [0.0]" << std::endl
		                      << "-q * Model of star: [3]" << std::endl
		                      << "      1 for Neutron/Hybrid quark star poly model" << std::endl
		                      << "      2 for CFL quark star poly model" << std::endl
		                      << "      3 for spherical model" << std::endl
		                      << "-r * Radius of star (at the spot), in km." << std::endl
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
    
    
	/**************************/
	/* WRITING THE FILE NAMES */
	/**************************/
	
	strcpy(param_out_file, "/Users/jasper/Documents/spot/run_data/");
	strcat(param_out_file, filenameheader);
	strcat(param_out_file, "_run_log.txt");
	
    strcpy(spot_run_out_file, "/Users/jasper/Documents/spot/contours/");
	strcat(spot_run_out_file, filenameheader);
	strcat(spot_run_out_file, "_spot_for_gnuplot.txt");
	
	strcpy(spot_run_out_file_pos, "/Users/jasper/Documents/spot/contours/");
	strcat(spot_run_out_file_pos, filenameheader);
	strcat(spot_run_out_file_pos, "_spot_for_gnuplot_pos.txt");
	
	strcpy(spot_run_out_file_neg, "/Users/jasper/Documents/spot/contours/");
	strcat(spot_run_out_file_neg, filenameheader);
	strcat(spot_run_out_file_neg, "_spot_for_gnuplot_neg.txt");
	
	
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
    
    if ( numtheta > 30 || numtheta < 1 ) {
        throw( Exception(" Illegal number of theta bins. Must be between 1 and 30, inclusive. Exiting.\n") );
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
    if ( !E_band_upper_1_set || !E_band_lower_1_set ) {
    	throw( Exception(" Need to set upper and lower bounds to energy band 1 for flux calculation. Exiting.\n") );
    	return -1;
    }
    if ( E_band_lower_1_set && E_band_upper_1_set && (E_band_lower_1 >= E_band_upper_1) ) {
    	throw( Exception(" Must have E_band_lower_1 < E_band_upper_1 (both in keV). Exiting.\n") );
    	return -1;
    }
    if ( !E_band_upper_2_set || !E_band_lower_2_set ) {
    	throw( Exception(" Need to set upper and lower bounds to energy band 2 for flux calculation. Exiting.\n") );
    	return -1;
    }
    if ( E_band_lower_2_set && E_band_upper_2_set && (E_band_lower_2 >= E_band_upper_2) ) {
    	throw( Exception(" Must have E_band_lower_2 < E_band_upper_2 (both in keV). Exiting.\n") );
    	return -1;
    }
    if ( numbins > MAX_NUMBINS || numbins <= 0 ) {
    	throw( Exception(" Illegal number of phase bins. Must be between 1 and MAX_NUMBINS, inclusive. Exiting.\n") );
    	return -1;
    }
    
    
    /*******************************************************/
    /* PRINT OUT INFORMATION ABOUT THE MODEL TO THE SCREEN */
    /*******************************************************/
        
    /*std::cout << "\nSpot: Stellar Model Parameters" << std::endl;
    std::cout << "Mass = " << mass << " Msun" << std::endl;
    std::cout << "Rspot = " << rspot << " km" << std::endl;
    std::cout << "2GM/Rc^2 = " << 2.0*Units::G*mass*Units::MSUN/(rspot*1e+5*Units::C*Units::C) << std::endl;
    std::cout << "Spin frequency = " << omega << " Hz" << std::endl;
    std::cout << "Theta = " << theta_1 << " degrees" << std::endl;
    std::cout << "Inclination = " << incl_1 << " degrees" << std::endl;
    std::cout << "Spot radius = " << rho << " degrees" << std::endl;
    std::cout << "Time shift = " << ts << std::endl;
    std::cout << "Distance = " << distance << " meters" << std::endl;
    if ( !T_mesh_in ) std::cout << "Spot temperature = " << spot_temperature << " keV" << std::endl;
    else std::cout << "Temperature mesh: " << T_mesh_file << std::endl;
    if ( two_spots ) std::cout << "Calculating flux from two antipodal hot spots." << std::endl;
    else std::cout << "Calculating flux from one hot spot." << std::endl;
 	*/
 	
 	/*****************************************************/
    /* UNIT CONVERSIONS -- MAKE EVERYTHING DIMENSIONLESS */
    /*****************************************************/
    
    incl_1 *= (Units::PI / 180.0);  // radians
    if ( only_second_spot ) incl_1 = Units::PI - incl_1; // for doing just the 2nd hot spot
    theta_1 *= (Units::PI / 180.0); // radians
    theta_2 = theta_1; // radians
    rho *= (Units::PI / 180.0);  // radians
    mu_1 = cos( theta_1 );
    mu_2 = mu_1; 
    mass = Units::cgs_to_nounits( mass*Units::MSUN, Units::MASS );
    rspot = Units::cgs_to_nounits( rspot*1.0e5, Units::LENGTH );
    omega = Units::cgs_to_nounits( 2.0*Units::PI*omega, Units::INVTIME );
    distance = Units::cgs_to_nounits( distance*100, Units::LENGTH );
	//std::cout << "Mass = " << mass << ",  Radius = " << rspot << ",  Freq = " << omega << std::endl;
	
    /**********************************/
    /* PASS VALUES INTO THE STRUCTURE */
    /**********************************/    	
    
    curve.para.mass = mass;
    curve.para.omega = omega;
    curve.para.radius = rspot;
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
    curve.para.rsc = r_sc;
    curve.para.Isc = I_sc;

    numphi = numtheta; // code currently only handles a square mesh over the hotspot
  
    curve.flags.ignore_time_delays = ignore_time_delays;
    curve.flags.spectral_model = spectral_model;
    curve.flags.beaming_model = beaming_model;
	//curve.flags.shift_t = ts; // ts has previously been pre-defined as 0 and 0.5/numbins

    
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
 		for (unsigned int y(0); y < NCURVES; y++) {
 			obsdata.t = new double[numbins];
 			obsdata.f[y] = new double[numbins];
 			obsdata.err[y] = new double[numbins];
 		}
 
		/****************************************/
 	  	/* READING IN FLUXES FROM THE DATA FILE */
    	/****************************************/
    	
    	while ( data.getline(line,265) ) {
        	i = numLines;
        	//std::cout << "Line = " << line << std::endl;
        	//std::cout << "i = " << i << std::endl;
			double get_f1;
			double get_err1;
			double get_f2;
			double get_err2;
			
        	sscanf( line, "%lf %lf %lf %lf %lf", &obsdata.t[i], &get_f1, &get_err1, &get_f2, &get_err2 );
      		//std::cout << "Obsdata.f[p=1][i="<<i<<"] = " << obsdata.f[1][i] << std::endl;
			obsdata.f[1][i] = get_f1;
			obsdata.err[1][i] = get_err1;
			obsdata.f[2][i] = get_f2;
			obsdata.err[2][i] = get_err2;
			//std::cout << "obsdata.f1 = " << obsdata.f[1][i] << std::endl;
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
	}
	
	//std::cout << "ts = " << ts << std::endl;
	
  	/***************************/
	/* START SETTING THINGS UP */
	/***************************/ 

    // Calculate the Equatorial Radius of the star.
    req = calcreq( omega, mass, theta_1, rspot );  // implementation of MLCB11

    /*********************************************************************************/
	/* Set up model describing the shape of the NS; oblate, funky quark, & spherical */
	/*********************************************************************************/
	
    OblModelBase* model;
    if ( NS_model == 1 ) { // Oblate Neutron Hybrid Quark Star model
        // Default model for oblateness
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
        // Use a spherical star
        model = new SphericalOblModel( rspot );
    }
    else {
        throw(Exception("\nInvalid NS_model parameter. Exiting.\n"));
        return -1;
    }

    // defltoa is a structure that "points" to routines in the file "OblDeflectionTOA.cpp"
    // used to compute deflection angles and times of arrivals 
    // defltoa is the deflection time of arrival
    OblDeflectionTOA* defltoa = new OblDeflectionTOA(model, mass); // defltoa is a pointer (of type OblDeclectionTOA) to a group of functions

    // Values we need in some of the formulas.
    cosgamma = model->cos_gamma(mu_1);   // model is pointing to the function cos_gamma
    curve.para.cosgamma = cosgamma;


    /**********************************************************/
	/* Compute maximum deflection for purely outgoing photons */
	/**********************************************************/
	
    double  b_mid;  // the value of b, the impact parameter, at 90% of b_max
    curve.defl.b_max =  defltoa->bmax_outgoing(rspot); // telling us the largest value of b
    curve.defl.psi_max = defltoa->psi_max_outgoing(curve.defl.b_max,rspot,&curve.problem); // telling us the largest value of psi

    /********************************************************************/
	/* COMPUTE b VS psi LOOKUP TABLE, GOOD FOR THE SPECIFIED M/R AND mu */
	/********************************************************************/
	
    b_mid = curve.defl.b_max * 0.9; // since we want to split up the table between coarse and fine spacing. 0 - 90% is large spacing, 90% - 100% is small spacing. b_mid is this value at 90%.
    curve.defl.b_psi[0] = 0.0; // definitions of the actual look-up table for b when psi = 0
    curve.defl.psi_b[0] = 0.0; // definitions of the actual look-up table for psi when b = 0
    
    for ( unsigned int i(1); i < NN+1; i++ ) { /* compute table of b vs psi points */
        curve.defl.b_psi[i] = b_mid * i / (NN * 1.0);
        curve.defl.psi_b[i] = defltoa->psi_outgoing(curve.defl.b_psi[i], rspot, curve.defl.b_max, curve.defl.psi_max, &curve.problem); // calculates the integral
    }
    
    // For arcane reasons, the table is not evenly spaced.
    for ( unsigned int i(NN+1); i < 3*NN; i++ ) { /* compute table of b vs psi points */
        curve.defl.b_psi[i] = b_mid + (curve.defl.b_max - b_mid) / 2.0 * (i - NN) / (NN * 1.0); // spacing for the part where the points are closer together
        curve.defl.psi_b[i] = defltoa->psi_outgoing(curve.defl.b_psi[i], rspot, curve.defl.b_max, curve.defl.psi_max, &curve.problem); // referenced same as above
    }
    
    curve.defl.b_psi[3*NN] = curve.defl.b_max;   // maximums
    curve.defl.psi_b[3*NN] = curve.defl.psi_max;
    // Finished computing lookup table

    /****************************/
	/* Initialize time and flux */
	/****************************/
	
    for ( unsigned int i(0); i < numbins; i++ ) {
        curve.t[i] = i / (1.0 * numbins);// + ts;  // defining the time used in the lightcurves
        for ( unsigned int p(0); p < NCURVES; p++ ) {
            Flux[p][i] = 0.0;                   // initializing flux to 0
            curve.f[p][i] = 0.0;
        }
    } 
    
    /****************************************/
	/* Set up the star and calculate fluxes */
	/****************************************/
	
    
	if ( rho == 0.0 ) { // Default is an infinitesmal spot, defined to be 0 degrees in radius.
		phi_0_1 = 0.0;
		theta_0_1 = theta_1;
     	curve.para.dS = trueSurfArea = 0.0; // dS is the area of the particular bin of spot we're looking at right now
   	}   
    else {   // for a nontrivially-sized spot
     	dtheta = 2.0 * rho / (numtheta * 1.0);    // center of spot at theta, rho is angular radius of spot; starts at theta_1-rho to theta_1+rho
        theta_edge_1 = theta_1 - rho;
        trueSurfArea = 2 * Units::PI * pow(rspot,2) * (1 - cos(rho));
    }
    
    /************************************/
	/* LOCATION OF THE SPOT ON THE STAR */
	/* FOR ONE OR FIRST HOT SPOT        */
	/************************************/
    bool negative_theta(false);

    /*********************************************/
	/* SPOT IS TRIVIALLY SIZED ON GEOMETRIC POLE */
	/*********************************************/
		
	if ( theta_1 == 0 && rho == 0 ) {
		for ( unsigned int i(0); i < numbins; i++ ) {
			for ( unsigned int p(0); p < NCURVES; p++ )
		       	curve.f[p][i] = 0.0;
	    }
	    
	    // Add curves
	    //if (numtheta != 1) { // don't think this matters...
	        for (unsigned int i(0); i < numbins; i++) {
		        for (unsigned int p(0); p < NCURVES; p++) {
		        	//std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
		            Flux[p][i] += curve.f[p][i];
		            //if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
		            //if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
		        }
	        }
	    //}
	} // ending trivial spot on pole
	
	/*****************************************/
	/* SPOT IS SYMMETRIC OVER GEOMETRIC POLE */
	/*****************************************/
		
	else if ( theta_1 == 0 && rho != 0 ) {
   	 	// Looping through the mesh of the spot
		for ( unsigned int k(0); k < numtheta; k++ ) {
			theta_0_1 = theta_1 - rho + k * dtheta + 0.5 * dtheta;   // note: theta_0_1 changes depending on how we've cut up the spot
			curve.para.theta = theta_0_1;
			// don't need a phi_edge the way i'm doing the phi_0_1 calculation
			dphi = Units::PI / numphi;
			for ( unsigned int j(0); j < numphi; j++ ) {
				phi_0_1 = -Units::PI/2 + dphi * j + 0.5 * dphi;
				if ( !T_mesh_in ) {
					T_mesh[k][j] = spot_temperature;
				}
				//std::cout << "\nk = " << k << ", j = " << j << std::endl;
	            /*  std::cout << "\ntheta_0_1 = " << theta_0_1*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "dtheta = " << dtheta*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "phi_0_1 = " << phi_0_1*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "dphi = " << dphi*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "fabs(theta_0_1) = " << fabs(theta_0_1)*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "sin(fabs(theta_0_1)) = " << sin(fabs(theta_0_1)) << std::endl;
	            //std::cout << "r spot = " << rspot << std::endl;
				std::cout << "theta_edge_1 = " << theta_edge_1*180.0/Units::PI << " degrees" << std::endl;
			 	std::cout << std::endl;*/
				
				curve.para.dS = pow(rspot,2) * sin(fabs(theta_0_1)) * dtheta * dphi;
				// Need to multiply by R^2 here because of my if numtheta=1 statement, 
	            // which sets dS = true surface area
	            
	            // alternative method for computing dS: dS=truesurfarea/pow(numtheta,2);
					
	            if ( numtheta == 1 ) 
	            	curve.para.dS = trueSurfArea;
	            if ( NS_model == 1 || NS_model == 2)
	            	curve.para.dS /= curve.para.cosgamma;
				curve.para.temperature = T_mesh[k][j];
	            curve.para.phi_0 = phi_0_1;
	            curve = ComputeAngles(&curve, defltoa);  // Computing the parameters it needs to compute the light curve; defltoa has the routines for what we want to do
	        	curve = ComputeCurve(&curve);
			
				if (curve.para.temperature == 0.0) { 
	        		for (unsigned int i(0); i < numbins; i++) {
		        		for (unsigned int p(0); p < NCURVES; p++)
		        			curve.f[p][i] = 0.0;  // if temperature is 0, then there is no flux!
	        		}
	        	}
				// Add curves
	       		//if (numtheta != 1) { // don't think this matters...
	        		for (unsigned int i(0); i < numbins; i++) {
		        		for (unsigned int p(0); p < NCURVES; p++) {
		        			//std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
		         			Flux[p][i] += curve.f[p][i];
		        			//if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
		         			//if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
		         		}
	        		}
	        	//}
			}
		}
	} // ending symmetric spot over pole
		
	/*********************************************/
	/* SPOT IS ASYMMETRIC OVER GEOMETRIC POLE */
	/*********************************************/
	
	else if ( (theta_1 - rho) <= 0 ) {
	
		if ( numtheta == 1 ) {
		
			curve.para.theta = theta_1;
			curve.para.phi_0 = 0.0;
			curve.para.dS = trueSurfArea;
			if ( !T_mesh_in ) {
	           		T_mesh[0][0] = spot_temperature;
	        }
			curve.para.temperature = T_mesh[0][0];

	        curve = ComputeAngles(&curve, defltoa);  // Computing the parameters it needs to compute the light curve; defltoa has the routines for what we want to do
	        curve = ComputeCurve(&curve);
	        	
	        if ( curve.para.temperature == 0.0 ) { 
	        	for ( unsigned int i(0); i < numbins; i++ ) {
		        	for ( unsigned int p(0); p < NCURVES; p++ )
		        		curve.f[p][i] = 0.0;
	        	}
	        }
	       	// Add curves
	        for ( unsigned int i(0); i < numbins; i++ ) {
		        for ( unsigned int p(0); p < NCURVES; p++ ) {
		        	//std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
		         	Flux[p][i] += curve.f[p][i];
		          	//if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
		          	//if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
		        }
	        }  		
	    } // ending if numtheta == 1
			
		else { // if numtheta != 1
			
		// Looping through the mesh of the spot
			for (unsigned int k(0); k < numtheta; k++) { // looping through the theta divisions
		// the way we are doing this, to avoid majorly messing up the spherical trig for getting phi_0,
		// is to step our way around the edge of the spot instead of slicing up the spot. 
		// we are chopping up the circumference into numtheta segments and evaluating at the center
		// of each of those segments, on the edge of the spot, instead of chopping up the area.
		// naturally, this approximation is only moderately ok for very small spots.
				dtheta = 2 * theta_1 / static_cast <double> (numtheta);
				theta_0_1 = theta_1 - rho + 0.5 * dtheta + k * dtheta;   // note: theta_0_1 changes depending on how we've cut up the spot
				if (theta_0_1 < 0.0) {
					negative_theta = true;
					theta_0_1 = fabs(theta_0_1);
				}
				if (theta_0_1 == 0.0) 
					theta_0_1 = 0.0 + DBL_EPSILON;
				curve.para.theta = theta_0_1;   // passing this value into curve theta, so that we can use it in another routine; somewhat confusing names
	        	double cos_phi_edge = (cos(rho) - cos(theta_1)*cos(theta_0_1))/(sin(theta_1)*sin(theta_0_1));
	        	if ( cos_phi_edge > 1.0 ) { 
	    			cos_phi_edge = 1.0;
	    			std::cout << "cos_phi_edge was greater than 1.0. Setting cos_phi_edge = 1.0" << std::endl;
	    		}

	    		if ( fabs( sin(theta_1) * sin(theta_0_1) ) > 0.0) { // checking for a divide by 0
	           		phi_edge_1 = acos( cos_phi_edge );   // value of phi (a.k.a. azimuth projected onto equatorial plane) at the edge of the circular spot at some latitude theta_0_1
					dphi = 2.0 * phi_edge_1 / ( numphi * 1.0 );
					phi_0_1 = acos(cos_phi_edge);
				}
	    		else {  // trying to divide by zero
    	    		throw( Exception(" Tried to divide by zero in calculation of phi_edge_1. Likely, theta_0_1 = 0. Exiting.") );
            		return -1;
        		}
				if ( negative_theta )
					phi_0_1 = -phi_0_1;
				
	          	//std::cout << "theta_0_1 = " << theta_0_1 * 180.0 / Units::PI << std::endl;
	          	//std::cout << "phi_0_1 = " << phi_0_1 * 180.0 / Units::PI << std::endl;
	        	if ( T_mesh_in ) {
	        		std::cout << "WARNING: code can't handle a spot asymmetric over the pole with a temperature mesh." << std::endl;
	          		spot_temperature = 2;
	         	}

	       	/*  std::cout << "\nk = " << k << ", j = " << j << std::endl;
	          	std::cout << "\ntheta_0_1 = " << theta_0_1*180.0/Units::PI << " degrees" << std::endl;
	          	std::cout << "inclination = " << incl_1*180.0/Units::PI << " degrees" << std::endl;
	          	std::cout << "dtheta = " << dtheta*180.0/Units::PI << " degrees" << std::endl;
	          	std::cout << "phi_0_1 = " << phi_0_1*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "dphi = " << dphi*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "fabs(theta_0_1) = " << fabs(theta_0_1)*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "sin(fabs(theta_0_1)) = " << sin(fabs(theta_0_1)) << std::endl;
	            //std::cout << "r spot = " << rspot << std::endl;
				std::cout << "theta_edge_1 = " << theta_edge_1*180.0/Units::PI << " degrees" << std::endl;
				std::cout << "phi_edge_1 = " << phi_edge_1*180.0/Units::PI << " degrees" << std::endl; 
			 	std::cout << std::endl;*/
			 	
			 	curve.para.dS = pow(rspot,2) * sin(fabs(theta_0_1)) * dtheta * dphi;
			 	// Need to multiply by R^2 here because of my if numtheta=1 statement, 
	            // which sets dS = true surface area
	            // alternative method for computing dS: dS=truesurfarea/pow(numtheta,2);

	            	
	            if ( NS_model == 1 || NS_model == 2 )
	            	curve.para.dS /= curve.para.cosgamma;
	            //std::cout << "dS = " << curve.para.dS << std::endl;
	            curve.para.temperature = spot_temperature;
	            curve.para.phi_0 = phi_0_1;

	            curve = ComputeAngles(&curve, defltoa);  // Computing the parameters it needs to compute the light curve; defltoa has the routines for what we want to do
	        	curve = ComputeCurve(&curve);
	        	
	        	if ( curve.para.temperature == 0.0 ) { 
	        		for ( unsigned int i(0); i < numbins; i++ ) {
		        		for ( unsigned int p(0); p < NCURVES; p++ )
		        			curve.f[p][i] = 0.0;
	        		}
	        	}
	       		// Add curves
	        	for ( unsigned int i(0); i < numbins; i++ ) {
		        	for ( unsigned int p(0); p < NCURVES; p++ ) {
		        		//std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
		         		Flux[p][i] += curve.f[p][i];
		          		//if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
		          		//if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
		        	}
	        	}
        	
	    	} // closing for loop through theta divisions
	    } // closing if numtheta != 1
	} // ending antisymmetric spot over pole

	/****************************************/
	/* SPOT DOES NOT GO OVER GEOMETRIC POLE */
	/****************************************/
		
	else {
   	 	// Looping through the mesh of the spot
    	for ( unsigned int k(0); k < numtheta; k++ ) { // looping through the theta divisions
			theta_0_1 = theta_1 - rho + 0.5 * dtheta + k * dtheta;   // note: theta_0_1 changes depending on how we've cut up the spot
			curve.para.theta = theta_0_1;   // passing this value into curve theta, so that we can use it in another routine; somewhat confusing names
        
	    	double cos_phi_edge = (cos(rho) - cos(theta_1)*cos(theta_0_1))/(sin(theta_1)*sin(theta_0_1));
	    	if ( cos_phi_edge > 1.0 || cos_phi_edge < -1.0 ) 
	    		cos_phi_edge = 1.0;
        
	    	if ( fabs( sin(theta_1) * sin(theta_0_1) ) > 0.0 ) { // checking for a divide by 0
	        	phi_edge_1 = acos( cos_phi_edge );   // value of phi (a.k.a. azimuth projected onto equatorial plane) at the edge of the circular spot at some latitude theta_0_1
		    	dphi = 2.0 * phi_edge_1 / ( numphi * 1.0 );
			}
	     	else {  // trying to divide by zero
    	    	throw( Exception(" Tried to divide by zero in calculation of phi_edge_1. Likely, theta_0_1 = 0. Exiting.") );
    	    	return -1;
    		}

    		for ( unsigned int j(0); j < numphi; j++ ) {   // looping through the phi divisions
	          	phi_0_1 = -phi_edge_1 + 0.5 * dphi + j * dphi;   // phi_0_1 is at the center of the piece of the spot that we're looking at; defined as distance from the y-axis as projected onto the equator
	          	if ( only_second_spot ) 
	          		phi_0_1 = Units::PI - phi_edge_1 + 0.5 * dphi + j * dphi;
	            if ( !T_mesh_in ) {
	        		T_mesh[k][j] = spot_temperature;
	        	}

	    	 /* std::cout << "\nk = " << k << ", j = " << j << std::endl;
	            std::cout << "\ntheta_0_1 = " << theta_0_1*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "inclination = " << incl_1*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "dtheta = " << dtheta*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "phi_0_1 = " << phi_0_1*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "dphi = " << dphi*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "fabs(theta_0_1) = " << fabs(theta_0_1)*180.0/Units::PI << " degrees" << std::endl;
	            std::cout << "sin(fabs(theta_0_1)) = " << sin(fabs(theta_0_1)) << std::endl;
	            //std::cout << "r spot = " << rspot << std::endl;
				std::cout << "theta_edge_1 = " << theta_edge_1*180.0/Units::PI << " degrees" << std::endl;
				std::cout << "phi_edge_1 = " << phi_edge_1*180.0/Units::PI << " degrees" << std::endl; 
			 	std::cout << std::endl;*/
           
	            curve.para.dS = rspot * rspot * sin(fabs(theta_0_1)) * dtheta * dphi;
	        	// Need to multiply by R^2 here because of my if numtheta=1 statement, 
	        	// which sets dS = true surface area
	        	
	        	// alternative method for computing dS: dS=truesurfarea/pow(numtheta,2);

	        	if( numtheta == 1 ) 
	          		curve.para.dS = trueSurfArea;
	            	
	         	if ( NS_model == 1 || NS_model == 2 )
	        		curve.para.dS /= curve.para.cosgamma;
	        	//std::cout << "dS = " << curve.para.dS << std::endl;
	        	curve.para.temperature = T_mesh[k][j];
	           	curve.para.phi_0 = phi_0_1;
	           	curve = ComputeAngles(&curve, defltoa);  // Computing the parameters it needs to compute the light curve; defltoa has the routines for what we want to do
	    		curve = ComputeCurve(&curve);            // Compute Light Curve, for each separate mesh bit
        	
	    		/*for (unsigned int i(0); i < numbins; i++) {
		       		for (unsigned int p(0); p < NCURVES; p++)
		        		if( std::isnan(curve.f[p][i]) || curve.f[p][i] == 0) std::cout << "NAN or 0 HERE" << std::endl;
	    		}*/
	    		if ( curve.para.temperature == 0.0 ) { //if the temperature of the spot is 0, then there is no x-ray emission
	    			for ( unsigned int i(0); i < numbins; i++ ) {
		      			for ( unsigned int p(0); p < NCURVES; p++ )
		          			curve.f[p][i] = 0.0;
	    			}
	    		}
	       		// Add curves
	       		//if (numtheta != 1) { // don't think this matters...
	        		for ( unsigned int i(0); i < numbins; i++ ) {
		    			for ( unsigned int p(0); p < NCURVES; p++ ) {
		    				//std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
		            		Flux[p][i] += curve.f[p][i];
		            		//if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
		            		//if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
		            	}
	        		}	
	        	//}
        	
	        } // closing for loop through phi divisions
	    } // closing for loop through theta divisions
    
    } // closing spot does not go over the pole or is antisymmetric over the pole
    
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
		 		for ( unsigned int p(0); p < NCURVES; p++ )
		     		curve.f[p][i] = 0.0;
	    	}
	    
			// Add curves
	       	//if (numtheta != 1) { // don't think this matters...
	      		for (unsigned int i(0); i < numbins; i++) {
		      		for (unsigned int p(0); p < NCURVES; p++) {
		      			//std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
		          		Flux[p][i] += curve.f[p][i];
		          		//if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
		          		//if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
		          	}
	      		}
	     	//}
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
					//std::cout << "\nk = " << k << ", j = " << j << std::endl;
	          	/*  std::cout << "\ntheta_0_2 = " << theta_0_2*180.0/Units::PI << " degrees" << std::endl;
	            	std::cout << "dtheta = " << dtheta*180.0/Units::PI << " degrees" << std::endl;
	            	std::cout << "phi_0_2 = " << phi_0_2*180.0/Units::PI << " degrees" << std::endl;
	            	std::cout << "dphi = " << dphi*180.0/Units::PI << " degrees" << std::endl;
	            	std::cout << "fabs(theta_0_2) = " << fabs(theta_0_2)*180.0/Units::PI << " degrees" << std::endl;
	            	std::cout << "sin(fabs(theta_0_2)) = " << sin(fabs(theta_0_2)) << std::endl;
	            	//std::cout << "r spot = " << rspot << std::endl;
					std::cout << "theta_edge_2 = " << theta_edge_2*180.0/Units::PI << " degrees" << std::endl;
			 		std::cout << std::endl;*/
				
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
		    				for (unsigned int p(0); p < NCURVES; p++)
	            				curve.f[p][i] = 0.0;  // if temperature is 0, then there is no flux!
        				}
	        		}
					// Add curves
	       			//if (numtheta != 1) { // don't think this matters...
	     				for (unsigned int i(0); i < numbins; i++) {
	        				for (unsigned int p(0); p < NCURVES; p++) {
	        					//std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
		         				Flux[p][i] += curve.f[p][i];
		        				//if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
	            				//if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
	            			}
        				}
        			//}
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
	 	    	//std::cout << cos_phi_edge << std::endl;
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
	         			//std::cout << "Only in here for full?" << std::endl;
	          		}

	     			/*  std::cout << "\nk = " << k << ", j = " << j << std::endl;
	         		std::cout << "\ntheta_0_2 = " << theta_0_2*180.0/Units::PI << " degrees" << std::endl;
	         		std::cout << "inclination = " << incl_1*180.0/Units::PI << " degrees" << std::endl;
	        		std::cout << "dtheta = " << dtheta*180.0/Units::PI << " degrees" << std::endl;
	        		std::cout << "phi_0_2 = " << phi_0_2*180.0/Units::PI << " degrees" << std::endl;
            		std::cout << "dphi = " << dphi*180.0/Units::PI << " degrees" << std::endl;
            		std::cout << "fabs(theta_0_2) = " << fabs(theta_0_2)*180.0/Units::PI << " degrees" << std::endl;
            		std::cout << "sin(fabs(theta_0_2)) = " << sin(fabs(theta_0_2)) << std::endl;
            		//std::cout << "r spot = " << rspot << std::endl;
					std::cout << "theta_edge_2 = " << theta_edge_2*180.0/Units::PI << " degrees" << std::endl;
					std::cout << "phi_edge_2 = " << phi_edge_2*180.0/Units::PI << " degrees" << std::endl; 
		 			std::cout << std::endl;*/
			 	
		 			curve.para.dS = pow(rspot,2) * sin(fabs(theta_0_2)) * dtheta * dphi;
		 			// Need to multiply by R^2 here because of my if numtheta=1 statement, 
            		// which sets dS = true surface area
            	
            		if( numtheta == 1 ) 
            			curve.para.dS = trueSurfArea;
            		if ( NS_model == 1 || NS_model == 2 )
            			curve.para.dS /= curve.para.cosgamma;
            		//std::cout << "dS = " << curve.para.dS << std::endl;
            		curve.para.temperature = T_mesh[k][j];
            		curve.para.phi_0 = phi_0_2;
	            	
	            	curve = ComputeAngles(&curve, defltoa);  // Computing the parameters it needs to compute the light curve; defltoa has the routines for what we want to do
        			curve = ComputeCurve(&curve);
	        	
	    			if ( curve.para.temperature == 0.0 ) { 
	    				for ( unsigned int i(0); i < numbins; i++ ) {
	        				for ( unsigned int p(0); p < NCURVES; p++ )
	            				curve.f[p][i] = 0.0;
        				}
        			}
       				// Add curves
       				//if (numtheta != 1) {
        				for ( unsigned int i(0); i < numbins; i++ ) {
	        				for ( unsigned int p(0); p < NCURVES; p++ ) {
	        					//std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
	            				Flux[p][i] += curve.f[p][i];
	            				//if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
	            				//if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
	            			}
        				}
        			//}
       	
       			} // closing for loop through phi divisions
	    	} // closing for loop through theta divisions
		} // ending antisymmetric spot over pole

		/***********************************************************/
		/* SECOND HOT SPOT -- SPOT DOES NOT GO OVER GEOMETRIC POLE */
		/***********************************************************/
    
		else {
       		// Looping through the mesh of the second spot
   			for ( unsigned int k(0); k < numtheta; k++ ) { // looping through the theta divisions
				theta_0_2 = theta_2 - rho + 0.5 * dtheta + k * dtheta;   // note: theta_0_2 changes depending on how we've cut up the spot
       			//std::cout << "Theta 2 = " << theta_0_2 * 180 / Units::PI << std::endl; // so that it prints to the screen in degrees for us!
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
           				//std::cout << "Only in here for full?" << std::endl;
           			}
           			/*std::cout << "\nk = " << k << ", j = " << j << std::endl;
           			std::cout << "inclination = " << incl_2 * 180 / Units::PI << std::endl;
           			std::cout << "theta_0_2 = " << theta_0_2*180/Units::PI << std::endl;
           			std::cout << "dtheta = " << dtheta*180/Units::PI << std::endl;
					std::cout << "phi_0_2 = " << phi_0_2*180/Units::PI << std::endl;
           			std::cout << "dphi = " << dphi*180/Units::PI << std::endl;
           			//std::cout << "fabs(theta_0_2) = " << fabs(theta_0_2) *180./Units::PI << std::endl;
           			//std::cout << "sin(fabs(theta_0_2)) = " << sin(fabs(theta_0_2)) << std::endl;
           			//std::cout << "r spot = " << rspot << std::endl;
					std::cout << "theta_edge_2 = " << theta_edge_2*180/Units::PI << std::endl;
					std::cout << "phi_edge_2 = " << phi_edge_2*180/Units::PI << std::endl;
					std::cout << std::endl;*/
           
           			curve.para.temperature = T_mesh[k][j];
           			curve.para.phi_0 = phi_0_2;
           			curve.para.dS = pow(rspot,2) * sin(fabs(theta_0_2)) * dtheta * dphi; // assigning partial dS here
           			if( numtheta == 1 ) 
           				curve.para.dS = trueSurfArea;
           			if ( NS_model == 1 || NS_model == 2 )
            			curve.para.dS /= curve.para.cosgamma;
           			//std::cout << "dS = " << curve.para.dS << std::endl;
           			
					curve = ComputeAngles(&curve, defltoa);  // Computing the parameters it needs to compute the light curve; defltoa has the routines for what we want to do
       				curve = ComputeCurve(&curve);            // Compute Light Curve, for each separate mesh bit

       				/*for (unsigned int i(0); i < numbins; i++) {
        				for (unsigned int p(0); p < NCURVES; p++)
            				if( (std::isnan(curve.f[p][i]) || curve.f[p][i] == 0) && p == 1) std::cout << "NAN or 0 HERE, i = " << i << std::endl;
       				}*/
        	
        			if ( curve.para.temperature == 0.0 ) { 
        				for ( unsigned int i(0); i < numbins; i++ ) {
	       					for ( unsigned int p(0); p < NCURVES; p++ )
	           					curve.f[p][i] += 0.0;
        				}
        			}
       				// Add curves
       				//if (numtheta != 1) { // don't think this matters...
        				for ( unsigned int i(0); i < numbins; i++ ) {
	       					for ( unsigned int p(0); p < NCURVES; p++ ) {
	       						//std::cout << Flux[p][i] << ", " << curve.f[p][i] << std::endl;
	           					Flux[p][i] += curve.f[p][i];
	           					//if (p == NCURVES-1) std::cout << "Time = " << curve.t[i] << ", i = " << i << ", Flux = " << Flux[p][i] << std::endl;
	           					//if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is NAN or 0 at p="<<p<<", i="<<i << std::endl;
	           				}
        				}
        			//}
        
        		} // closing for loop through phi divisions
    		} // closing for loop through theta divisions
		} // closing spot doesn't go over geometric pole
    		
	} // closing if two spots
	
	//std::cout << "SPOT dOmega = " << curve.dOmega_s[0] << std::endl;
	//std::cout << "spot1 curve.f[3][28] = " << curve.f[3][28] << std::endl;
	
	// printing the temperature mesh
	 /* std::cout << "Temperature Mesh of Spot:" << std::endl;
    	for (unsigned int i(0); i < numtheta; i++) {
    		for (unsigned int j(0); j < numphi; j++) {
    	    	std::cout << T_mesh[i][j] << "  ";
    		}
    		std::cout << std::endl;
    	}*/
     
     
     
    // read in from prebinned, put that into Flux proper/associated column, copy its output into a different file 
    // plot this output with outputFerret flux for same case 
     
 /*    
     
     	std::ifstream prebinned_in;
		prebinned_in.open( "run_data/05-Mar-2013_prebinned_flux_high.txt" );  // opening the

		char line[265]; // line of the data file being read in
		unsigned int numLines(0), z(0);

 		if ( prebinned_in.fail() || prebinned_in.bad() || !prebinned_in ) {
   			throw( Exception("Couldn't open data file."));
   			return -1;
 		}
 		
 		double prebinned_time[MAX_NUMBINS];
		double prebinned_Flux[MAX_NUMBINS]; 	
*/
  /*  	
    	while ( prebinned_in.getline(line,265) ) {
        	z = numLines;			
        	sscanf( line, "%lf %lf", &prebinned_time[z], &prebinned_Flux[z] );
        	numLines++;
        } 

        if ( numLines != numbins ) {
        	throw (Exception( "Numbins indicated in command-line not equal to numbins in data file."));
        	std::cout << "Command-line numbins = " << numbins <<", data file numbins = " << numLines << std::endl;
        	return -1;
        }

 		prebinned_in.close();
 		
 		
 		for (unsigned int k(0); k < numbins; k++) {
 			Flux[3][k] = prebinned_Flux[k];
 		}
 		
 		
     
     
     */
     
     
     
     
     
     
 	/*******************************/
	/* NORMALIZING THE FLUXES TO 1 */
	/*******************************/
	
	// Initializing normcurve to 0
	for ( unsigned int i(0); i < numbins; i++ ) {
     	for ( unsigned int p(0); p < NCURVES; p++ ) {
     		normcurve.f[p][i] = 0.0;
	           //if( std::isnan(normcurve.f[p][i]) || normcurve.f[p][i] == 0) std::cout << "normcurve.f is NAN or 0 at p="<<p<<", i="<<i << std::endl;
     	}
    }  
    
    // Normalizing the flux to 1
 	if ( normalize_flux ) {
    	normcurve = Normalize( Flux, numbins );
     	// fun little way around declaring the Chi.cpp method Normalize as returning a matrix! so that we can still print out Flux
     	for ( unsigned int i(0); i < numbins; i++ ) {
     		for ( unsigned int p(0); p < NCURVES; p++ ) {
     			Flux[p][i] = normcurve.f[p][i];
     			curve.f[p][i] = Flux[p][i]; // loading flux into curve structure, since chisquared routine gets passed curve not Flux

	            	//if( std::isnan(normcurve.f[p][i]) || normcurve.f[p][i] == 0) std::cout << "normcurve.f is NAN or 0 at p="<<p<<", i="<<i << std::endl;
     		}
     	}  
    }
     
    /*
    std::ofstream test2;
    char test2_file[265] = "./testout_2.txt";
    test2.open(test2_file, std::ios_base::trunc);
    if( test2.bad() || test2.fail() ){
    	throw (Exception("Testout 2 failed."));
    	return -1;
    }
    test2 << "# time \t Flux" << std::endl;
    for ( unsigned int i(0); i < numbins; i++ ) {
	    for ( unsigned int p(0); p < NCURVES; p++ ) {
	        test2 << curve.t[i] << "\t" << Flux[p][i] << std::endl;
	    }
    }
    test2.close();
    */    
    
    /************************************************************/
	/* If data file is set, calculate chi^2 fit with simulation */
	/************************************************************/
	
    if ( datafile_is_set ) {
    	chisquared = ChiSquare ( &obsdata, &curve );
    	//std::cout << "X^2 = " << chisquared << std::endl;
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
    for ( unsigned int j(0); j < NCURVES; j++ ) {
    	sumMaxFlux += curve.maxFlux[j];
    	sumMinFlux += curve.minFlux[j];
        avgPulseFraction += curve.pulseFraction[j];
        //std::cout << "Maximum = " << curve.maxFlux[j] << std::endl;
        //std::cout << "Minimum = " << curve.minFlux[j] << std::endl; 
        //if (j == 2 || j == 3) std::cout << "Pulse Fraction ["<<j<<"] = " << curve.pulseFraction[j] << std::endl; 
    }
	overallPulseFraction = (sumMaxFlux - sumMinFlux) / (sumMaxFlux + sumMinFlux);
    avgPulseFraction /= (NCURVES);
    //std::cout << "Averaged Pulse Fraction = " << avgPulseFraction << std::endl;
    //std::cout << "Overall Pulse Fraction = " << overallPulseFraction << std::endl;
	
	double U(0.0), Q(0.0), A(0.0);   // as defined in PG19 and PG20
    U = (1 - 2 * mass / rspot) * sin(incl_1) * sin(theta_1); // PG20
    Q = 2 * mass / rspot + (1 - 2 * mass / rspot) * cos(incl_1) * cos(theta_1); // PG20
    A = U / Q; // PG19
    B = rspot * sin(incl_1) * sin(theta_1) * omega / ( sqrt( 1 - ((2*mass) / rspot) ) ); // param_degen/equations.pdf 2
	//std::cout << "B ~ v/c ~ " << B << std::endl;
    //std::cout << "U = " << U << std::endl;
    //std::cout << "Q = " << Q << std::endl;
    //std::cout << "A = U/Q = " << A << "\n" << std::endl;
    
	
    std::ofstream spot_run_out;
    spot_run_out.open(spot_run_out_file, std::ios_base::app);
    if ( spot_run_out.bad() || spot_run_out.fail() ) {
    	std::cerr << "Spot_run_out had problems. Exiting." << std::endl;
    	return -1;
    }
    spot_run_out << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN << "\t"
    		 	 << Units::nounits_to_cgs(rspot, Units::LENGTH )*1.0e-5 << "\t"
    		 	 << mass / rspot << "\t"
    			 << incl_1 * 180.0/Units::PI << "\t"
    			 << theta_1 * 180.0/Units::PI << "\t"
    			 << ts << "\t"
    			 << A << "\t"
    			 << B << "\t"
    			 << curve.pulseFraction[2] << "\t"
    			 << curve.pulseFraction[3] << "\t"
    			 << chisquared << std::endl;
    
    spot_run_out.close();
    
    if (pd_neg_soln) {
		std::ofstream spot_run_out_neg;
    	spot_run_out_neg.open(spot_run_out_file_neg, std::ios_base::app);
    	if ( spot_run_out_neg.bad() || spot_run_out_neg.fail() ) {
    		std::cerr << "Spot_run_out_neg had problems. Exiting." << std::endl;
    		return -1;
    	}
    	spot_run_out_neg << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN << "\t"
    		 	 		<< Units::nounits_to_cgs(rspot, Units::LENGTH )*1.0e-5 << "\t"
    		 	 		<< mass / rspot << "\t"
    			 		<< incl_1 * 180.0/Units::PI << "\t"
    			 		<< theta_1 * 180.0/Units::PI << "\t"
    			 		<< ts << "\t"
    			 		<< A << "\t"
    			 		<< B << "\t"
    			 		<< curve.pulseFraction[2] << "\t"
    			 		<< curve.pulseFraction[3] << "\t"
    			 		<< chisquared << std::endl;
    
    	spot_run_out_neg.close();
    }
    else {
    	std::ofstream spot_run_out_pos;
    	spot_run_out_pos.open(spot_run_out_file_pos, std::ios_base::app);
    	if ( spot_run_out_pos.bad() || spot_run_out_pos.fail() ) {
    		std::cerr << "Spot_run_out_pos had problems. Exiting." << std::endl;
    		return -1;
    	}
    	spot_run_out_pos << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN << "\t"
    		 	 		<< Units::nounits_to_cgs(rspot, Units::LENGTH )*1.0e-5 << "\t"
    		 	 		<< mass / rspot << "\t"
    			 		<< incl_1 * 180.0/Units::PI << "\t"
    			 		<< theta_1 * 180.0/Units::PI << "\t"
    			 		<< ts << "\t"
    			 		<< A << "\t"
    			 		<< B << "\t"
    			 		<< curve.pulseFraction[2] << "\t"
    			 		<< curve.pulseFraction[3] << "\t"
    			 		<< chisquared << std::endl;
    
    	spot_run_out_pos.close();
    }


    /********************************************/
    /* WRITING THE SIMULATION TO AN OUTPUT FILE */
    /********************************************/ 
    
    
    // MKDIR TO MAKE THE OUTPUT FILE??
	
    out.open(out_file, std::ios_base::trunc);
    if ( out.bad() || out.fail() ) {
        std::cerr << "Couldn't open output file. Exiting." << std::endl;
        return -1;
    }

    if ( rho == 0.0 ) rho = Units::PI/180.0;

    out << "# Photon Flux for hotspot on a NS. \n"
        << "# R = " << Units::nounits_to_cgs(rspot, Units::LENGTH )*1.0e-5 << " km; "
        << "# M = " << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN << " Msun; "
        << "# Spin = " << Units::nounits_to_cgs(omega, Units::INVTIME)/(2.0*Units::PI) << " Hz \n"
        << "# Gravitational Redshift 1+z = " << 1.0/sqrt(1.0 - 2*mass/rspot) << "\n"
        << "# Inclination Angle = " << incl_1 * 180.0/Units::PI << " degrees \n"
        << "# Spot 1 at Angle = " << theta_1 * 180.0/Units::PI << " degrees \n"
        //<< "# Spot 2 at Angle = " << theta_2 * 180.0/Units::PI << " degrees \n"
        << "# Angular Radius of Spot = " << rho * 180.0/Units::PI << " degrees \n"
        << "# Number of spot bins = " << numtheta << " \n"
        << "# Distance to NS = " << Units::nounits_to_cgs(distance, Units::LENGTH)*.01 << " m \n" 
        << "# Phase shift or time lag = " << ts << " \n"
        << "# Pulse fractions: bolo = " << curve.pulseFraction[0] <<", mono = " << curve.pulseFraction[1] << ", \n"
        << "#                  low E band = " << curve.pulseFraction[2] << ", high E band = " << curve.pulseFraction[3] << " \n"
        << "# B = " << B << " \n"
        << "# U = " << U << " \n"
        << "# Q = " << Q << " \n"
        << "# A = U/Q = " << A 
        << std::endl;
        
    if (datafile_is_set)
    	out << "# Data file " << data_file << ", chisquared = " << chisquared << std::endl;

    if ( T_mesh_in )
    	out << "# Temperature mesh input: "<< T_mesh_file << std::endl;
    else
    	out << "# Spot Temperature, (star's frame) kT = " << spot_temperature << " keV " << std::endl;
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
    	<< "# Column 1: phase bins (0 to 1)\n"
        << "# Column 2: Bolometric Number flux, photons/(cm^2 s) \n"
        << "# Column 3: Monochromatic Number flux (photons/(cm^2 s keV) measured at energies (at infinity) of 2 keV\n" 
        << "# Column 4: Number flux (photons/(cm^2 s)) in the energy band " << E_band_lower_1 << " keV to " << E_band_upper_1 << " keV \n"
    	<< "# Column 5: Number flux (photons/(cm^2 s)) in the energy band " << E_band_lower_2 << " keV to " << E_band_upper_2 << " keV \n"
        << "#"
        << std::endl;

    for ( unsigned int i(0); i < numbins; i++ ) {
        out << curve.t[i]<< "\t" ;
        for ( unsigned int p(0); p < NCURVES; p++ ) { // bolometric + 1 monochromatic + 2 energy bands
        	if( std::isnan(Flux[p][i]) || Flux[p][i] == 0) std::cout << "Flux is "<< Flux[p][i]<<" at p="<<p<<", i="<<i << std::endl;
            out << Flux[p][i] << "\t" ;
            //std::cout << " i = " << i << ", p = "<< p << ", Flux = " << Flux[p][i] << std::endl;
        }
        //std::cout << std::endl;
        out << std::endl;
    }
    out.close();
    
    /***************************************************************************/
    /* WRITING PARAMETERS TO AN OUTPUT FILE, FOR KEEPING TRACK OF LOTS OF RUNS */
    /***************************************************************************/

    param_out.open(param_out_file, std::ios_base::app | std::ios_base::out);
    if( param_out.bad() || param_out.fail() ){
    	throw (Exception("Testout failed."));
    	return -1;
    }
    param_out << "Comparing with input file "<<data_file<<" \n"
    		  << "Chi^2 = " << chisquared << " \n"
    		  << "R = " << Units::nounits_to_cgs(rspot, Units::LENGTH )*1.0e-5 << " km;  "
        	  << "M = " << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN << " Msun;  "
        	  << "1+z = " << 1.0/sqrt(1.0 - 2*mass/rspot) << ";  "
        	  << "ts = " << ts << " \n"
        	  << "Incl = " << incl_1 * 180.0/Units::PI << " degrees;  "
        	  << "Theta = " << theta_1 * 180.0/Units::PI << " degrees;  "
              << "Rho = " << rho * 180.0/Units::PI << " degrees \n"
        	  << "Spin = " << Units::nounits_to_cgs(omega, Units::INVTIME)/(2.0*Units::PI) << " Hz \n"
        	  << "Numtheta = " << numtheta << " \n"
        	  << "Distance to NS = " << Units::nounits_to_cgs(distance, Units::LENGTH)*.01 << " m \n"
        	  << "B = " << B << " \n"
              << "U = " << U << " \n"
              << "Q = " << Q << " \n"
              << "A = U/Q = " << A
        	  << std::endl;
        	  
    if ( T_mesh_in )
    	param_out << "Temperature mesh input: "<< T_mesh_file << std::endl;
    else
    	param_out << "Spot Temperature = " << spot_temperature << " keV " << std::endl;
	if ( NS_model == 1)
    	param_out << "Oblate NS model;  ";
    else if (NS_model == 3)
    	param_out << "Spherical NS model;  ";
    if ( beaming_model == 0 )
        param_out << "Isotropic Emission;  ";
    else
        param_out << "Graybody;  ";
    if ( normalize_flux )
    	param_out << "Flux normalized to 1.  " << std::endl;
    else
    	param_out << "Flux not normalized.  " << std::endl;
    	
    param_out << std::endl;
    
    param_out.close();
    
    /*******************************************************/
    /* WRITING TO A TEST OUTPUT FILE, WITH FAKE ERROR BARS */
    /*******************************************************/

    testout.open(testout_file, std::ios_base::trunc);
    if( testout.bad() || testout.fail() ){
    	throw (Exception("Testout failed."));
    	return -1;
    }
    testout << "# Photon Flux for hotspot on a NS. \n"
            << "# R = " << Units::nounits_to_cgs(rspot, Units::LENGTH )*1.0e-5 << " km; "
        	<< "# M = " << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN << " Msun; "
        	<< "# Spin = " << Units::nounits_to_cgs(omega, Units::INVTIME)/(2.0*Units::PI) << " Hz \n"
        	<< "# Gravitational Redshift 1+z = " << 1.0/sqrt(1.0 - 2*mass/rspot) << "\n"
        	<< "# Inclination Angle = " << incl_1 * 180.0/Units::PI << " degrees \n"
        	<< "# Spot 1 at Angle = " << theta_1 * 180.0/Units::PI << " degrees \n"
        	//<< "# Spot 2 at Angle = " << theta_2 * 180.0/Units::PI << " degrees \n"
        	<< "# Angular Radius of Spot = " << rho * 180.0/Units::PI << " degrees \n"
        	<< "# Number of spot bins = " << numtheta << " \n"
        	<< "# Distance to NS = " << Units::nounits_to_cgs(distance, Units::LENGTH)*.01 << " m \n"
        	<< "# Spot Temperature, (star's frame) kT = " << spot_temperature << " keV \n"
       		<< "# Phase shift or time lag = " << ts << " \n"
        	<< "# Avg Pulse Fraction = " << avgPulseFraction
        	<< std::endl;
        	
	if ( NS_model == 1)
    	testout << "# Oblate NS model " << std::endl;
    else if (NS_model == 3)
    	testout << "# Spherical NS model " << std::endl;
    if ( beaming_model == 0 )
        testout << "# Isotropic Emission " << std::endl;
    else
        testout << "# Limb Darkening for Gray Atmosphere (Hopf Function) " << std::endl;
    if ( normalize_flux )
    	testout << "# Flux normalized to 1 " << std::endl;
    else
    	testout << "# Flux not normalized " << std::endl;
    
    	
    testout << "#\n"
    		<< "# Column 1: Phase bins (0 to 1)\n"
        	<< "# Column 2: Normalized energy band flux, " << E_band_lower_1 << " - " << E_band_upper_1 << " keV \n" //, in photons/(cm^2 s) \n" 
        	<< "# Column 3: Normalized energy band flux, " << E_band_lower_2 << " - " << E_band_upper_2 << " keV \n" //, in photons/(cm^2 s) \n" 
        	<< "# Column 4: dOmega_s \n" 
        	<< "# Column 5: eta \n"
        	<< "#"
        	<< std::endl;
    for (unsigned int i(0); i < numbins; i++) {
    	testout << curve.t[i] << " " << Flux[NCURVES-2][i] << " " 
    	<< Flux[NCURVES-1][i] << " " << curve.dOmega_s[i] << " " << curve.eta[i]
    	<< std::endl;
    }
    testout.close();
        	

    delete defltoa;
    delete model;
    return 0;

} 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
