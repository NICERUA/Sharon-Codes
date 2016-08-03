/************************************************************************** 
*                            int.c                                        *
*        A simple integration routine.
*     
*                      gcc -lm -o ode ode.c
* 
*
*                      ode -n 2
*                           
*                                                                         *
*  Author:  Sharon Morsink	                                          *
*                                                                         *   
*  Address:                                                               *
*                                                                         *
*           Department of Physics                                         *
*           University of Alberta                                         *    
*           Edmonton, AB, Canada, T6G 2J1                                 *  
*                                                                         *
*           E-mail: morsink@phys.ualberta.ca                              *
*                   
*                                                                         *
*  Date:    Jan. 2001                                                     * 
*                                                                         * 
***************************************************************************/

/* The following standard libraries have to be included */

#include <stdio.h>
#include <string.h> 
#include <math.h>




/*************************************************************************/
/* Main program.                                                         */
/*************************************************************************/
int main(int argc,                    /* Number of command line arguments */ 
         char **argv)                 /* Command line arguments */
{

  /* Define the variables */


  FILE *input;
  FILE  *output;

  int i;

  double phase[33], low[33], lerr[33], high[33], herr[33];

  double lnorm, hnorm;

  char out_file[180]="No file specified", in_file[180]="No file specified";

  lnorm=0.0;
  hnorm=0.0;

    for ( i=1; i < argc; i++ ) {
        if ( argv[i][0] == '-' ) {  // the '-' flag lets the computer know that we're giving it information from the cmd line
            switch ( argv[i][1] ) {

	    case 'o':  // Name of output file
	      sscanf(argv[i+1], "%s", out_file);
	      break;

	    case 'i':  // Name of input file
	      sscanf(argv[i+1], "%s", in_file);
	      break;


	    }
	}
    }



 input = fopen(in_file,"r");
 output = fopen(out_file,"w");


   for(i=1;i<=(32);i++) {  
     
     fscanf(input,"%lf %lf %lf %lf %lf \n",&phase[i],&low[i],&lerr[i],&high[i],&herr[i]) ;
      
     lnorm += low[i];
     hnorm += high[i];

    }

   lnorm /= 32.0;
   hnorm /= 32.0;

   for (i=1;i<=32;i++){

     fprintf(output,"%lf %lf %lf %lf %lf \n",
	     phase[i],low[i]/lnorm,lerr[i]/lnorm,high[i]/hnorm,herr[i]/hnorm);

   }




  
}




