/* The sor06.C main program */

// $Id: sor06.C,v 1.13 2010-10-20 15:38:28 mette Exp $

#include "sor06.h"
#include "tools.h"
#include "arteries.h"
#include "junction.h"

extern "C"  void impedance_init_driver_(int *tmstps);

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>

using namespace std;

// The vessel-network is specified and initialized, and the flow and
// pressures are to be determed. All global constants must be defined
// in the header file. That is the number of dots per cm, and the number
// of vessels.
int main(int argc, char *argv[])
{
    double tstart, tend, finaltime;
    // Declare parameters: MJC 9/2019
    double f1, f2, f3, fs1, fs2, fs3;
    double alpha,beta, lrr, rm, Z0, vaso_L, vaso_R;
    int total_vessels, total_terminal, total_conn, num_pts;
    int ST_calc, ves_stiff;
    double vaso_INPUT=0.0;
    
    if (argc < 14) {
        fprintf(stdout,"Not enough entries: Only %d inputs. Exiting now.\n",argc);
        return 1;
    }
    ////// Load in fluids parameters
    f1    = atof(argv[1]);
    f2    = atof(argv[2]);
    f3    = atof(argv[3]);
    
    fs1   = atof(argv[4]);
    fs2   = atof(argv[5]);
    fs3   = atof(argv[6]);

    alpha = atof(argv[7]);
    beta  = atof(argv[8]);
    lrr   = atof(argv[9]);
    rm    = atof(argv[10]);
    Z0    = atof(argv[11]);
    vaso_L  = atof(argv[12]);
    vaso_R  = atof(argv[13]);
    
    
    int cycles = 1; // = 1; // NOTE: Can make this a parameter to pass in

    /////// Load in network parameters
    total_vessels  = 21; // Number of vessels
    total_terminal = 11; // Number of terminal vessels
    num_pts        = 40; // Spatial resolution (e.g. delta X = 1/40)
    total_conn     = total_vessels-total_terminal;
    nbrves         = total_vessels;
    
    // A parameter that can alleviate calculating the impedance for each iteration
    ST_calc = 1;
        
    
    // Workspace used by bound_right
    for(int i=0; i<4; i++) fjac_ST[i] = new double[4];
    // Workspace used by bound_bif
    for(int i=0; i<24; i++) fjac[i] = new double[24];
    
//  clock_t c1 = clock();        // Only used when timing the program.
    tstart    = 0.0;            // Starting time.
    
    // The number of vessels in the network is given when the governing array of
    // vessels is declared.
    
    impedance_init_driver_(&tmstps);
    Tube   *Arteries[nbrves];                    // Array of blood vessels.
    
    
    
    
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
////////////////// Load in network ////////////////////////////////////////
    // MJC: This code is typically for bifurcatig branches. We will have to
    // also load in the number of stenosed vessels so that the counting is
    // not off.
//    int conn_rows = (total_vessels-1+total_sten)/2;
    int conn_rows = total_vessels;//(total_vessels-1)/2;
    int connectivity_matrix[conn_rows][max_D];
    int terminal_vessels[total_terminal];
    FILE *conn;
    conn = fopen("connectivity.txt","rt");
    int parent, daughter1, daughter2, daughter3, r_in;
    int conn_id = 0;
    
    // Check to see if we have the connecitivity file
    if (conn == NULL)
    {
        fprintf(stdout,"Error: Connectivity File Does Not Exist \n");
        return 1;
    }
    while ((r_in = fscanf(conn, "%d %d %d %d", &parent, &daughter1, &daughter2, &daughter3)) != EOF)
    {
        connectivity_matrix[conn_id][0] = parent;
        connectivity_matrix[conn_id][1] = daughter1;
        connectivity_matrix[conn_id][2] = daughter2;
        connectivity_matrix[conn_id][3] = daughter3;
        //fprintf(stdout,"parent: %d d1: %d d2: %d d3: %d\n",parent,daughter1,daughter2,daughter3);
        conn_id++;
    }
    fclose(conn);
    
    
    /////////////// Load in terminal vessels
    
    FILE *terminal_file;
    terminal_file = fopen("terminal_vessels.txt","rt");
    for (int i=0; i<total_terminal; i++){
        fscanf(terminal_file, "%d", &terminal_vessels[i]);
    }
    fclose(terminal_file);
    
    
    /////////////// Load in Dimensions
    FILE *dim_file;
    dim_file = fopen("Dimensions.txt","rt");
    double dimensions_matrix[total_vessels][3];
    for (int dim_ID = 0; dim_ID < total_vessels; dim_ID++){
        fscanf(dim_file, "%6lf %6lf %6lf", &dimensions_matrix[dim_ID][0],&dimensions_matrix[dim_ID][1],&dimensions_matrix[dim_ID][2]);
    }
    fclose(dim_file);
    
//    for (int dim_ID = 0; dim_ID < total_vessels; dim_ID++){
//        fprintf(stdout,"Dimensions: %5.4f %5.4f %5.4f \n", dimensions_matrix[dim_ID][0],dimensions_matrix[dim_ID][1],dimensions_matrix[dim_ID][2]);
//    } 
    
    /* Initialization of the Arteries.
    // Definition of Class Tube:
//     Tube (double Length,
//               double topradius, double botradius,
//               Tube *LeftDaughter, Tube *RightDaughter,
//               double rmin_in, double points, int init, double K,
//               double f1, double f2, double f3, double fs1, double fs2, double fs3, double alpha_in, double beta_in, double lrr_in, double terminal_ID):
     Note: In this version, we pass the terminal ID to save impdance files from the ST*/
    /*                 FORWARD VERSION             */
    
    // MJC: For trying to pass junction condition around
    int *daughter_ptr = &connectivity_matrix[0][0];
    int term_id = total_terminal-1;
    conn_id = total_conn-1;//+total_sten;

    
///// A hard code for a single vessel
    if (total_vessels == 1){
        Arteries[0] = new Tube( dimensions_matrix[0][0], dimensions_matrix[0][1], dimensions_matrix[0][2],
                               daughter_ptr, rm, num_pts, 1,0,f1,f2,f3,fs1,fs2,fs3, alpha, beta, lrr,Z0,0,ST_calc,vaso_INPUT);
    }
    else if (total_vessels > 1){
        for (int i=total_vessels-1; i>=0; i--) {
            if (i==0) {
               fprintf(stdout,"At first vessel %d with dimensions %lf %lf %lf\n",i,dimensions_matrix[i][0],dimensions_matrix[i][1],dimensions_matrix[i][2]);
                
                Arteries[i] = new Tube( dimensions_matrix[0][0], dimensions_matrix[0][1], dimensions_matrix[0][2],
                                       daughter_ptr, 0, num_pts, 1,0, f1,f2,f3,fs1,fs2,fs3, 0, 0, 0,0, 0,ST_calc,vaso_INPUT);

            }
            else{
                if (i == terminal_vessels[term_id] && term_id>=0){
                    if (i==4 || i==8 || i==10 || i==11 || i==12){
                        vaso_INPUT = vaso_R;
                    }
                    else if (i==6 || i==15 || i==16 || i==18 || i==19 || i==20){
                        vaso_INPUT = vaso_L;
                    }
//                   fprintf(stdout,"At terminal vessel %d with dimensions %lf %lf %lf\n",i,dimensions_matrix[i][0],dimensions_matrix[i][1],dimensions_matrix[i][2]);
                    
                        Arteries[i] = new Tube( dimensions_matrix[i][0], dimensions_matrix[i][1], dimensions_matrix[i][2],
                        daughter_ptr, rm, num_pts, 0, 0,f1,f2,f3,fs1,fs2,fs3, alpha, beta, lrr,Z0,term_id,ST_calc,vaso_INPUT);
                    
                    term_id--;
                }
                else
                {
                
//                      fprintf(stdout,"At vessel %d with dimensions %lf %lf\n",i,dimensions_matrix[i][0],    dimensions_matrix[i][1]);
                    
                    Arteries[i] = new Tube( dimensions_matrix[i][0], dimensions_matrix[i][1], dimensions_matrix[i][2],
                                           daughter_ptr, 0, num_pts, 0, 0, f1,f2,f3,0,0,0,0,0,0,0,0,ST_calc,vaso_INPUT);
                    
                }
            }
            
        }
        
    }
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
    
    
    
    // Solves the equations until time equals tend.///////////////////////////
    /* ADDED BY Mitchel J. Colebank
     * Rather than specifying the number of cycles as an input to the function,
     * we want to test to see if the solution has converged. If so, we should exit.*/
    
    int period_counter = 1; // Count the number of periods you have solved for
    double norm_sol = 1e+6;
    double sol_tol  = 1e-3;
    fflush(stdout);
    printf("NORM_SOL: %f\n",norm_sol);
    double sol_p1[tmstps],sol_p2[tmstps];
    tend      = Deltat;
    
    
    // SOLVE THE MODEL ONCE
    // Note: Only want to test the pressure at the inlet
    int sol_ID = 0;
    while (tend<=period_counter*Period)
    {
        solver (Arteries, tstart, tend, k, WALL_MODEL);
        sol_p1[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0],WALL_MODEL); // for printing
        sol_p1[sol_ID] *= rho*g*Lr/conv;
        tstart = tend;
        tend   = tend + Deltat; // The current ending time is increased by Deltat.
        sol_ID++;
    }
    
    
    // LOOP FOR CONVERGENCE
    double sse;
    while (norm_sol>=sol_tol && period_counter<60)
    {
        sol_ID = 0;
        sse    = 0;
        period_counter++;
        while (tend<=period_counter*Period)
        {
            solver (Arteries, tstart, tend, k, WALL_MODEL);
            sol_p2[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0],WALL_MODEL); // for printing
            sol_p2[sol_ID] *= rho*g*Lr/conv;
            sse = sse+ sq(sol_p1[sol_ID]-sol_p2[sol_ID]);
            tstart = tend;
            tend   = tend + Deltat; // The current ending time is increased by Deltat.
            sol_ID++;
        }
        norm_sol = sse;
        memcpy (sol_p1, sol_p2, sizeof(sol_p2));
        fflush(stdout);
        printf("NORM_SOL:%f\n",norm_sol);
    }
    fflush(stdout);
    printf("num_cylces:%d\n",period_counter);
    

  // The loop is continued until the final time
  // is reached. If one wants to make a plot of
  // the solution versus x, tend is set to final-
  // time in the above declaration.
    
  period_counter++;
    
    const char *namepuALL  = "pu_ALL.2d";
    
    FILE *fpALL = fopen (namepuALL, "w");
    
  finaltime = (period_counter+(cycles-1))*Period;
  while (tend <= finaltime)
  {
    for (int j=0; j<nbrves; j++)
    {
      int ArtjN = Arteries[j]->N;
      for (int i=0; i<ArtjN; i++)
      {
        Arteries[j]->Qprv[i+1] = Arteries[j]->Qnew[i+1];
        Arteries[j]->Aprv[i+1] = Arteries[j]->Anew[i+1];
      }
    }

    // Solves the equations until time equals tend.
    solver (Arteries, tstart, tend, k, WALL_MODEL);
	
     for (int save_id=0; save_id<nbrves; save_id++)
     {
         Arteries[ save_id] -> printAllt (fpALL, tend, 0,WALL_MODEL);
     }
    for (int save_id=0; save_id<nbrves; save_id++)
    {
        Arteries[ save_id] -> printAllt (fpALL, tend, Arteries[save_id]->N/2,WALL_MODEL);
    }
    for (int save_id=0; save_id<nbrves; save_id++)
    {
        Arteries[ save_id] -> printAllt (fpALL, tend, Arteries[save_id]->N,WALL_MODEL);
    }
    // The time within each print is increased.
    tstart = tend;
    tend   = tend + Deltat; // The current ending time is increased by Deltat.
  }
  
    // In order to termate the program correctly the vessel network and hence
    // all the vessels and their workspace are deleted.
    for (int i=0; i<nbrves; i++) delete Arteries[i];

    // Matrices and arrays are deleted
    for (int i=0; i<18; i++) delete[] fjac[i];
         
     fclose(fpALL);

    return 0;
}
