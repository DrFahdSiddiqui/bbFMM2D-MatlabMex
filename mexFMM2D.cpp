/******************************************************************************
 *                                                                            *
 *                   BLACK BOX FAST MULTIPOLE METHOD 2D                       *
 *                             Version 2.0                                    *
 *               Written for C++ by : Sivaram Ambikasaran, Ruoxi Wang         *
 *        Written for MATLAB-Mex by : Judith Yue Li                           *
 *       Modified for MATLAB-Mex by : Fahd Siddiqui, Ali Rezaei               *
 *           https://github.com/DrFahdSiddiqui/bbFMM2D-MatlabMex              *
 *                                                                            *
 * ========================================================================== *
 * LICENSE: MOZILLA 2.0                                                       *
 *   This Source Code Form is subject to the terms of the Mozilla Public      *
 *   License, v. 2.0. If a copy of the MPL was not distributed with this      *
 *   file, You can obtain one at http://mozilla.org/MPL/2.0/.                 *
 *****************************************************************************/

/******************************************************************************
 * DOCUMENTATION                                                              *
 *   Main entry point for MATLAB to C/C++ interface                           *
 *   No changes needed                                                        *
 *   Run make.m if any changes are made                                       *
 *****************************************************************************/

/*****************************************************************************/

#include <iostream>
#include "math.h"
#include "mex.h"
#include "matrix.h"
#include <Eigen/Core>
#include "H2_2D_Tree.hpp"
#include "H2_2D_Node.hpp"
#include "kernel_Base.hpp"
#include "kernelfun.hpp"
#include "class_handle.hpp"
#include <omp.h>
#include <ctime>
#include <ratio>
#include <chrono>
//const int NUM_THREADS=16;

using namespace Eigen;
using namespace std;
using namespace std::chrono;


extern void _main();

#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_SCALAR(P) (IS_REAL_2D_FULL_DOUBLE(P) && mxGetNumberOfElements(P) == 1)

/* ------------------------------------------------------------------------- */
// Pass location from matlab to C
void read_location(const mxArray* x, const mxArray* y, vector<Point>& location){
    unsigned long N;
    double *xp, *yp;
    N = mxGetM(x);
    xp = mxGetPr(x);
    yp = mxGetPr(y);
    for (unsigned long i = 0; i < N; i++){
        Point new_Point;
        new_Point.x = xp[i];
        new_Point.y = yp[i];
        location.push_back(new_Point);
    }
}

/* ------------------------------------------------------------------------- */
//MATLAB mexFunction interface
void mexFunction(int nlhs,mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	omp_set_num_threads(2);
	// Get the command string
	char cmd[64];
	if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");


    // BUILDING THE FMM TREE -----------------------------------------
    if (!strcmp("new", cmd)) {
		// Output Arguments
		// Tree instance          plhs[0]			   new H2_2D_Tree(nChebNodes, location, N, m)

		// Input Argument
		// String case       	  		prhs[0] 			 New
		//Location x              		prhs[1]				location
		//Location y         		    prhs[2] 			location
		// Sets of charges           prhs[3]			m
		//nCheb Nodes   			prhs[4]				nChebNodes
		//PrintFlag					     prhs[5]			print


		unsigned long N;     		unsigned m;
		N = mxGetM(prhs[1]);		m = *mxGetPr(prhs[4]);
		vector<Point> location;	  read_location(prhs[1],prhs[2],location);
		unsigned short nChebNodes = *mxGetPr(prhs[3]);
		bool print = *mxGetPr(prhs[5]);

		cout <<  "Number of charges:"  << N  << endl;
		cout <<  "Number of sets of charges:" << m << endl;
		cout <<  "Number of Chebyshev Nodes: " << nChebNodes << endl;

		// 1. Build Tree
		high_resolution_clock::time_point tic_BuildTree = high_resolution_clock::now();
			plhs[0] = convertPtr2Mat<H2_2D_Tree>(new H2_2D_Tree(nChebNodes, location, N, m));
		high_resolution_clock::time_point toc_BuildTree = high_resolution_clock::now();

    	duration<double> time_BuildTree = duration_cast<duration<double>>(toc_BuildTree - tic_BuildTree);
        if(print)
			mexPrintf("\nTime taken for FMM(build tree) is: %.4g\n", time_BuildTree.count());
        return;
    }

    // DELETE THE FMM TREE -----------------------------------------
    if (!strcmp("delete", cmd)) {
        // Destroy the C++ object
        destroyObject<H2_2D_Tree>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }


	H2_2D_Tree *Atree = convertMat2Ptr<H2_2D_Tree>(prhs[1]);

    // FMM MATRIX VECTOR PRODUCT -----------------------------------------
	if (!strcmp("CalcPot", cmd)) {

		// Output Arguments
		// QH_OUT          plhs[0]			   QHp

		// Input Argument
		// String case       prhs[0] 			 CalcPot
		// Tree					prhs[1] 			ATree
		// Kernel Name     prhs[2] 			   kernelNameC
		// Charges   		 prhs[3] 			 charges
		// Print Flag 		  prhs[4] 			   print

		bool print = *mxGetPr(prhs[4]);

		char* kernelNameC; kernelNameC=mxArrayToString(prhs[2]);
		std::string kernelName(kernelNameC);
		double *charges;			charges = mxGetPr(prhs[3]);

		high_resolution_clock::time_point tic_CalcPot = high_resolution_clock::now();

		plhs[0] = mxCreateDoubleMatrix(Atree->N, Atree->m, mxREAL);
		double *QHp;  					QHp = mxGetPr(plhs[0]);

		//kernelList kernelChoice=kernelName;
		kernel_Base* kernel=NULL;
		selectKernel(kernelName, kernel);

		kernel->calculate_Potential(*Atree, QHp, charges);
		high_resolution_clock::time_point toc_CalcPot = high_resolution_clock::now();
    	duration<double> time_CalcPot = duration_cast<duration<double>>(toc_CalcPot - tic_CalcPot);

		if(print)
        	mexPrintf("Time taken for FMM(calculating potential) is: %.4g\n",time_CalcPot.count());

		delete kernel;
		return;
    }

    // EXACT MATRIX VECTOR PRODUCT -----------------------------------------
    if (!strcmp("ExactPot", cmd)) {

		// Output Arguments
		// QH_OUT          plhs[0]			   QHp

		// Input Argument
		// String case       prhs[0] 			 ExactPot
		// Kernel Name    prhs[2] 			   kernelNameC
		// Location x        prhs[3]				location
		// Location y        prhs[4] 			location
		// Charges   		 prhs[5] 			 charges
		// Print Flag 		  prhs[6] 			   print


		unsigned long N;
		unsigned long m;
		vector<Point> location;
		char* kernelNameC;
		double *charges;

		N = mxGetM(prhs[5]);
		m = mxGetN(prhs[5]);
		read_location(prhs[3],prhs[4],location);
		bool print = *mxGetPr(prhs[6]);

		kernelNameC=mxArrayToString(prhs[2]);
		std::string kernelName(kernelNameC);

		charges = mxGetPr(prhs[5]);

		high_resolution_clock::time_point tic_ExactPot = high_resolution_clock::now();

		plhs[0] = mxCreateDoubleMatrix(N, m, mxREAL);
		double *QHp;  					QHp = mxGetPr(plhs[0]);

		//kernelList kernelChoice=kernelName;
		kernel_Base* kernel=NULL;
		selectKernel(kernelName, kernel);

        clock_t start = clock();
        MatrixXd Q;
        kernel->kernel_2D(N, location, N, location, Q);// Q is initialized inside function A.kernel_2D

        // Compute exact Matrix vector product
        Map<MatrixXd> QHT(QHp,N,m);
		MatrixXd H = Map<MatrixXd>(charges, N, m);
        QHT = Q*H;
		high_resolution_clock::time_point toc_ExactPot = high_resolution_clock::now();
    	duration<double> time_ExactPot = duration_cast<duration<double>>(toc_ExactPot - tic_ExactPot);
        if(print)
            mexPrintf("\nThe total exact computation time is: %.4g\n",time_ExactPot.count());
		delete kernel;
    }

    return;

}

/*****************************************************************************/
