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
 *   Define multiple kernel functions here in C/C++ syntax                    *
 *   Update the list in the function selectKerne                              *
 *   Run make.m if any changes are made                                       *
 *****************************************************************************/

/*****************************************************************************/
#include "matrix.h"

/* ------------------------------------------------------------------------- */
//Class ex1
class ex2: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1){
        double X =   (r0.x - r1.x) ;
        double Y =   (r0.y - r1.y) ;
        double r = sqrt((X*X) + (Y*Y));// Define kernel functio  here
        return r;
    }
};

/* ------------------------------------------------------------------------- */
//Class pdn
class pdn: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1){
        double X =   (r0.x - r1.x) ;
        double Y =   (r0.y - r1.y) ;
        double r =   ((X*X) + (Y*Y));// Define kernel functio  here
        return r;
    }
};


/* ------------------------------------------------------------------------- */
//Update the kernel list here
void selectKernel(string name, kernel_Base*& kernel){
    if (name=="ex2")     kernel=new ex2;
    if (name== "pdn")    kernel=new pdn;
    // Update kernel function  here to be used in FMM_Main.m
}

/*****************************************************************************/
