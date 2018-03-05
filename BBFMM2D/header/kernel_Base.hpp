/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file kernel_Base.hpp
 Green's function or kernel
*/
#ifndef __kernel_Base_hpp__
#define __kernel_Base_hpp__

#include <iostream>
#include"H2_2D_Tree.hpp"

using namespace Eigen;

/*! Functions for calculating potential and virtual kernel function */
class kernel_Base {
public:
    void calculate_Potential(H2_2D_Tree& tree, double* potential , double* const charges);

    void calculate_Potential(H2_2D_Node*& node, MatrixXd& potential,H2_2D_Tree& tree);

    void set_Tree_Potential_Zero(H2_2D_Node* node);

    /*!	Obtains charge to node when needed;*/
    	void get_Charge(H2_2D_Node*& node , H2_2D_Tree& tree);

    /*! Obtains Chebyshev node potential from well separated clusters;*/
	  void calculate_NodePotential_From_Wellseparated_Clusters(H2_2D_Node*& node, unsigned short rank,unsigned short nChebNodes);

    /*!	Tranfers potential from node to final potential matrix when needed;*/
	void tranfer_Potential_To_Potential_Tree(H2_2D_Node*& node, MatrixXd& potential);

    /*!	Evaluate kernel at Chebyshev nodes;*/
	void kernel_Cheb_2D(const unsigned short& M, const vector<Point>& x, const unsigned short& N, const vector<Point>& y, MatrixXd& K);

    /*!	Evaluate the kernel;*/
    void kernel_2D(const unsigned long M, const vector<Point>& x, const unsigned long N, const vector<Point>& y, MatrixXd& kernel);

    /*!	Tranfers potential from Chebyshev node of parent to Chebyshev node of children;*/
	void transfer_NodePotential_To_Child(H2_2D_Node*& node,MatrixXd R[]);

    //! A pure virtural function
    /*! This function is used to define different types of kernel */
    virtual double kernel_Func(Point r0, Point r1) = 0;

    // Functions created by US
    void update_Charge (H2_2D_Tree& Tree , H2_2D_Node*& node);

    // Function created by US
    void set_Node_Charge_Zero (H2_2D_Node* node);

};


#endif //(__kernel_Base_hpp__)
