%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                   BLACK BOX FAST MULTIPOLE METHOD 2D                        %
%                             Version 2.0                                     %
%               Written for C++ by : Sivaram Ambikasaran, Ruoxi Wang          %
%        Written for MATLAB-Mex by : Judith Yue Li                            %
%       Modified for MATLAB-Mex by : Fahd Siddiqui, Ali Rezaei                %
%           https://github.com/DrFahdSiddiqui/bbFMM2D-MatlabMex               %
%                                                                             %
% =========================================================================== %
% LICENSE: MOZILLA 2.0                                                        %
%   This Source Code Form is subject to the terms of the Mozilla Public       %
%   License, v. 2.0. If a copy of the MPL was not distributed with this       %
%   file, You can obtain one at http://mozilla.org/MPL/2.0/.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOCUMENTATION                                                               %
%   Compiles the Mex file from the C++ code to be interfaced with MATLAB      %
%   Change the makefile parameters if needed                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename='FMM_MatVec';
srcTree = 'BBFMM2D/src/H2_2D_Tree.cpp';
srcNode = 'BBFMM2D/src/H2_2D_Node.cpp';
srcKernelBase = 'BBFMM2D/src/kernel_Base.cpp';
main = 'mexFMM2D.cpp';
eigenDIR = '/usr/include/eigen3/';
header = 'BBFMM2D/header/';
compiler_option = '-c -std=c++11 -fopenmp -w -Wall -DNDEBUG -O3 -ffast-math -ffinite-math-only ';
linker_option = '-g -fopenmp';
cxxoptim_option = ' -fopenmp -O3 -ffast-math -ffinite-math-only';
mex('-O',main,srcTree,srcNode,srcKernelBase,'-largeArrayDims', ...
    ['-I',eigenDIR],['-I',header],'-output',filename,...
    ['COMPFLAGS="$COMPFLAGS ' compiler_option '"'],...
    ['CXXOPTIMFLAGS="$CXXOPTIMFLAGS ' cxxoptim_option '"'],...
    ['LDFLAGS="$LDFLAGS ' linker_option '"']);
fprintf('A MEX function %s.mex is ready for use\n',filename);
fprintf('Run FMM_Main.m to perform computation\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
