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
%   Main script to be called in MATLAB                                        %
%   Set the flags, load input data and choose the kernel name from kernel.hpp %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% Set Parameters
PrintFlag = true;       % Print computation time and error
TestingMode = true;    % Perfom exact computation (set to true)

%% GET INPUT DATA FROM DATA FILES ------------------------------------------- %
% Read data from input file
Data       = load('Input/input.txt');
location   = Data(:,1:2);           % Locations of the charges matrix
charges    = Data(:,3:end);         % Sets of Charges
nChebNodes = 5;                     % Number of Chebyshev nodes( >= 3)
m=size(charges,2);
xloc = location(:,1);
yloc = location(:,2);


%% FAST MATRIX VECTOR PRODUCT ----------------------------------------------- %
% Calls the MEX Function FMM_MatVec

Tree=FMMTree(location, nChebNodes, m, PrintFlag);

kernelName='ex2';  % Choose kernel name from kernel.hpp

[QH] = FMMCalcPot(Tree, kernelName, charges, PrintFlag);

charges=ones(size(charges))*2000;
[QH] = FMMCalcPot(Tree, kernelName, charges, PrintFlag);

if TestingMode
    [QHE] = Tree.FMMExactPot( kernelName, location, charges, PrintFlag );
    fprintf('\n Maximum Error is: %0.3e \n', ...
               norm(QHE-QH)/norm(QHE));
end
%[QH,QHexact] = FMM_MatVec(location(:,1), location(:,2), charges,nChebNodes,kernelName,PrintFlag);

%clear Tree;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
