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

%CLASS_INTERFACE Example MATLAB class wrapper to an underlying C++ class
classdef FMMTree < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance
        function this = FMMTree(location, nChebNodes, m, print)
            this.objectHandle = FMM_MatVec('new', location(:,1), location(:,2), nChebNodes, m, print);
        end

        %% Destructor - Destroy the C++ class instance
        function delete(this)
            FMM_MatVec('delete', this.objectHandle);
        end

        %% Calculate Potential - an example class method call
        function [QH] = FMMCalcPot(this, kernelName, charges, print)
            [QH] = FMM_MatVec('CalcPot', this.objectHandle, kernelName, charges, print);
        end


        %% Calculate Potential - an example class method call
        function [QH] = FMMExactPot(this, kernelName, location, charges, print )
            [QH] = FMM_MatVec('ExactPot', this.objectHandle, kernelName, location(:,1), location(:,2), charges, print);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
