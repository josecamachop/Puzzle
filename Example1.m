%% First example: ADIEN. Arteaga et al. Solving the Puzzle
% simulation of multivariate data with control over the structure of 
% columns and rows using two-sided orthogonal Procrustes.
%
% Dependencies: 
%
%   - MEDA Toolbox v1.7 at https://github.com/codaslab/MEDA-Toolbox  
%
% designed by: Fransciso Arteaga (francisco.arteaga@ucv.es)   
%
% coded by: Jose Camacho (josecamacho@ugr.es)   
%
% last modification: 25/Mar/2025
%
% Copyright (C) 2025  
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Example matrix

clear
close all
clc

Y = [0.706	1.765	-1.267	-0.217	-0.121
1.347	0.711	0.586	1.524	1.280
-0.707	0.267	-0.325	-0.862	-1.504
-0.227	0.644	0.490	0.323	0.179
-0.577	-0.850	-0.064	-1.455	-1.037
0.532	0.156	-0.380	0.074	1.521
1.442	0.631	-1.220	0.574	-1.135
0.067	-1.047	-0.188	0.071	0.216
-1.063	-0.964	2.239	-1.305	0.574
-1.519	-1.312	0.129	1.274	0.026]

Cy = cov(Y)


%% Compute rounded eigenvalues and create covariance

eiV = sort(round(eig(Cy),1),'descend');
D = diag(eiV)


%% ADIEN

X = afamily(D,Y,'Method','ADIEN')


%% Validation

eiV_ADIEN = eig(cov(X))
MSE = norm(Y-X,'fro')^2/prod(size(Y))


%% ADIEN for 1 eigenvalue

eiV(2:end) = 0;
D = diag(eiV)

X = afamily(D,Y,'Method','ADIEN')


%% Validation

eiV_ADIEN_1e = eig(cov(X))
MSE_1e = norm(Y-X,'fro')^2/prod(size(Y))