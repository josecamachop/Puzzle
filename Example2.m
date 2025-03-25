%% Second example: ARDICO and ACORDI. Arteaga et al. Solving the Puzzle
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

C = [1.0	0.6	0.1	0.4	-0.3
0.6	1.0	-0.4	0.0	-0.6
0.1	-0.4	1.0	0.1	0.6
0.4	0.0	0.1	1.0	0.1
-0.3	-0.6	0.6	0.1	1.0]


%% ARDICO, ACORDI and ADICOV

X1 = afamily(C,Y,'Method','ARDICO')

X2 = afamily(C,Y,'Method','ACORDI')

disp('Method is ADICOV')
X3 = ADICOV(9*C,Y,5)


%% Validation

MSE_ARDICO = norm(Y-X1,'fro')^2/prod(size(Y));
MSGMD_ARDICO = norm(Y*Y'-X1*X1','fro')^2/(size(Y,1)^2);
MSCD_ARDICO = norm(cov(X1)-C,'fro')^2/prod(size(C));
MSDD_ARDICO = norm(dist(Y').^2-dist(X1').^2,'fro')^2/(size(Y,1)*(size(Y,1)-1)/2);
ARDICO = [MSE_ARDICO MSGMD_ARDICO MSCD_ARDICO MSDD_ARDICO];

MSE_ACORDI = norm(Y-X2,'fro')^2/prod(size(Y));
MSGMD_ACORDI = norm(Y*Y'-X2*X2','fro')^2/(size(Y,1)^2);
MSCD_ACORDI = norm(cov(X2)-C,'fro')^2/prod(size(C));
MSDD_ACORDI = norm(dist(Y').^2-dist(X2').^2,'fro')^2/(size(Y,1)*(size(Y,1)-1)/2);
ACORDI = [MSE_ACORDI MSGMD_ACORDI MSCD_ACORDI MSDD_ACORDI];

MSE_ADICOV = norm(Y-X3,'fro')^2/prod(size(Y));
MSGMD_ADICOV = norm(Y*Y'-X3*X3','fro')^2/(size(Y,1)^2);
MSCD_ADICOV = norm(cov(X3)-C,'fro')^2/prod(size(C));
MSDD_ADICOV = norm(dist(Y').^2-dist(X3').^2,'fro')^2/(size(Y,1)*(size(Y,1)-1)/2);
ADICOV = [MSE_ADICOV MSGMD_ADICOV MSCD_ADICOV MSDD_ADICOV];

T = table({'MSE','MSGMD','MSCD','MSDD'}', ARDICO', ACORDI', ADICOV','VariableNames', {'Metric','ARDICO','ACORDI','ADICOV'})


