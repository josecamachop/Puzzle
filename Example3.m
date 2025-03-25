%% Third example: continuum between ARDICO and ACORDI. Arteaga et al. Solving the Puzzle
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
% Copyright (C) 2024
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


%% Continuum for 0.25 and 0.75

X025 = afamily(C,Y,'Method','Continuum','Gamma',0.25)

X075 = afamily(C,Y,'Method','Continuum','Gamma',0.75)


%% Validation

MSE_025 = norm(Y-X025,'fro')^2/prod(size(Y))
MSGMD_025 = norm(Y*Y'-X025*X025','fro')^2/(size(Y,1)^2)
MSCD_025 = norm(cov(X025)-C,'fro')^2/prod(size(C))
MSDD_025 = norm(dist(Y').^2-dist(X025').^2,'fro')^2/(size(Y,1)*(size(Y,1)-1)/2)

MSE_075 = norm(Y-X075,'fro')^2/prod(size(Y))
MSGMD_075 = norm(Y*Y'-X075*X075','fro')^2/(size(Y,1)^2)
MSCD_075 = norm(cov(X075)-C,'fro')^2/prod(size(C))
MSDD_075 = norm(dist(Y').^2-dist(X075').^2,'fro')^2/(size(Y,1)*(size(Y,1)-1)/2)


%% For a range of values

MSE = [];
MSGMD = [];
MSCD = [];
MSDD = [];
for gamma = 0:0.1:1
    X = afamily(C,Y,'Method','Continuum','Gamma',gamma);
    MSE = [MSE norm(Y-X,'fro')^2/prod(size(Y))];
    MSGMD = [MSGMD norm(Y*Y'-X*X','fro')^2/(size(Y,1)^2)];
    MSCD = [MSCD norm(cov(X)-C,'fro')^2/prod(size(C))];
    MSDD = [MSDD norm(dist(Y').^2-dist(X').^2,'fro')^2/(size(Y,1)*(size(Y,1)-1)/2)];
end

plotVec(MSE,'XYLabel', {'Gamma','MSE'}, 'EleLabel', 0:0.1:1, 'PlotType', 'Lines');
axis([0 1 0.56 0.57])
saveas(gcf,'Fig/MSE.jpg','jpg')

plotVec(MSGMD', 'XYLabel', {'Gamma','FDD'}, 'EleLabel', 0:0.1:1, 'PlotType', 'Lines');
saveas(gcf,'Fig/MSGMD.jpg','jpg')

plotVec(MSCD', 'XYLabel', {'Gamma','MSCD'}, 'EleLabel', 0:0.1:1, 'PlotType', 'Lines');
saveas(gcf,'Fig/MSCD.jpg','jpg')

plotVec(MSDD', 'XYLabel', {'Gamma','MSDD'}, 'EleLabel', 0:0.1:1, 'PlotType', 'Lines');
saveas(gcf,'Fig/MSDD.jpg','jpg')