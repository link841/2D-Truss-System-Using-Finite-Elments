function [K,theta,L]=getKbar(node1,node2,A,E)
%This function returns the stiffness matrix for a truss element connecting
%nodes 1 and node 2 such that
%[Fx1;Fy1;Fx2;Fy2]=K*[u1;v1;u2;v2]
%
%[K]=getKmat(node1,node2,A,E)
%
%outputs: K     -   stiffness matrix of dimensions 4x4
%inputs:  node1 -   coordinates of node 1 [x1,y1] (m)
%         node2 -   coordinates of node 2 [x2,y2] (m)
%         A     -   bar cross-sectional area (m^2)
%         E     -   bar Youngs modulus (GPa)
%Written by E Levis, February 2016

dist=node2-node1; 
L=sqrt(sum(dist.^2)); %caluclate bar length
theta=atan2(dist(2),dist(1)); %caluclate bar angle
tempA=[cos(theta)^2,sin(theta)*cos(theta);sin(theta)*cos(theta),sin(theta)^2];
K=E*A*[tempA,-tempA;-tempA,tempA]*10^9/L; %stiffness matrix for single bar

%end of function