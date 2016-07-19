%2D Truss(f=ku) 
%Kevin(Yinuo) Huang CID:01051134 04:20-07:32 10/04/2016

%******** 1.IMPORT DATA(feel no need for checking data format 'cause users are intelligent)
%******** 2.DATA PROCESSOR+CHECK IF THE TRUSS STATISTICALLY DETERMINATES
%******** 3.GET K-BAR
%******** 4.ASSEMBLE GLOBAL K
%******** 5.INITILISATION
%******** 6.FREE NODE DISPLACEMENTS
%******** 7.REACTIONS AND BAR FORCES
%******** 8.OUTPUT

clear
clc
close all;

inputfilename=uigetfile('*.txt','Select the input file please.');
%------------------------------------------------------------------%
fid=fopen(inputfilename,'r');

%n=number of nodes
%m=number of bars
n=fscanf(fid,'%f',1);
m=fscanf(fid,'%f\n',1);

%[trussnodes]=data for all nodes
%[trussbars]=data for all bars
trussnodes=fscanf(fid,'%f%f%f%f%f%f%f',[7,n]);
trussnodes=trussnodes';
trussbars=fscanf(fid,'%f%f%f%f%f',[5,m]);
trussbars=trussbars';

fclose(fid);
%--------------check statistically (in)determinates----------------------------------------%
r=sum(trussnodes(:,4))+sum(trussnodes(:,5));
%maxwell equation
M=m*r-2*n;
switch M
    case M<0 
        error('This is a mechanism rather than a truss!');
    case M>0
        error('The truss is statistically indeterminates.(redundant)');
end

%-----------------------------------------------------------------------%
%initilise the stiffness
globalk=zeros(2*n,2*n);

for i=1:m
    %prepare index to get k bar
    q=trussbars(i,2);
    w=trussbars(i,3);
    node1=trussnodes(q,2:3);
    node2=trussnodes(w,2:3);
    A=trussbars(i,4);
    E=trussbars(i,5);
    [k,theta,L]=getKbar(node1,node2,A,E);
    
    %lengths of the bars
    trussbars(i,6)=L;
    %angles
    trussbars(i,7)=theta;
    %basic k needed to calculate the bar forces
    trussbars(i,8)=trussbars(i,4)*trussbars(i,5)/L;
    
    %define x-y location of stiffness
    klocation1=2*(trussbars(i,2)-1)+[1,2];
    klocation2=2*(trussbars(i,3)-1)+[1,2];
    
    %assemble the stiffness
    globalk(klocation1,klocation1)=globalk(klocation1,klocation1)+k(1:2,1:2);
    globalk(klocation1,klocation2)=globalk(klocation1,klocation2)+k(1:2,3:4);
    globalk(klocation2,klocation1)=globalk(klocation2,klocation1)+k(3:4,1:2);
    globalk(klocation2,klocation2)=globalk(klocation2,klocation2)+k(3:4,3:4);    
end
%-------------initialise some matrix-------------------------------------------------%
%initialise the displacement force matrix
u=zeros(2*n,1);
f=zeros(m,1);
p=[-1 0 1 0];

%assign constraints,pre-loads
constraints=reshape(trussnodes(:,4:5)',[],1);
constraintloads=reshape(trussnodes(:,6:7)',[],1);

%initialise support reactions and fixed node forces
supports=zeros(2*n,1);
fixednodeforces=zeros(2*n,1);

%-----free node displacements+support reactions--------------------%

%find locations of free and supp.
freenodelocations=find(~constraints);
supportnodelocations=find(constraints);

%assign supp f, fixed f and corresponding stiffness 
suppf=constraintloads(freenodelocations);
fixednodef=fixednodeforces(freenodelocations);
kfreecurrent=globalk(freenodelocations,freenodelocations);

%free node displacements
freenodeu=kfreecurrent\(suppf-fixednodef);
u(freenodelocations)=freenodeu;

%support reactions
ksuppcurrent=globalk(supportnodelocations,freenodelocations);
currentsuppreaction=ksuppcurrent*freenodeu;
supports(supportnodelocations)=currentsuppreaction;

%-----bar forces------------------------%

for i=1:m
    %define two ends of a bar
    startlocation=trussbars(i,2);
    endlocation=trussbars(i,3);
    
    %same way as assembly the stiffness
    klocation1=2*(trussbars(i,2)-1)+[1,2];
    klocation2=2*(trussbars(i,3)-1)+[1,2];
    
    %grab angle and temp matrix
    angle=trussbars(i,7);
    temp=[cos(angle),sin(angle);-sin(angle),cos(angle)];
    tempzero=[temp,temp.*0;temp.*0,temp];
    
    %compute the force according to handout
    f(i)=trussbars(i,8)*p*tempzero*u([klocation1,klocation2]);
end

%-----------------output results table-----------------------------&=%
trussbars=[trussbars,f];
trussnodes(:,4:5)=reshape(supports,2,[])';
trussnodes=[trussnodes,reshape(u,2,[])'];

Columnheadersofnodes={'Bar','coordinateX','coordinateY','Xreaction','Yreacction','Xpreload','Ypreload','Xdisplacement','Ydisplacement'};
Columnheadersofbars={'Bar','Start','End','A','E','Length','Angle','XBarForce','YBarForce'};

nodetable=array2table(trussnodes,'VariableNames',Columnheadersofnodes);
bartable=array2table(trussbars,'VariableNames',Columnheadersofbars);

disp(nodetable);
disp(bartable);

%USE THE FOLLOWING IF YOU WANT A TABLE IN TXT FORMAT
%writetable(nodetable);
%type nodetable.txt;
%writetable(bartable);
%type bartable.txt;














































