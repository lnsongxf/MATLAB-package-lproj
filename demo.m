% Code to estimate Smooth Local Projections
% R Barnichon and C Brownlees, 05/2017
% Preliminary version

data = csvread('data.csv',1,1);

T = size(data,2);
P  = 4; %number of lags in LP

% IRF 1
h1 = 1; %start LP at h=0 or h=1 (h=1 if impose no contemporanous impact)
H  = 20; %#horizon of IR

y  = data(:,1); %endogenous variable
x  = data(:,3); %shock variable
w  = [ data(:,1:2) , lagmatrix( data , 1:P ) ]; %control variables (contemporaneous vars, lagged vars)
w( ~isfinite(w) ) = 0;

lp = locproj(y,x,w,h1,H,'reg'); %IR from Local Projection

r = 3; %(r-1)=order of the limit polynomial (so r=3 implies the IR is shrunk towards a quadratic polynomial)

lambda = 100; %value of penalty
slp    = locproj(y,x,w,h1,H,'smooth',r,lambda); %IR from Smooth Local Projection

slp_lim= locproj(y,x,w,h1,H,'smooth',r,1e10); % Limit IR in Smooth Local Projection

figure(1)
hold on,
plot( 0:H , [ lp.IR slp.IR slp_lim.IR] )
plot( 0:H , zeros(H+1,1) , '-k' , 'LineWidth' , 2 )
grid
xlim([0 H])
legend('IR_{lp}','IR_{slp}','IR_{slp,max pen}','Location','Best')

% IRF2
h1 = 0; 
H  = 20; %#horizon of IR

y  = data(:,3); %endogenous variable
x  = data(:,1); %shock variable
w  = [ lagmatrix( data , 1:P ) ]; %control variables (contemporaneous vars, lagged vars)
w( ~isfinite(w) ) = 0;

lp = locproj(y,x,w,h1,H,'reg'); %IR from Local Projection

r = 3; %(r-1)=order of the limit polynomial (so r=3 implies the IR is shrunk towards a quadratic polynomial)

lambda = 1000; %value of penalty
slp    = locproj(y,x,w,h1,H,'smooth',r,lambda); %IR from Smooth Local Projection

slp_lim= locproj(y,x,w,h1,H,'smooth',r,1e10); % Limit IR in Smooth Local Projection

figure(2)
hold on,
plot( 0:H , [ lp.IR slp.IR slp_lim.IR] )
plot( 0:H , zeros(H+1,1) , '-k' , 'LineWidth' , 2 )
grid
xlim([0 H])
legend('IR_{lp}','IR_{slp}','IR_{slp,max pen}','Location','Best')

