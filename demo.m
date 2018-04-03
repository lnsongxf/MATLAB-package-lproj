% Smooth Local Projections Demo
% R. Barnichon and C. Brownlees, 04/2018

%% LOAD DATA
data = csvread('data.csv',1,1);

T = size(data,2);
P  = 4; % number of lags used in LP for controls

%% ESTIMATING THE IR OF GDP TO A MONETARY SHOCK

H_min = 1; % start LP at H_min=0 or 1 (H_min=1 if impose no contemporanous impact)
H_max = 20; 

y  = data(:,1); % endogenous variable
x  = data(:,3); % endoegnous variable related to the shock 
w  = [ data(:,1:2) , lagmatrix( data , 1:P ) ]; % control variables (contemporaneous vars, lagged vars)
w( ~isfinite(w) ) = 0;

lp = locproj(y,x,w,H_min,H_max,'reg'); % IR from (standard) Local Projection

r = 2; %(r-1)=order of the limit polynomial (so r=2 implies the IR is shrunk towards a line polynomial)
lambda = 100; % value of penalty

slp    = locproj(y,x,w,H_min,H_max,'smooth',r,lambda); %IR from Smooth Local Projection
slp_lim= locproj(y,x,w,H_min,H_max,'smooth',r,1e10); % Limit IR in Smooth Local Projection

figure(1)
hold on,
plot( 0:H_max , [ lp.IR slp.IR slp_lim.IR] )
plot( 0:H_max , zeros(H_max+1,1) , '-k' , 'LineWidth' , 2 )
grid
xlim([0 H_max])
legend('IR_{lp}','IR_{slp}','IR_{slp,max pen}','Location','Best')

%% ESTIMATING THE IR OF FFR TO A MONETARY SHOCK
H_min  = 0; 
H_max  = 20;

y  = data(:,3); % endogenous variable
x  = data(:,1); % shock variable
w  = [ lagmatrix( data , 1:P ) ]; % control variables (contemporaneous vars, lagged vars)
w( ~isfinite(w) ) = 0;

lp = locproj(y,x,w,H_min,H_max,'reg'); % IR from Local Projection

r      = 3; % (r-1)=order of the limit polynomial (so r=3 implies the IR is shrunk towards a quadratic polynomial)
lambda = 1000; % value of penalty

slp    = locproj(y,x,w,H_min,H_max,'smooth',r,lambda); % IR from Smooth Local Projection
slp_lim= locproj(y,x,w,H_min,H_max,'smooth',r,1e10); % Limit IR in Smooth Local Projection

figure(2)
hold on,
plot( 0:H_max , [ lp.IR slp.IR slp_lim.IR] )
plot( 0:H_max , zeros(H_max+1,1) , '-k' , 'LineWidth' , 2 )
grid
xlim([0 H_max])
legend('IR_{lp}','IR_{slp}','IR_{slp,max pen}','Location','Best')

%% Confidence Intervals

H_min = 1;  % start LP at H_min=0 or 1 (H_min=1 if impose no contemporanous impact)
H_max = 20; 

y  = data(:,1); % endogenous variable
x  = data(:,3); % endoegnous variable related to the shock 
w  = [ data(:,1:2) , lagmatrix( data , 1:P ) ]; % control variables (contemporaneous vars, lagged vars)
w( ~isfinite(w) ) = 0;

lp = locproj(y,x,w,H_min,H_max,'reg'); % IR from (regular) Local Projection
lp = locproj_conf(lp,H_max); % it take a bit too run! please be patiante

figure(3)
hold on,
plot( 0:H_max , lp.IR   , 'r' , 'LineWidth' , 2 )
plot( 0:H_max , lp.conf , 'r' )
plot( 0:H_max , zeros(H_max+1,1) , '-k' , 'LineWidth' , 2 )
grid
xlim([0 H_max])

r      = 2;
lambda = 100;
slp    = locproj(y,x,w,H_min,H_max,'smooth',r,lambda); %IR from Smooth Local Projection
slp    = locproj_conf(slp,H_max,lambda);

figure(4)
hold on,
plot( 0:H_max , slp.IR   , 'r' , 'LineWidth' , 2 )
plot( 0:H_max , slp.conf , 'r' )
plot( 0:H_max , zeros(H_max+1,1) , '-k' , 'LineWidth' , 2 )
grid
xlim([0 H_max])