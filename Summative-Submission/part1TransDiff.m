% This script declares the variables and solves the transiet diffusion
% equation (as per the conditions in Part 1 of the summative CW)

% Set the title of the plot
titlestr = "Transient-Diffusion";

% Material/equation coefficients
D = @(x,t) 1; % Declare the diffusion function in x and t
lambda = @(x,t) 0; % Declare the reaction function in x and t
f = @(x,t) 0; % Declare the source function in x and t

% Initial and Boundary conditions
ICs = @(x) 0; % Declare the initials conditions as a function of x
BC0 = @(t) 0; % Declare the LHS boundary condition as function of t
BC0type = "Dirichlet"; % Declare the LHS boundary condition type
BC1 = @(t) 1; % Declare the RHS boundary condition as function of t
BC1type = "Dirichlet"; % Declare the RHS boundary condition type

% Time and mesh parameters
xmin = 0; % The lower x bound of the domain to solve over
xmax = 1; % The upper x bound of the domain to solve over
ne = 10; % The  number  of elements in the 1D mesh
dt = 0.05; % The size of the transient time steps
time = 1; % The time over which the solution should be determined
theta = 1; % Declare which time stepping scheme to be used

% Solve this problem
mesh = OneDimMeshGen(xmin,xmax,ne,1); % Generate the mesh
[solution,tvec,xvec] = transientFEMSolver(titlestr,mesh,theta,dt,time, ...
    D,lambda,f,ICs,BC0,BC0type,BC1,BC1type);

% Generate any other plots
figure;
hold on
treq = [0.05,0.1,0.3,1.0];
selectivePlot(treq,xvec,solution,dt,titlestr,"n","X","t");
% Plot the analytical solution
astep = 0.01;
atvec = 0:astep:1;
axvec = 0:astep:1;
asolution = zeros(length(axvec),length(atvec));
for i = 1:length(atvec)
    for j = 1:length(axvec)
        asolution(j,i) = TransientAnalyticSoln(axvec(j),atvec(i));
    end
end
selectivePlot(treq,axvec,asolution,astep,titlestr,"a","X","t");
hold off

% Plot the analytical and numerical at x=0.8
figure;
hold on
xreq=[0.8];
xstep=(xmax-xmin)/ne;
selectivePlot(xreq,tvec,(solution)',xstep,titlestr,"n","Time [s]","x");
selectivePlot(xreq,atvec,(asolution)',astep,titlestr,"a","Time [s]","x");
hold off



