function [solution,tvec,xvec] = transientFEMSolver(titles,mesh,theta,dt,...
    endtime,Df,lambdaf,ff,ICs,BC0,BC0type,BC1,BC1type)
%TRANSIENTFEMSOLVER solves a transient, 1D FEM problem in 1 Dimension
%   This function returns the numerical solution for a diffusion-reaction
%   equation, in a transient, 1D case. This function takes in the following
%   inputs:
%       -title: a string containing the 'title' of the problem being solved
%               for example: "Transient-Diffusion"
%       -mesh: in instance of the mesh data type, created by calling either
%              OneDimLinearMeshGen or OneDimRefinedMeshGen functions
%       -theta: declare the time stepping scheme, the values are:
%                   -theta=1/2: Crank-Nicolson
%                   -theta=0: Forward Euler
%                   -theta=1: Backward Euler
%       -dt: the size of the time step for the transient problem
%       -time: the total time over which the problem should be solved
%       -Df: the diffusion term, as a function handle (D = @(x,t))
%       -lambdaf: the reaction term, as function handle (lamdba = @(x,t))
%       -ff: the source term, as a function handle (f = @(x,t))
%       -ICs: the initial system conditions, as a function handle in x
%             (ICs = f(x))
%       -BC0: the boundary conditions at the left most x position as a 
%             function handle in t (BC0 = @(t))
%       -BC0type: the type of the boundary condition at the left most x
%                 position (either "Neumann" or "Dirichlet")
%       -BC1: the boundary conditions at the right most x position as a 
%             function handle in t (BC1 = @(t))
%       -BC1type: the type of the boundary condition at the right most x
%                 position (either "Neumann" or "Dirichlet")
%   The function returns 3 outputs:
%       -solution: a matrix of the solution at each mesh x position for 
%                  each time step. Dimensions (mesh.ngn x nt) being the 
%                  number of global nodes and the number of time steps
%       -tvect: a column vector of all the time values 
%       -xvect: a column vector of all the x values
%   This function also outputs a surface plot for the solution

% Create the first ccurr vector corresponding to the initial conditions
ccurr = initialConditions(ICs,mesh);
% Define cnext as a column vector of zeros
cnext = zeros(mesh.ngn,1);

% Initialise the time variable and number of time steps
time = 0;
N = (endtime/dt);
tvec = zeros(N+1,1);

% Set the output xvec
xvec = mesh.nvec;

% Set the order of the Lagrange basis functions being used 
order = mesh.order;

% Initialise the solution matrix 
solution = zeros(mesh.ngn,N+1);
% Add the first solution into the solution matrix 
solution(1:mesh.ngn,1) = ccurr;

% Loop over all the time steps 
for tstep = 1:N
    
    % Set all the material coefficients at this timestep
    D = @(x) Df(x,time);
    lambda = @(x) lambdaf(x,time);
    ffcurr = @(x) ff(x,time);
    ffnext = @(x) ff(x,time+dt);

    % Create the global mass matrix
    M = globalMassMatrix(mesh,order);

    % Create the global stiffness matrix 
    K = globalStiffnessMatrix(D,lambda,mesh,order);

    % Create global matrix ( M + theta * dt * K )
    gM = M + (theta*dt*K);

    % Create the matrix to multiply prev cvec ( M - (1-theta)*dt*K )
    P = M - (1-theta)*dt*K;

    % Multiply current sulition (M-(1-theta)*dt*K)*ccurr store in global
    % vector! )
    gV = P * ccurr;
    
    % Create the global vector ( dt * ( theta * Fnext + (1-theta)*Fcurr ) )
    Fcurr = globalSourceVector(ffcurr,mesh,order);
    Fnext = globalSourceVector(ffnext,mesh,order);
    gV = gV + dt*(theta*Fnext + (1-theta)*Fcurr);

    % If Neumann add to the gvec dt*(theta*NBCnext + (1-theta)*NBCcurr) 
    % Apply Dirichlet BCs in normal way 
    [gM,gV] = boundaryConditions(BC0,BC0type,BC1,BC1type,gM,gV,mesh,...
                                    time,dt,theta);

    % Solve final matrix system 
    cnext = gM \ gV;

    % Set new cnext to ccurr
    ccurr = cnext;

    % Add ccurr into the solution matrix
    solution(1:mesh.ngn,tstep+1) = ccurr;

    % Iterate the time step
    time = time + dt;
    % Add the time to tvec
    tvec(tstep+1,1) = time;

end

% Surface plot of the solution
surf(mesh.nvec,tvec,(solution)');
xmin = min(mesh.nvec,[],'all');
xmax = max(mesh.nvec,[],'all');
hold on
title(strcat(titles,"-t[0,",num2str(time),"]-x[",...
    num2str(xmin),",",num2str(xmax),"]"));
xlabel("x");
ylabel("time");
zlabel("Solution");
hold off

end

