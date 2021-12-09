% Steps innit bludren

% Initialise mesh
xmin = 0;
xmax = 1;
ne = 10;
mesh = OneDimLinearMeshGen(xmin,xmax,ne);

% Initialise theta, detla t and total time 
theta = 0;
dt = 0.1;
N = 10;

% Define material coefficients (D, lambda, f)
% These are currently being defined as function in x and t
% If you wish for a constant value, D = @(x,t) 5; for example
% For a varying value, use the format D = @(x,t) (x^2)*t + 3; for example
D = @(x,t) 0;
lambda = @(x,t) 0;
f = @(x,t) 0;

% Define initial conditions as a function of x
% Define constant and varying values as described above
ICs = @(x) 0;
% Create the first ccurr vector corresponding to the initial conditions
ccurr = initialConditions(ICs,mesh);
% Define cnext as a column vector of zeros
cnext = zeros(mesh.ngn,1);

% Define the boundary condition values as functions of time
BC0 = @(t) 0;
BC1 = @(t) 0;
% Declare the types of the two boundaries ("Neumann" or "Dirichlet")
BC0type = "Neumann";
BC1type = "Dirichlet";

% Initialise the time variable
t = 0;

% Loop over all the time steps 
for tstep = 1:N

    % Create the global mass matrix
    M = globalMassMatrix(mesh);

    % Create the global stiffness matrix 
    K = globalStiffnessMatrix(D,lambda,mesh,t);

    % Create global matrix ( M + theta * dt * K )
    gM = M + (theta*dt*K);

    % Create the matrix to multiply prev cvec ( M - (1-theta)*dt*K )
    P = M - (1-theta)*dt*K;

    % Multiply current sulition (M-(1-theta)*dt*K)*ccurr store in global
    % vector! )
    gV = P * ccurr;
    
    % Create the global vector ( dt * ( theta * Fnext + (1-theta)*Fcurr ) )
    Fcurr = globalSourceVector(f,mesh,t);
    Fnext = globalSourceVector(f,mesh,t+dt);
    gV = gV + dt*(theta*Fnext + (1-theta)*Fcurr);

    % If Neumann add to the gvec dt*(theta*NBCnext + (1-theta)*NBCcurr) 
    % Apply Dirichlet BCs in normal way 
    [gV,gM] = boundaryConditions(BC0,BC0type,BC1,BC1type,gM,gV,mesh,t,dt,theta);

    % Solve final matrix system 
    cnext = gM \ gV;

    % Set new cnext to ccurr
    ccurr = cnext;

    % Iterate the time step
    t = t + dt;

end

% Plot the solution

