%% Test 1: Test 2 Dirichlet Boundary Conditions
% % This test will confirm whether the function correctly returns the
% % modified global matrix and vector for the condition

% Set the Dirichlet values at the 2 boundaries
bc0 = 3;
bc1 = 0;
t = "Dirichlet";

tol = 1e-14; % Assign a tolerance
mesh = OneDimLinearMeshGen(0,1,4); % Create the mesh

% Create global matrix/vector of 1's in what would be non-zero positions
gmatrix = [1,1,0,0,0;
          1,1,1,0,0;
          0,1,1,1,0;
          0,0,1,1,1;
          0,0,0,1,1];

gvector = [1;
           1;
           1;
           1;
           1];
       
% Apply boundary conditions
[gmatrix,gvector] = BoundaryConditions(bc0,t,bc1,t,gmatrix,mesh,gvector);

% Create the expected matrix/vector
expmat = [1,0,0,0,0;
          1,1,1,0,0;
          0,1,1,1,0;
          0,0,1,1,1;
          0,0,0,0,1];

expvec = [bc0;
          1;
          1;
          1;
          bc1];
      
% Calculate differences and check against tolerance
diffv = gvector - expvec;
diffnormv = sum(sum(diffv.*diffv));
diffm = gmatrix - expmat;
diffnormm = sum(sum(diffm.*diffm));

% Assert the difference in both vector and matrix is below tolerance
assert((abs(diffnormv) <= tol) && (abs(diffnormm) <= tol))
      
%% Test 2: Test 2 Neumann Boundary Conditions
% % This test will confirm whether the function correctly returns the
% % modified global matrix and vector for the condition

% Set the Neumann values at the 2 boundaries
bc0 = 3;
bc1 = 0;
t = "Neumann";

tol = 1e-14; % Assign a tolerance
mesh = OneDimLinearMeshGen(0,1,4); % Create the mesh

% Create global matrix/vector of 1's in what would be non-zero positions
gmatrix = [1,1,0,0,0;
          1,1,1,0,0;
          0,1,1,1,0;
          0,0,1,1,1;
          0,0,0,1,1];

gvector = [1;
           1;
           1;
           1;
           1];
       
% Apply boundary conditions
[gmatrix,gvector] = BoundaryConditions(bc0,t,bc1,t,gmatrix,mesh,gvector);

% Create the expected matrix/vector
expmat = [1,1,0,0,0;
          1,1,1,0,0;
          0,1,1,1,0;
          0,0,1,1,1;
          0,0,0,1,1];

expvec = [1-bc0;
          1;
          1;
          1;
          1+bc1];
      
% Calculate differences and check against tolerance
diffv = gvector - expvec;
diffnormv = sum(sum(diffv.*diffv));
diffm = gmatrix - expmat;
diffnormm = sum(sum(diffm.*diffm));

% Assert the difference in both vector and matrix is below tolerance
assert((abs(diffnormv) <= tol) && (abs(diffnormm) <= tol))

%% Test 3: Test 1 Neuman and 1 Dirichlet
% % This test will confirm whether the function correctly returns the
% % modified global matrix and vector for the condition

% Set the Neumann values at the 2 boundaries
bc0 = 3;
bc1 = 0;
t0 = "Neumann";
t1 = "Dirichlet";

tol = 1e-14; % Assign a tolerance
mesh = OneDimLinearMeshGen(0,1,4); % Create the mesh

% Create global matrix/vector of 1's in what would be non-zero positions
gmatrix = [1,1,0,0,0;
          1,1,1,0,0;
          0,1,1,1,0;
          0,0,1,1,1;
          0,0,0,1,1];

gvector = [1;
           1;
           1;
           1;
           1];
       
% Apply boundary conditions
[gmatrix,gvector] = BoundaryConditions(bc0,t0,bc1,t1,gmatrix,mesh,gvector);

% Create the expected matrix/vector
expmat = [1,1,0,0,0;
          1,1,1,0,0;
          0,1,1,1,0;
          0,0,1,1,1;
          0,0,0,0,1];

expvec = [1-bc0;
          1;
          1;
          1;
          bc1];
      
% Calculate differences and check against tolerance
diffv = gvector - expvec;
diffnormv = sum(sum(diffv.*diffv));
diffm = gmatrix - expmat;
diffnormm = sum(sum(diffm.*diffm));

% Assert the difference in both vector and matrix is below tolerance
assert((abs(diffnormv) <= tol) && (abs(diffnormm) <= tol))


