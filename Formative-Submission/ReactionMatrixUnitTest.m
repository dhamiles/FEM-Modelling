%% Test 1: test symmetry of the matrix
% % Test that this matrix is symmetric
tol = 1e-14;
lambda = -5; %reaction coefficient
eID=1; %element ID
msh = OneDimLinearMeshGen(0,1,10);

elemat = ReactionElemMatrix(lambda,eID,msh); 
assert(abs(elemat(1,2) - elemat(2,1)) <= tol)

%% Test 2: test 2 different elements of the same size produce same matrix
% % Test that for two elements of an equispaced mesh, as described in the
% % lectures, the element matrices calculated are the same
tol = 1e-14;
lambda = 7; %reaction coefficient
eID=1; %element ID
msh = OneDimLinearMeshGen(0,1,10);

elemat1 = ReactionElemMatrix(lambda,eID,msh);

eID=2; %element ID

elemat2 = ReactionElemMatrix(lambda,eID,msh);

diff = elemat1 - elemat2;
diffnorm = sum(sum(diff.*diff));
assert(abs(diffnorm) <= tol)

%% Test 3: test that one matrix is evaluted correctly
% % Test that element 1's matrix from the example in the lectures (3
% % elements)evaluates correctly
tol = 1e-14;
lambda = 3; %reaction coefficient
eID=1; %element ID
msh = OneDimLinearMeshGen(0,1,3);

elemat1 = ReactionElemMatrix(lambda,eID,msh);

elemat2 = [ 1/3 1/6; 1/6 1/3];
diff = elemat1 - elemat2; %calculate the difference between the two matrices
diffnorm = sum(sum(diff.*diff)); %calculates the total squared error between the matrices
assert(abs(diffnorm) <= tol)

%% Test 4: test an element in non equi-spaced mesh is correct
% % Evaluates that the last element in a non equi-spaced mesh is evaluated
% % correctly
tol = 1e-14;
lambda = -9; %reaction coefficient
eID=4; %element ID
msh = OneDimSimpleRefinedMeshGen(0,1,4);

elemat1 = ReactionElemMatrix(lambda,eID,msh);

elemat2 = [ -0.3750 -0.1875; -0.1875 -0.3750];
diff = elemat1 - elemat2; %calculate the difference between the two matrices
diffnorm = sum(sum(diff.*diff)); %calculates the total squared error between the matrices
assert(abs(diffnorm) <= tol)