function [mesh] = OneDimMeshGen(xmin,xmax,Ne,order)
%%This function generates a one dimensional, equispaced, linear finite
%%element mesh, with Ne number of elements, between the points at x
%%position xmin and xmax.

mesh.ne = Ne; %set number of elements
dx = (xmax - xmin)/Ne; %calculate element size
mesh.order = order; % save the lagrange order in the mesh

if order==1

    mesh.ngn = Ne+1; %set number of global nodes
    mesh.nvec = zeros(mesh.ngn,1); %create node vector to store values
    mesh.nvec = xmin:dx:xmax;
    
    %loop over elements & set the element properties
    for i=1:Ne

        %set spatial positions of nodes
        mesh.elem(i).x(1) = xmin + (i-1)*dx;
        mesh.elem(i).x(2) = xmin + i*dx ;

        %set global IDs of the nodes
        mesh.elem(i).n(1) = i;
        mesh.elem(i).n(2) = i+1;

        %set element Jacobian based on mapping to standard element
        mesh.elem(i).J = 0.5*dx; %assuming standard element of -1 to 1

    end

elseif order==2

    mesh.ngn = (2*Ne)+1; %set number of global nodes
    mesh.nvec = zeros(mesh.ngn,1); %create node vector to store values
    mesh.nvec = xmin:(dx/2):xmax;
    
    %loop over elements & set the element properties
    for i=1:Ne

        %set spatial positions of nodes
        mesh.elem(i).x(1) = xmin + (i-1)*dx;
        mesh.elem(i).x(2) = xmin + i*(dx/2);
        mesh.elem(i).x(3) = xmin + i*dx;

        %set global IDs of the nodes
        mesh.elem(i).n(1) = 1 + (i-1)*2;
        mesh.elem(i).n(2) = 2 + (i-1)*2;
        mesh.elem(i).n(3) = 3 + (i-1)*2;

        %set element Jacobian based on mapping to standard element
        mesh.elem(i).J = 0.5*dx; %assuming standard element of -1 to 1

    end

else
    error("Invalid order for Lagrange basis functions, please choose" + ...
        "1 (linear) or 2 (quadratic)");

end
