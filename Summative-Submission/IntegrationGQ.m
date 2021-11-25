function [out] = IntegrationGQ(const,coeffs)
%INTEGRATIONGQ Calculates the definite integral of an equation
%   Function takes in the constant infront of equation, and an array
%   consisting of all the variable coefficients in the equation. This will
%   return then numerical value of the integral in the range -1 <= x <= 1. 
%   EXAMPLE, general function in x:
%
%           to integrate f(x) = x^3 + 2x^2 - 4x + 12
% 
%           call the function with the following inputs:
%           const = 1 and coeffs = [12,-4,2,1]
%
%   EXAMPLE 2, in context of the FEM solver:
%
%           to integrate the const Diffusion term ie Int(00) = D/2J
%           
%           call the function with the following inputs:
%           const = (1/2J) and coeffs = [D]
%
%   EXAMPLE 3, for a spacially varying diffusion term:
%
%           to integrate the const Diffusion term ie Int(00) = (2x/2J)
%           
%           call the function with the following inputs:
%           const = (1/2J) and coeffs = [0,2]
%
%   The function will automatically determine and generate the GQ scheme.
%   If not clear: the coeff array starts from the coefficient of x^0
%   onwards, ie the second term is the coefficient of the x term etc



end

