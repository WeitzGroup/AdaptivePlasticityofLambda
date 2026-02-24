function dydt = ODE_SELV_MOI2(t,y,params)

% ODE_SELV_MOI2: solves the SELV ODEs for a one host, multi-virus systems with max MOI of 2.
% Note: function returns an error if max MOI is greater than 2.

% Inputs:
% t - (double) time
% y - (vector of doubles) system state at time t [S, E1, E2, L, V]; E1 is
%      vector of singly infected states and E2 is vector of doubly infected
%      states, L is vector of lysogen states.
% params - (struct) contains life history parameters, voting rules for
%       mixed infections and details of the host input function. 
% params.numstrains - (scalar) number of viral strains in system
% params.r - (scalar) cell growth rate
% params.kappa - (scalar) inverse carrying capacity
% params.d - (scalar) cell death rate
% params.a - (scalar) adsorption rate
% params.b - (scalar) infection efficiency
% params.d_V - (scalar) viral death rate
% params.alpha - (scalar) lysogen induction rate
% params.mu - (scalar) commitment rate
% params.cost - (scalar) cost of pleiotropy if there is pleiotropy.
% params.B - (numstrainsx1 vector) burst size B(i) is the burst size for
%               strain i. B(i) = B0(1-cost*(xi^2+nu^2))
% params.phi1 - (numstrainsx1 vector) phi1(i) is the probability of lysogeny for singly infected class E1_i
% params.phi2 - (numstrainsxnumstrains matrix) phi2(i,j) is the probability
%               of lysogeny of doubly infected class E2_(i,j)
%       

% Outputs:
% dydt = (vector of doubles) derivative of system state at time t d[S, E1,
% E2, L, V]/dt

% Author: Tapan Goel
% Date: 2026-02-23


    
    y = y(:); %make sure y is a column vector

    n = params.numstrains; %extract number of viral strains in the population
    
    %% abort if max MOI is greater than 2
    [~,c] = size(params.z);
    if c > 2
        error('max MOI set to be greater than 2. Check params.z');
    end
    %% abort if y vector contains elements more than those needed for incorporating max MOI 2
    if length(y) ~= n*n+3*n+1
        error('state vector array size incompatible with MOI = 2 model');
    end
    
    %% isolate all vectors as column vectors or matrices
    S = y(1);
    E1 = y(2:n+1);
    E2 = y(n+2:n+n*n+1);
    L = y(n*n+n+2:n*n+2*n+1);
    V = y(n*n+2*n+2:n*n+3*n+1);
    E2 = reshape(E2,[n,n]); %this converts the vector into columns. so E2(1:n) becomes the first column, E2(n+1:2*n) becomes the second column and so on;

    N = sum(y(1:n*n+2*n+1)); %total number of cells
    Vtot = sum(V); %total number of free virus particles
    
        
    dS = params.theta(t) + params.r*S*(1-params.kappa*N) - (params.d + params.a*1*Vtot)*S;

    dE1 = params.a*1*S*V - (params.d + params.mu + params.a*1*Vtot)*E1;

    dE2 = params.a*1*E1*V' - (params.d+params.mu)*E2;

    dL = params.r*(1-params.kappa*N)*L - (params.d+params.alpha)*L + params.mu*params.phi1.*E1 + 0.5*params.mu*(sum(params.phi2.*E2,2) + sum((params.phi2.*E2)',2)); %the .5 factor in th last term is to account for double counting in doing the sum E_ij+Eji

    dV = params.B.*params.alpha.*L - (params.d_V+params.a*N)*V + params.B.*params.mu.*(1-params.phi1).*E1 + 0.5*params.B.*params.mu.*( sum((1-params.phi2).*E2,2) + sum(((1-params.phi2).*E2)',2) );


    dydt = [dS;dE1;dE2(:);dL;dV]; %reconstruct dydt as a column vector. Note that dE2(:) creates a col vector [dE2(:,1); dE2(:,2) ....]
end

    