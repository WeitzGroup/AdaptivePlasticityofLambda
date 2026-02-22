%% Multispecies dynamics in the plastic and pleiotropic case at MOI = 2.

% This function accepts the oscillation period, predicted coESS strategy
% for the resident, number of cycles and number of strains as parameters. It then computes trajectory and saves the output to a file

% The non-resident strategies are chosen at random using LHS sampling

%Author: Tapan Goel
%Date created: 6/19/2025

function MultiSpeciesTrajectories_YesPleiotropy(coESS, period, numstrains, numcycles)
%% Fixed life history parameters
params.r = 1; %cell growth rate in per hour
params.kappa = 1e-9; %inverse carrying capacity in mL/cells
params.a = 1e-7; % adsorption rate in mL/hr
params.d = 0.9; %per capita cell death rate in per hour
params.d_V = 0.9; %per capita viral decay rate in per hour.
params.alpha = 1e-2; %per capita induction rate in per hour
params.mu = 2; %per capita transition rate in per hour
params.B0 = 80; %burst size.
params.numstrains = numstrains;
params.cost = 1e-4;

%% External forcing
params.T = period; %duration of cycle in hours.
params.theta_max = 0.1/params.kappa; %max host input rate;
params.theta_min = 0; % min host input rate;
params.good = 0.1; %fraction of cycle time at max host input;
params.theta = @(t) 0.5*((params.theta_max - params.theta_min)*(square(t*2*pi/params.T,params.good*100)) + params.theta_max + params.theta_min);


%% Initialize viral strategies:
if numstrains == 1
    params.z = [coESS(1) coESS(2)];
    xi = 1/params.z(1)-1;
    nu = log2((1/params.z(1)-1)/(1/params.z(2)-1));
    params.B = max(0,params.B0*(1-params.cost*(xi^2+nu^2)));
    S0 = 0.1/params.kappa;
    V0 = 0.1/params.kappa;
    x0 = [S0 zeros(1,3) V0];
elseif numstrains == 2
    params.z = [coESS(1) coESS(1);...
                1-coESS(1) 1-coESS(2)];
    xi = 1./params.z(:,1)-1;
    nu = log2((1./params.z(:,1)-1)./(1./params.z(:,2)-1));
    params.B = max(0,params.B0*(1-params.cost*(xi.^2+nu.^2)));
    S0 = 0.1/params.kappa;
    V0 = (0.1/params.kappa)/numstrains*ones(1,numstrains);
    x0 = [S0 zeros(1,8) V0];
else
    params.z = [coESS(1) coESS(2);...
                1-coESS(1) 1-coESS(2);...
                lhsdesign(numstrains-2,2)];
    xi = 1./params.z(:,1)-1;
    nu = log2((1./params.z(:,1)-1)./(1./params.z(:,2)-1));
    params.B = max(0,params.B0*(1-params.cost*(xi.^2+nu.^2)));
    S0 = 0.1/params.kappa;
    V0 = (0.1/params.kappa)/numstrains*ones(1,numstrains);
    x0 = [S0 zeros(1,2*numstrains+numstrains*numstrains) V0];
end

params.phi1 = params.z(:,1);
params.phi2 = 0.5*(repmat(params.z(:,2),1,numstrains) + repmat(params.z(:,2),1,numstrains)'); %phi2(i,j) = .5*(z_i^(2) + z_j^(2));

%% Tolerance of simulations
params.tolerance = 1e-3;
%% Solve ODEs to resident steady state using the NoPleiotropy
tic;
[T,Y] = ode45(@(t,y)ODE_SELV_MOI2_YesPleiotropy(t,y,params), [0:0.1:params.T], x0);  %%Run ode dynamics for 1 cycle.

iter1 = 0;

while ~(sum( (abs(Y(end,:)-x0) < params.tolerance) | (abs(Y(end,:)-x0)./Y(end,:) < 1e-6) ) == length(x0)) & iter1 < numcycles
    x0 = Y(end,:);
    [t,y] = ode45(@(t,y)ODE_SELV_MOI2_YesPleiotropy(t,y,params), T(end) + [0:0.1:params.T], x0);
    T = [T;t(2:end)];
    Y = [Y;y(2:end,:)];
    iter1 = iter1+1;
   
end
toc;

S = Y(:,1);
E1 = Y(:,2:params.numstrains+1);
E2 = Y(:,params.numstrains+2:params.numstrains+params.numstrains*params.numstrains+1);
L = Y(:,params.numstrains*params.numstrains+params.numstrains+2:params.numstrains*params.numstrains+2*params.numstrains+1);
V = Y(:,params.numstrains*params.numstrains+2*params.numstrains+2:params.numstrains*params.numstrains+3*params.numstrains+1);
E2 = reshape(E2,[],params.numstrains,params.numstrains); %this converts the vector into columns. so E2(1:params.numstrains) becomes the first column, E2(params.numstrains+1:2*params.numstrains) becomes the second column and so on;

clear Y;
% calculate population fraction
popfractions = E1 + L + V + sum(E2,3) + sum(permute(E2,[1 3 2]),3);

popfractions = popfractions./sum(popfractions,2);

% rearrange pop fractions and z values by dominance
[~,indices] = sort(popfractions(end,:));
popfractions = popfractions(:,indices);
strategies_ordered = params.z(indices,:);

mean_phi = popfractions*strategies_ordered;
std_phi = sqrt(popfractions*(strategies_ordered.^2)-mean_phi.^2);

filename = sprintf('../Data/MultispeciesDynamics_YesPleiotropy_n=%d,Period=%d,z1=%.2f,z2=%.2f.mat',numstrains, period,coESS(1),coESS(2));

save(filename,"T","S","E1","E2","L","V","popfractions","strategies_ordered","mean_phi","params","iter1","std_phi");


end




