%% Multispecies dynamics when MOI = 2.

% This function accepts the oscillation period, predicted coESS strategy
% for the resident, number of cycles and number of strains as parameters.
% It then computes trajectory, saves the output to a file and returns the
% matrix of strategies (best performing to worst), population fractions 
% (from worst to best performing) and time series of
% mean and std of traitvalues

% The non-resident strategies are chosen at random using LHS sampling

%Author: Tapan Goel
%Date created: 2/23/2026

function [coESS, OrderedStrategies, PopFractions, Time, TraitMean, TraitStdev] = MultiSpeciesTrajectories(coESS, period, numstrains, numcycles, vals)

arguments
    coESS (1,2) double
    period double
    numstrains double
    numcycles uint16
    vals.Pleiotropy logical = false
    vals.Save logical  = false
end

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



%% External forcing
params.T = period; %duration of cycle in hours.
params.theta_max = 0.1/params.kappa; %max host input rate;
params.theta_min = 0; % min host input rate;
params.good = 0.1; %fraction of cycle time at max host input;
params.theta = @(t) 0.5*((params.theta_max - params.theta_min)*(square(t*2*pi/params.T,params.good*100)) + params.theta_max + params.theta_min);

rng(0); % so that sim results are reproducible

%% Set initial host and virus density.
S0 = 0.1/params.kappa;
V0 = (0.1/params.kappa)/numstrains*ones(1,numstrains);
x0 = [S0 zeros(1,2*numstrains+numstrains*numstrains) V0];

%% Initialize viral strategies
if numstrains == 1
    params.z = [coESS(1) coESS(2)];
elseif numstrains == 2
    params.z = [coESS(1) coESS(2);...
                1-coESS(1) 1-coESS(2)];
else
    params.z = [coESS(1) coESS(2);...
                1-coESS(1) 1-coESS(2);...
                lhsdesign(numstrains-2,2)];
end

% Assign pleitropy costs and burst size accordingly
if vals.Pleiotropy
    params.cost = 1e-4;
    xi = 1./params.z(:,1)-1;
    nu = log2((1./params.z(:,1)-1)./(1./params.z(:,2)-1));
    params.B = max(0,params.B0*(1-params.cost*(xi.^2+nu.^2)));
else
    params.cost = 0;
    params.B = params.B0*ones(numstrains,1);
end

params.phi1 = params.z(:,1);
params.phi2 = 0.5*(repmat(params.z(:,2),1,numstrains) + repmat(params.z(:,2),1,numstrains)'); %phi2(i,j) = .5*(z_i^(2) + z_j^(2));

%% Tolerance of simulations
params.tolerance = 1e-3;
%% Solve ODEs to resident steady state using the NoPleiotropy
tic;
[T,Y] = ode45(@(t,y)ODE_SELV_MOI2(t,y,params), 0:0.1:params.T, x0);  %%Run ode dynamics for 1 cycle.

iter1 = 0;

while ~(sum( (abs(Y(end,:)-x0) < params.tolerance) | (abs(Y(end,:)-x0)./Y(end,:) < 1e-6) ) == length(x0)) & iter1 < numcycles
    x0 = Y(end,:);
    [t,y] = ode45(@(t,y)ODE_SELV_MOI2(t,y,params), T(end) + (0:0.1:params.T), x0);
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
PopFractions = E1 + L + V + sum(E2,3) + sum(permute(E2,[1 3 2]),3);

PopFractions = PopFractions./sum(PopFractions,2);

% rearrange pop fractions and z values by dominance
TraitMean = PopFractions*params.z;
TraitStdev = sqrt(PopFractions*(params.z.*params.z)-TraitMean.^2);

[~,indices] = sort(PopFractions(end,:));
PopFractions = PopFractions(:,indices);
OrderedStrategies = params.z(indices,:);

Time = T;

% save file if save flag is true
if vals.Save
    if ~isfolder('Data')
        mkdir('Data')
    end
    if vals.Pleiotropy
        if ~isfolder('Data/YesPleiotropy')
            mkdir('Data/YesPleiotropy')
        end
        filename = sprintf('Data/YesPleiotropy/MultiSpeciesDynamics_Period=%.2f,n=%d,z1=%.3f,z2=%.3f.mat',period,numstrains,coESS(1),coESS(2));
    else
        if ~isfolder('Data/NoPleiotropy')
            mkdir('Data/NoPleiotropy')
        end
        filename = sprintf('Data/NoPleiotropy/MultiSpeciesDynamics_Period=%.2f,n=%d,z1=%.3f,z2=%.3f.mat',period,numstrains,coESS(1),coESS(2));
    end
    save(filename,"T","PopFractions","OrderedStrategies","TraitMean","TraitStdev","params");
end

end




