%% PairWiseInvasion dynamics in the plastic and pleiotropic case at MOI = 2.

% This function accepts the oscillation period and predicted coESS strategy
% for the resident as parameters. Then computes the steady state cycle for
% the resident. Then converts 1% of resident virus to mutant and computes
% trajectory from there for the first 20 cycles. 

% The trajectory up until the resident reaches steady state is stored in
% 'Data/YesPleiotropy/Resident_Period=<Period>,Xi=<ResidentStrategy(1)>,Nu=<ResidentStrategy(2)>.mat'.
% The trajectory after the addition of the mutant is stored in 
% 'Data/YesPleiotropy/Mutant_Period=<period>,z1=<z1_mutant>,z2=<z2_mutant>.mat'.

% Note that this function uses the parfor loop and requires the parallel
% computation toolbox.

%% Inputs:
% ResidentStrategy: 1x2 vector. First element is Xi and second is Nu
% Period: scalar. Period oscillation in hours.

%Author: Tapan Goel
%Date created: 4/22/2025
%Date modified: 6/19/2025

function PairInvasionDynamics_YesPleiotropy(ResidentStrategy, period)

%% Create subfolder to store the data in, based on period
if ~isfolder(sprintf('../Data/YesPleiotropy/Period=%d',period))
    mkdir(sprintf("../Data/YesPleiotropy/Period=%d",period));
end

%% Fixed life history parameters
params.r = 1; %cell growth rate in per hour
params.kappa = 1e-9; %inverse carrying capacity in mL/cells
params.a = 1e-7; % adsorption rate in mL/hr
params.d = 0.9; %per capita cell death rate in per hour
params.d_V = 0.9; %per capita viral decay rate in per hour.
params.alpha = 1e-2; %per capita induction rate in per hour
params.mu = 2; %per capita transition rate in per hour
params.B0 = 80; %burst size baseline
params.cost = 1e-4; %cost of pleiotropy
params.numstrains = 2;

%% External forcing
params.T = period; %duration of cycle in hours.
params.theta_max = 0.1/params.kappa; %max host input rate;
params.theta_min = 0; % min host input rate;
params.good = 0.1; %fraction of cycle time at max host input;
params.theta = @(t) 0.5*((params.theta_max - params.theta_min)*(square(t*2*pi/params.T,params.good*100)) + params.theta_max + params.theta_min);


%% Initial conditions
S0 = .1/params.kappa;
V0 = .1/params.kappa;
x0 = [S0 zeros(1,8) V0 0]; %initial condition with only resident

%% Define viral strategy
z1_resident = 1/(1+ResidentStrategy(1));
z2_resident = 2^ResidentStrategy(2)/(2^ResidentStrategy(2)+ResidentStrategy(1));
params.z = [z1_resident z2_resident; ...
            z1_resident z2_resident];

params.phi1 = params.z(:,1);
params.phi2 = 0.5*(repmat(params.z(:,2),1,2) + repmat(params.z(:,2),1,2)'); %phi2(i,j) = .5*(z_i^(2) + z_j^(2));

params.B = params.B0*(1- params.cost*([ResidentStrategy(1)^2 + ResidentStrategy(2)^2;ResidentStrategy(1)^2 + ResidentStrategy(2)^2]));
%% Tolerance of simulations
params.tolerance = 1e-3;
%% Solve ODEs to resident steady state using the YesPleiotropy
tic;
[T,Y] = ode45(@(t,y)ODE_SELV_MOI2_YesPleiotropy(t,y,params), [0:0.1:params.T], x0);  %%Run ode dynamics for 1 cycle.

iter1 = 0;

while ~(sum( (abs(Y(end,:)-x0) < params.tolerance) | (abs(Y(end,:)-x0)./Y(end,:) < 1e-6) ) == length(x0)) & iter1 < 5000
    x0 = Y(end,:);
    [t,y] = ode45(@(t,y)ODE_SELV_MOI2_YesPleiotropy(t,y,params), T(end) + [0:0.1:params.T], x0);
    T = [T;t(2:end)];
    Y = [Y;y(2:end,:)];
    iter1 = iter1+1;
   
end
toc;
residentdynamics_filename = sprintf('../Data/YesPleiotropy/Period=%d/Resident_Xi=%.2f,Nu=%.2f.mat',period,ResidentStrategy(1),ResidentStrategy(2));
save(residentdynamics_filename,"T","Y","params","iter1");

%% Add mutant
x0 = Y(end,:);
x0(8)= 0.99*Y(end,8); %reduce resident lysogens
x0(9) = 0.01*Y(end,8); %add mutant lysogens
x0(10) = 0.99*Y(end,10); %reduce resident virions
x0(11) = 0.01*Y(end,10); %add mutant virions

tic;
%% Define range of mutant strategies:
Z = 0:0.01:1;

poolobj = parpool;

parfor ii = 1:length(Z)*length(Z)
    
    [i,j] = ind2sub([length(Z) length(Z)],ii);
    X0 = x0;
    Params = params;

    %% Reset strategies
    z1_mutant = Z(i);
    z2_mutant = Z(j);

    MutantStrategy = [1/Z(i)-1 log2((1/Z(i)-1)/(1/Z(j)-1))]; %% Mutant strategy as (Xi,Nu);
   
    Params.z = [z1_resident z2_resident; ...
                z1_mutant z2_mutant];
    Params.phi1 = Params.z(:,1);
    Params.phi2 = 0.5*(repmat(Params.z(:,2),1,2) + repmat(Params.z(:,2),1,2)'); %phi2(i,j) = .5*(z_i^(2) + z_j^(2));

    Params.B = Params.B0*(1- Params.cost*([ResidentStrategy(1)^2 + ResidentStrategy(2)^2;MutantStrategy(1)^2 + MutantStrategy(2)^2]));
    Params.B(Params.B<0 | isnan(Params.B)) = 0; % Reset any negative burst size to zero
    %% Solve mutant to 20 cycles
    [T1,Y1] = ode45(@(t,y)ODE_SELV_MOI2_YesPleiotropy(t,y,Params), [0:0.1:Params.T], X0);  %%Run ode dynamics for 1 cycle.
    
    iter = 0;
    %x = [];
    while iter < 20
        X0 = Y1(end,:);
        [t,y] = ode45(@(t,y)ODE_SELV_MOI2_YesPleiotropy(t,y,Params), T1(end) + [0:0.1:Params.T], X0);
        T1 = [T1;t(2:end)];
        Y1 = [Y1;y(2:end,:)];
        iter = iter+1;
        
    end
    
    mutantdynamics_filename = sprintf('../Data/YesPleiotropy/Period=%d/Mutant_z1=%.2f,z2=%.2f.mat',period,z1_mutant,z2_mutant);
    parsave(mutantdynamics_filename,Params.T, MutantStrategy, T1, Y1, Params);

end
toc
delete(poolobj);
end