%% PairWiseInvasion dynamics in the plastic and pleiotropic case at MOI = 2.

% This function accepts the oscillation period and predicted coESS strategy
% for the resident as parameters. Then computes the steady state cycle for
% the resident. Then converts 1% of resident virus to mutant and computes
% trajectory from there for the first 20 cycles. 

% The trajectory up until the resident reaches steady state is stored in
% 'Data/YesPleiotropy/Period=<Period>/Resident_z1=<ResidentStrategy(1)>,z2=<ResidentStrategy(2)>.mat'.
% The trajectory after the addition of the mutant is stored in 
% 'Data/YesPleiotropy/Period=<Period>/Mutant_z1=<MutantStrategy(1)>,z2=<MutantStrategy(2)>.mat'.

% Also stores all parameters and all resident and mutant trajectories in 
% Data/YesPleiotropy/Period=<Period>/Resident_z1=<ResidentStrategy(1)>,z2=<ResidentStrategy(2)>_Summary.mat

% Note that this function uses the parfor loop and requires the parallel
% computation toolbox.

%% Inputs:
% ResidentStrategy: 1x2 vector. First element is z1 and second is z2
% Period: scalar. Period oscillation in hours.

%Author: Tapan Goel
%Date created: 4/22/2025
%Date modified: 6/19/2025

function [coESS, ResidentTrajectory, MutantTrajectories, GrowthRateTable] = PairInvasionDynamics(coESS, period, vals)

arguments
    coESS (1,2) double
    period double
    vals.Pleiotropy logical = false
    vals.Save logical  = false
    vals.MutantFraction double = 0.01
    vals.MutantStrategies (:,2) double = table2array(combinations(0:.01:1,0:.01:1))
end

%% If save flag is true, create folder structure and filename to to save data
% save file if save flag is true
if vals.Save
    if ~isfolder('Data')
        mkdir('Data')
    end
    if vals.Pleiotropy
        if ~isfolder('Data/YesPleiotropy')
            mkdir('Data/YesPleiotropy')
        end
        subfolder = 'Data/YesPleiotropy/';
    else
        if ~isfolder('Data/NoPleiotropy')
            mkdir('Data/NoPleiotropy')
        end
        subfolder = 'Data/NoPleiotropy/';
    end
    savefolder = sprintf('%sInvasion_Period=%.2f',subfolder,period);
    if ~isfolder(savefolder)
        mkdir(savefolder);
    end
    filename = sprintf('%s/Resident_z1=%.3f,z2=%.3f,MutantFrac=%.2f.mat',savefolder,coESS(1),coESS(2),vals.MutantFraction);
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
params.numstrains = 2; % number of strains for invasion analysis is 2.
%% External forcing
params.T = period; %duration of cycle in hours.
params.theta_max = 0.1/params.kappa; %max host input rate;
params.theta_min = 0; % min host input rate;
params.good = 0.1; %fraction of cycle time at max host input;
params.theta = @(t) 0.5*((params.theta_max - params.theta_min)*(square(t*2*pi/params.T,params.good*100)) + params.theta_max + params.theta_min);


%% Resident dynamics
% Initial conditions
S0 = .1/params.kappa;
V0 = .1/params.kappa;
x0 = [S0 zeros(1,8) V0 0]; %initial condition with only resident

% Viral strategy
params.z = [1;1]*coESS;
params.phi1 = params.z(:,1);
params.phi2 = 0.5*(repmat(params.z(:,2),1,2) + repmat(params.z(:,2),1,2)'); %phi2(i,j) = .5*(z_i^(2) + z_j^(2));

% Cost of pleiotropy
if vals.Pleiotropy
    params.cost = 1e-4;
    xi = 1./params.z(:,1)-1;
    nu = log2((1./params.z(:,1)-1)./(1./params.z(:,2)-1));
    params.B = max(0,params.B0*(1-params.cost*(xi.^2+nu.^2)));
else
    params.cost = 0;
    params.B = params.B0*ones(2,1);
end

% Tolerance of simulations
params.tolerance = 1e-3;

% Solve ODEs to resident steady state using the YesPleiotropy
tic;
[T,Y] = ode45(@(t,y)ODE_SELV_MOI2(t,y,params), 0:0.1:params.T, x0);  %%Run ode dynamics for 1 cycle.

iter1 = 0;

while ~(sum( (abs(Y(end,:)-x0) < params.tolerance) | (abs(Y(end,:)-x0)./Y(end,:) < 1e-6) ) == length(x0)) & iter1 < 5000
    x0 = Y(end,:);
    [t,y] = ode45(@(t,y)ODE_SELV_MOI2(t,y,params), T(end) + (0:0.1:params.T), x0);
    T = [T;t(2:end)];
    Y = [Y;y(2:end,:)];
    iter1 = iter1+1;
   
end
toc;

% Save resident trajectory
ResidentTrajectory = cell2table({params,T,Y},'VariableNames',{'Params','Time','Densities'});

%% Mutant Invasions

Z = vals.MutantStrategies; % mutant strategies

% Initialize mutant data storage
GrowthRateTable = array2table(zeros(length(Z),3),'VariableNames',{'Mutant_z1','Mutant_z2','GrowthRate'});
MutantTrajectories = cell(length(Z),5);

% Add mutants
x0 = Y(end,:);
x0(8)= (1-vals.MutantFraction)*Y(end,8); %reduce resident lysogens
x0(9) = vals.MutantFraction*Y(end,8); %add mutant lysogens
x0(10) = (1-vals.MutantFraction)*Y(end,10); %reduce resident virions
x0(11) = vals.MutantFraction*Y(end,10); %add mutant virions

poolobj = parpool;
tic;
parfor ii = 1:length(Z)
    
    X0 = x0; % set initial condition for invasion
    Params = params; % copy parameters before modifying to add mutant parameters.

    % Reset strategies
    Params.z = [coESS(1) coESS(2); ...
                Z(ii,1) Z(ii,2)];

    % Cost of pleiotropy
    if vals.Pleiotropy
        Params.cost = 1e-4;
        xi = 1./Params.z(:,1)-1;
        nu = log2((1./Params.z(:,1)-1)./(1./Params.z(:,2)-1));
        Params.B = max(0,Params.B0*(1-Params.cost*(xi.^2+nu.^2)));
    else
        Params.cost = 0;
        Params.B = Params.B0*ones(2,1);
    end
   
    Params.phi1 = Params.z(:,1);
    Params.phi2 = 0.5*(repmat(Params.z(:,2),1,2) + repmat(Params.z(:,2),1,2)'); %phi2(i,j) = .5*(z_i^(2) + z_j^(2));

    % Solve invasion dynamics to 20 cycles
    [T1,Y1] = ode45(@(t,y)ODE_SELV_MOI2(t,y,Params), 0:0.1:Params.T, X0);  %%Run ode dynamics for 1 cycle.
    
    iter = 0;
    while iter < 20
        X0 = Y1(end,:);
        [t,y] = ode45(@(t,y)ODE_SELV_MOI2(t,y,Params), T1(end) + (0:0.1:Params.T), X0);
        T1 = [T1;t(2:end)];
        Y1 = [Y1;y(2:end,:)];
        iter = iter+1;   
    end
    
    % Calculate and store mutant growth rate
    [r,~] = MutantGrowthRate(T1,Y1);
    GrowthRateTable(ii,:) = array2table([Z(ii,1) Z(ii,2) r],'VariableName',{'Mutant_z1','Mutant_z2','GrowthRate'});
  
    % Store mutant trajectory
    MutantTrajectories(ii,:) = {Z(ii,1), Z(ii,2), Params.B, T1, Y1};
    % clear variables
    T1 = []; 
    Y1 = [];
        
end
toc
delete(poolobj);

MutantTrajectories = cell2table(MutantTrajectories,'VariableNames',{'Mutant_z1', 'Mutant_z2','BurstSizes','Time','Densities'});


if vals.Save
    save(filename,"ResidentTrajectory","MutantTrajectories","GrowthRateTable");
end

end