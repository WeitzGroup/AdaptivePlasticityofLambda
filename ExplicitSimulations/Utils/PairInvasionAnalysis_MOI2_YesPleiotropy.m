% This script reads in the resident-mutant dynamics time series, calculates
% the exponential growthrate of the mutant frequency and generates a
% heatmap of this growthrate as a function of mutant strategy.

clear all;
close all;

files = dir('../Data/YesPleiotropy_Mutant_Period=T=3zbase/*Mutant*.mat');

datatable = array2table(zeros(0,3),'VariableNames',{'z1','z2','growth_rate'});

poolobj = parpool;

parfor i = 1:length(files)

   
    S = load([files(i).folder '\' files(i).name]);
    [r,~] = MutantGrowthRate(S.T1,S.Y1);

    datatable(i,:) = array2table([S.Params.z(2,1) S.Params.z(2,2) r]);
end

delete(poolobj);

%% Create figure
f = figure;
t = tiledlayout(1,3);

%% Heatmap of growth rate
M = unstack(datatable,'growth_rate','z2');

matrix = table2array(M(:,2:end)); %% each row corresponds to a value of z1 and each column corresponds to a value of z2.

nexttile(2);
imagesc(0:0.01:1,0:.01:1,matrix');
axis xy;
hold on;
[X,Y] = meshgrid(.01:.01:1,.01:.01:1);
contour(X,Y,matrix',[0 0],'-k');
scatter(0,1,'filled','MarkerFaceColor','k');
xlabel('$z_1$','Interpreter','latex');
ylabel('$z_2$','Interpreter','latex');
c = colorbar;
c.Label.String = 'Growth rate (hr$^{-1}$)';
c.Label.Interpreter = 'latex';
c.Label.Rotation = -90;
title('Mutant growth rate for T = 3 hr','Interpreter','latex');


%% Single trajectory.
S_resident = load('Data\NoPleiotropy_T=3\NoPleiotropy_Resident_Period=3,z1=0.00,z2=1.00.mat');
S_mutant = load('Data\NoPleiotropy_T=3\NoPleiotropy_Mutant_Period=3,z1=0.50,z2=0.50.mat');

temp_index = S_resident.T > (S_resident.T(end) - 3*S_resident.params.T);
T = S_resident.T(temp_index);
Y = S_resident.Y(temp_index,:);

T = [T;T(end)+S_mutant.T1];
Y = [Y;S_mutant.Y1];

T = T(T < T(1)+9*S_resident.params.T);
Y = Y(T < T(1)+9*S_resident.params.T,:);

% obtain densities of resident and mutant genomes
resident = Y(:,2)+ Y(:,4) + 0.5*(Y(:,5)+Y(:,6)) + Y(:,8) + Y(:,10);
mutant = Y(:,3)+ Y(:,7) + 0.5*(Y(:,5)+Y(:,6)) + Y(:,9) + Y(:,11);

% obtain mutant fraction
mutantfraction = mutant./(resident+mutant);
mutantlogfraction = log(mutant)./log(resident.*mutant);
residentlogfraction = 1-mutantlogfraction;
nexttile(1);
semilogy(T,1-mutantfraction,'-k','LineWidth',1.5,'DisplayName','Resident fraction'); hold on;
semilogy(T,mutantfraction,'--k','LineWidth',1.5,'DisplayName','Mutant fraction');
xline(S_resident.T(end));
legend('Resident fraction','Mutant fraction','Interpreter','latex');
xlabel('Time (hr)','Interpreter','latex');
ylabel('Population log fraction','Interpreter','latex');


%% Trait trajectory.
nexttile(3);
load('MultispeciesDynamics_n=25,Period=3,z1=0.00,z2=1.00.mat');
plot(T,mean_phi(:,1),'LineWidth',1,'DisplayName','$\langle \phi_1 \rangle$');
hold on
plot(T,mean_phi(:,2),'LineWidth',1,'DisplayName','$\langle \phi_2 \rangle$');
legend('Interpreter','latex')
xlabel('Time (hr)','Interpreter','latex');
ylabel('Mean trait value','Interpreter','latex');