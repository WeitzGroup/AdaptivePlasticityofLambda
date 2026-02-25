% Ths script generates figure S2 of the manuscript which creates the
% 1. sample resident attractor dynamics
% 2. sample invasion dynamics
% 3. Heatmap of invasion fitness
% for both, pleiotropic and non-pleiotropic cases. 
% The data for generating these plots is stored in 'Data/FigureData.mat'

%Date: July 16, 2025
%Author: Tapan Goel

clear;
close all;

%% Generate figure
addpath('Utils\');
%load('Data/FigureData.mat'); %load file that contains data for the figure.
f = figure('Position',[160,190,1904,935]);
t = tiledlayout(2,3);


%% No Pleiotropy case %%%%%%%%%%%%
coESS = [0 1];

% sample resident attractor -- no pleiotropy
if ~isfile('Data/NoPleiotropy/Invasion_Period=3.00/Resident_z1=0.000,z2=1.000,MutantFrac=0.50.mat')
    [coESS, ResidentTrajectory, MutantTrajectories, GrowthRateTable] = PairInvasionDynamics([0 1], 3, MutantStrategies= [.5 .5],...
                                                                        MutantFraction = 0.5, Save = true, Pleiotropy = false);
end
load('Data/NoPleiotropy/Invasion_Period=3.00/Resident_z1=0.000,z2=1.000,MutantFrac=0.50.mat');

%Sample resident attractor -- no pleiotropy
nexttile(1);
T = [ResidentTrajectory.Time{:}];
Y = [ResidentTrajectory.Densities{:}];

semilogy(T,Y(:,1),'LineWidth',1,'Color','b','DisplayName','S'); hold on;
semilogy(T,Y(:,2),'LineWidth',1,'Color','r','LineStyle','-','DisplayName','E$_1$');
semilogy(T,Y(:,4),'LineWidth',1,'Color','r','LineStyle','--','DisplayName','E$_2$');
semilogy(T,Y(:,8),'LineWidth',1,'Color',[165 42 42]/255,'LineStyle','-','DisplayName','L');
semilogy(T,Y(:,10),'LineWidth',1,'Color','g','LineStyle','-','DisplayName','V');
set(gca,'YMinorTick','off','Box','off','TickLabelInterpreter','latex','FontSize',16,'YLim',[1e4 5e8], 'YTick',10.^(4:1:8));
yyaxis right;
area(T,ResidentTrajectory.Params.theta(T),'LineStyle','none','FaceColor',[.8 .8 .8],'FaceAlpha',.5);
set(gca,'YTick',[]);
yyaxis left;
xlabel('Time (hr)','FontSize',22, 'Interpreter','latex');
ylabel('Density (mL$^{-1}$)','FontSize',22,'Interpreter','latex');
legend('S','E$_1$','E$_2$','L','V','Location','best','Box','off','Interpreter','latex','FontSize',18);
xlim(ResidentTrajectory.Time{:}(end)+ [-3*ResidentTrajectory.Params.T 0]);

% Sample invasion trajectory -- no pleiotropy
nexttile(2);
T = [[ResidentTrajectory.Time{:}]; ResidentTrajectory.Time{:}(end)+[MutantTrajectories.Time{:}]]; % concatenate resident and mutant time series
Y = [[ResidentTrajectory.Densities{:}]; [MutantTrajectories.Densities{:}]]; % concatenate resident and mutant time series

resident = Y(:,2)+ Y(:,4) + 0.5*(Y(:,5)+Y(:,6)) + Y(:,8) + Y(:,10); % get resident density
mutant = Y(:,3)+ Y(:,7) + 0.5*(Y(:,5)+Y(:,6)) + Y(:,9) + Y(:,11); % get mutant density
mutantfraction = mutant./(mutant+resident); %get mutant fraction time series

[rate,prefactor] = MutantGrowthRate([MutantTrajectories.Time{:}],[MutantTrajectories.Densities{:}]); %get invasion growth rate of mutant
T_fit = (ResidentTrajectory.Time{:}(end)-ResidentTrajectory.Params.T):.1:(T(end)+ResidentTrajectory.Params.T); %define time sequence over which to plot fit
mutantfraction_fit = prefactor*exp(rate*(T_fit-ResidentTrajectory.Time{:}(end))); % calculate fit values

plot(T,mutantfraction,'-k','LineWidth',1.5); hold on; %plot mutant fraction
plot(T_fit,mutantfraction_fit,'-','LineWidth',2,'Color','r') % plot fit
set(gca,'TickLabelInterpreter','latex','FontSize',16,'YTick',0:.1:1,'Box','off');
xlabel('Time (hr)','FontSize',22,'Interpreter','latex');
ylabel('Mutant fraction','FontSize',22,'Interpreter','latex');
ylim([0 1]);
xlim([ResidentTrajectory.Time{:}(end)-4*ResidentTrajectory.Params.T T_fit(end)]);
xline(ResidentTrajectory.Time{:}(end),'LineWidth',1,'LineStyle','--','Color','k');
annotation_y = [.7 .7];
x_length = T_fit(end)+4*ResidentTrajectory.Params.T - ResidentTrajectory.Time{:}(end);
x_startloc = 12*ResidentTrajectory.Params.T;
x_endloc = 13*ResidentTrajectory.Params.T;
annotation_x = ([x_startloc x_endloc])/x_length;
annotation('textarrow','String','mutant added','Interpreter','latex',...
    'FontSize',16,'X', annotation_x,'Y',annotation_y);
text(.5*(T_fit(end)+ResidentTrajectory.Time{:}(end)),0.4,'$y = ae^{rt}$','Interpreter','latex','FontSize',16,'Color','r');
legend('Mutant fraction','Exponential Fit','Interpreter','Latex','Box','off','FontSize',18,'Location','best');

% Heatmaps of mutant growth rate -- no pleiotropy;
nexttile(3);
if ~isfile('Data/NoPleiotropy/Invasion_Period=3.00/Resident_z1=0.000,z2=1.000,MutantFrac=0.01.mat')
    [coESS, ResidentTrajectory, MutantTrajectories, GrowthRateTable] = PairInvasionDynamics([0 1], 3, ...
                                                                        MutantFraction = 0.01, Save = true, Pleiotropy = false);
end
load('Data/NoPleiotropy/Invasion_Period=3.00/Resident_z1=0.000,z2=1.000,MutantFrac=0.01.mat');

z1 = unique(GrowthRateTable.Mutant_z1);
z2 = unique(GrowthRateTable.Mutant_z2);

M = unstack(GrowthRateTable,"GrowthRate","Mutant_z2");
growthratematrix = table2array(M(:,2:end));

[Z1,Z2] = meshgrid(z1,z2);

imagesc(z1,z2,growthratematrix');
axis xy;
axis equal;
hold on;
contour(Z1,Z2,growthratematrix',[0 0]);
set(gca,'FontSize',16,'TickLabelInterpreter','latex','XTick',0:.2:1,'YTick',0:.2:1);
scatter(coESS(1),coESS(2),'filled','MarkerFaceColor','k');
xlabel('$\phi^1$','FontSize',22,'Interpreter','latex');
ylabel('$\phi^2$','FontSize',22,'Interpreter','latex');
cmap = colormap('hot');
c = colorbar;
c.Label.String = 'Growth rate (hr$^{-1}$)';
c.Label.Interpreter = 'latex';
c.Label.Rotation = -90;
c.TickLabelInterpreter = 'latex';
c.FontSize = 16;
c.Label.FontSize = 22;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% With Pleiotropy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coESS = [0.061 0.881];

% sample resident attractor -- yes pleiotropy
if ~isfile('Data/YesPleiotropy/Invasion_Period=3.00/Resident_z1=0.061,z2=0.881,MutantFrac=0.50.mat')
    [coESS, ResidentTrajectory, MutantTrajectories, GrowthRateTable] = PairInvasionDynamics(coESS, 3, MutantStrategies= [.5 .5],...
                                                                        MutantFraction = 0.5, Save = true, Pleiotropy = true);
end
load('Data/YesPleiotropy/Invasion_Period=3.00/Resident_z1=0.061,z2=0.881,MutantFrac=0.50.mat');

nexttile(4);
T = [ResidentTrajectory.Time{:}];
Y = [ResidentTrajectory.Densities{:}];

semilogy(T,Y(:,1),'LineWidth',1,'Color','b','DisplayName','S'); hold on;
semilogy(T,Y(:,2),'LineWidth',1,'Color','r','LineStyle','-','DisplayName','E$_1$');
semilogy(T,Y(:,4),'LineWidth',1,'Color','r','LineStyle','--','DisplayName','E$_2$');
semilogy(T,Y(:,8),'LineWidth',1,'Color',[165 42 42]/255,'LineStyle','-','DisplayName','L');
semilogy(T,Y(:,10),'LineWidth',1,'Color','g','LineStyle','-','DisplayName','V');
set(gca,'YMinorTick','off','Box','off','TickLabelInterpreter','latex','FontSize',16,'YLim',[1e4 5e8], 'YTick',10.^(4:1:8));
yyaxis right;
area(T,ResidentTrajectory.Params.theta(T),'LineStyle','none','FaceColor',[.8 .8 .8],'FaceAlpha',.5);
set(gca,'YTick',[]);
yyaxis left;
xlabel('Time (hr)','FontSize',22, 'Interpreter','latex');
ylabel('Density (mL$^{-1}$)','FontSize',22,'Interpreter','latex');
legend('S','E$_1$','E$_2$','L','V','Location','best','Box','off','Interpreter','latex','FontSize',18);
xlim(ResidentTrajectory.Time{:}(end)+ [-3*ResidentTrajectory.Params.T 0]);

% Sample invasion trajectory -- yes pleiotropy
nexttile(5);
T = [[ResidentTrajectory.Time{:}]; ResidentTrajectory.Time{:}(end)+[MutantTrajectories.Time{:}]]; % concatenate resident and mutant time series
Y = [[ResidentTrajectory.Densities{:}]; [MutantTrajectories.Densities{:}]]; % concatenate resident and mutant time series

resident = Y(:,2)+ Y(:,4) + 0.5*(Y(:,5)+Y(:,6)) + Y(:,8) + Y(:,10); % get resident density
mutant = Y(:,3)+ Y(:,7) + 0.5*(Y(:,5)+Y(:,6)) + Y(:,9) + Y(:,11); % get mutant density
mutantfraction = mutant./(mutant+resident); %get mutant fraction time series

[rate,prefactor] = MutantGrowthRate([MutantTrajectories.Time{:}],[MutantTrajectories.Densities{:}]); %get invasion growth rate of mutant
T_fit = (ResidentTrajectory.Time{:}(end)-ResidentTrajectory.Params.T):.1:(T(end)+ResidentTrajectory.Params.T); %define time sequence over which to plot fit
mutantfraction_fit = prefactor*exp(rate*(T_fit-ResidentTrajectory.Time{:}(end))); % calculate fit values

plot(T,mutantfraction,'-k','LineWidth',1.5); hold on; %plot mutant fraction
plot(T_fit,mutantfraction_fit,'-','LineWidth',2,'Color','r') % plot fit
set(gca,'TickLabelInterpreter','latex','FontSize',16,'YTick',0:.1:1,'Box','off');
xlabel('Time (hr)','FontSize',22,'Interpreter','latex');
ylabel('Mutant fraction','FontSize',22,'Interpreter','latex');
ylim([0 1]);
xlim([ResidentTrajectory.Time{:}(end)-4*ResidentTrajectory.Params.T T_fit(end)]);
xline(ResidentTrajectory.Time{:}(end),'LineWidth',1,'LineStyle','--','Color','k');
annotation_y = [.7 .7];
x_length = T_fit(end)+4*ResidentTrajectory.Params.T - ResidentTrajectory.Time{:}(end);
x_startloc = 4*ResidentTrajectory.Params.T;
x_endloc = 6*ResidentTrajectory.Params.T;
annotation_x = ([x_startloc x_endloc])/x_length;
annotation('textarrow','String','mutant added','Interpreter','latex',...
    'FontSize',16,'X', annotation_x,'Y',annotation_y);
text(.5*(T_fit(end)+ResidentTrajectory.Time{:}(end)),0.4,'$y = ae^{rt}$','Interpreter','latex','FontSize',16,'Color','r');
legend('Mutant fraction','Exponential Fit','Interpreter','Latex','Box','off','FontSize',18,'Location','best');

% Heatmaps of mutant growth rate -- yes pleiotropy;
nexttile(6);
if ~isfile('Data/YesPleiotropy/Invasion_Period=3.00/Resident_z1=0.061,z2=0.881,MutantFrac=0.01.mat')
    [coESS, ResidentTrajectory, MutantTrajectories, GrowthRateTable] = PairInvasionDynamics(coESS, 3, ...
                                                                        MutantFraction = 0.01, Save = true, Pleiotropy = true);
end
load('Data/YesPleiotropy/Invasion_Period=3.00/Resident_z1=0.061,z2=0.881,MutantFrac=0.01.mat');

z1 = unique(GrowthRateTable.Mutant_z1);
z2 = unique(GrowthRateTable.Mutant_z2);

M = unstack(GrowthRateTable,"GrowthRate","Mutant_z2");
growthratematrix = table2array(M(:,2:end));

[Z1,Z2] = meshgrid(z1,z2);

imagesc(z1,z2,growthratematrix');
axis xy;
axis equal;
hold on;
contour(Z1,Z2,growthratematrix',[0 0]);
set(gca,'FontSize',16,'TickLabelInterpreter','latex','XTick',0:.2:1,'YTick',0:.2:1);
scatter(coESS(1),coESS(2),'filled','MarkerFaceColor','k');
xlabel('$\phi^1$','FontSize',22,'Interpreter','latex');
ylabel('$\phi^2$','FontSize',22,'Interpreter','latex');
cmap = colormap('hot');
c = colorbar;
c.Label.String = 'Growth rate (hr$^{-1}$)';
c.Label.Interpreter = 'latex';
c.Label.Rotation = -90;
c.TickLabelInterpreter = 'latex';
c.FontSize = 16;
c.Label.FontSize = 22;


