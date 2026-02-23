% Ths script generates figure S2 of the manuscript which creates the
% 1. sample resident attractor dynamics
% 2. sample invasion dynamics
% 3. Heatmap of invasion fitness
% for both, pleiotropic and non-pleiotropic cases. 
% The data for generating these plots is stored in 'Data/FigureData.mat'

%Date: July 16, 2025
%Author: Tapan Goel

clear all;
close all;
%% Generate figure
addpath('Utils\');
load('Data/FigureData.mat'); %load file that contains data for the figure.
f = figure('Position',[160,190,1904,935]);
t = tiledlayout(2,3);

%Sample resident attractor -- no pleiotropy
nexttile(1)
T = NoPleiotropyInvasionTrajectory.T;
Y = NoPleiotropyInvasionTrajectory.Y;

semilogy(T,Y(:,1),'LineWidth',1,'Color','b','DisplayName','S'); hold on;
semilogy(T,Y(:,2),'LineWidth',1,'Color','r','LineStyle','-','DisplayName','E$_1$');
semilogy(T,Y(:,4),'LineWidth',1,'Color','r','LineStyle','--','DisplayName','E$_2$');
semilogy(T,Y(:,8),'LineWidth',1,'Color',[165 42 42]/255,'LineStyle','-','DisplayName','L');
semilogy(T,Y(:,10),'LineWidth',1,'Color','g','LineStyle','-','DisplayName','V');
set(gca,'YMinorTick','off','Box','off','TickLabelInterpreter','latex','FontSize',16,'YLim',[1e4 5e8], 'YTick',10.^(4:1:8));
yyaxis right;
area(T,NoPleiotropyInvasionTrajectory.params.theta(T),'LineStyle','none','FaceColor',[.8 .8 .8],'FaceAlpha',.5);
set(gca,'YTick',[]);
yyaxis left;
xlabel('Time (hr)','FontSize',22, 'Interpreter','latex');
ylabel('Density (mL$^{-1}$)','FontSize',22,'Interpreter','latex');
legend('S','E$_1$','E$_2$','L','V','Location','best','Box','off','Interpreter','latex','FontSize',18);
xlim(T(end)+ [-3*NoPleiotropyInvasionTrajectory.params.T 0]);

% Sample invasion trajectory -- no pleiotropy
nexttile(2);
T1 = [NoPleiotropyInvasionTrajectory.T; NoPleiotropyInvasionTrajectory.T1+NoPleiotropyInvasionTrajectory.T(end)];
Y1 = [NoPleiotropyInvasionTrajectory.Y; NoPleiotropyInvasionTrajectory.Y1];
T = T1(T1>NoPleiotropyInvasionTrajectory.T(end)-4*NoPleiotropyInvasionTrajectory.params.T);
Y = Y1(T1>NoPleiotropyInvasionTrajectory.T(end)-4*NoPleiotropyInvasionTrajectory.params.T,:);

resident = Y(:,2)+ Y(:,4) + 0.5*(Y(:,5)+Y(:,6)) + Y(:,8) + Y(:,10);
mutant = Y(:,3)+ Y(:,7) + 0.5*(Y(:,5)+Y(:,6)) + Y(:,9) + Y(:,11);
mutantfraction = mutant./(mutant+resident);

% Obtain mutant growth rate -- no pleiotropy
[rate,prefactor] = MutantGrowthRate(NoPleiotropyInvasionTrajectory.T1,NoPleiotropyInvasionTrajectory.Y1);

%plot mutant fraction
plot(T,mutantfraction,'-k','LineWidth',1.5); hold on;
temp_time = [(-3:.1:-.1)';NoPleiotropyInvasionTrajectory.T1;NoPleiotropyInvasionTrajectory.T1(end)+(.1:.1:3)']+NoPleiotropyInvasionTrajectory.T(end);
plot(temp_time,prefactor*exp((temp_time-NoPleiotropyInvasionTrajectory.T(end))*rate),'-','LineWidth',2,'Color',0.5*[1 1 1]);
set(gca,'TickLabelInterpreter','latex','FontSize',16,'YTick',0:.1:1,'Box','off');
xlabel('Time (hr)','FontSize',22,'Interpreter','latex');
ylabel('Mutant fraction','FontSize',22,'Interpreter','latex');
ylim([0 1]);
xline(NoPleiotropyInvasionTrajectory.T(end),'LineWidth',1,'LineStyle','--','Color','k');
annotation('textarrow','String','mutant added','Interpreter','latex',...
    'FontSize',16,'X',[NoPleiotropyInvasionTrajectory.T(end)/T(end) NoPleiotropyInvasionTrajectory.T(end)/T(end)+.1],...
    'Y',[.7 .7]);
legend('Mutant fraction','Exponential Fit','Interpreter','Latex','Box','off','FontSize',18,'Location','best');

% Heatmaps of mutant growth rate -- no pleiotropy;
nexttile(3);
z1 = unique(NoPleiotropyDataTable.z1);
z2 = unique(NoPleiotropyDataTable.z2);

M = unstack(NoPleiotropyDataTable,'growth_rate','z2');
growthratematrix = table2array(M(:,2:end));

[Z1,Z2] = meshgrid(z1,z2);

imagesc(z1,z2,growthratematrix');
axis xy;
axis equal;
hold on;
contour(Z1,Z2,growthratematrix',[0 0],'-k');
set(gca,'FontSize',16,'TickLabelInterpreter','latex','XTick',0:.2:1,'YTick',0:.2:1);
scatter(NoPleiotropyESS(1),NoPleiotropyESS(2),'filled','MarkerFaceColor','k');
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

%sample resident trajectory -- yes pleiotropy
nexttile(4)
T = YesPleiotropyInvasionTrajectory.T;
Y = YesPleiotropyInvasionTrajectory.Y;

semilogy(T,Y(:,1),'LineWidth',1,'Color','b','DisplayName','S'); hold on;
semilogy(T,Y(:,2),'LineWidth',1,'Color','r','LineStyle','-','DisplayName','E$_1$');
semilogy(T,Y(:,4),'LineWidth',1,'Color','r','LineStyle','--','DisplayName','E$_2$');
semilogy(T,Y(:,8),'LineWidth',1,'Color',[225 42 42]/255,'LineStyle','-','DisplayName','L');
semilogy(T,Y(:,10),'LineWidth',1,'Color','g','LineStyle','-','DisplayName','V');
set(gca,'YMinorTick','off','Box','off','TickLabelInterpreter','latex','FontSize',16,'YLim',[1e4 5e8], 'YTick',10.^(4:1:8));
yyaxis right;
area(T,YesPleiotropyInvasionTrajectory.params.theta(T),'LineStyle','none','FaceColor',[.8 .8 .8],'FaceAlpha',.5);
set(gca,'YTick',[]);
yyaxis left;
xlabel('Time (hr)','FontSize',22, 'Interpreter','latex');
ylabel('Density (mL$^{-1}$)','FontSize',22,'Interpreter','latex');
legend('S','E$_1$','E$_2$','L','V','Location','best','Box','off','Interpreter','latex','FontSize',18);
xlim(T(end)+ [-3*YesPleiotropyInvasionTrajectory.params.T 0]);

% Sample invasion trajectory -- Yes pleiotropy
nexttile(5);
T1 = [YesPleiotropyInvasionTrajectory.T; YesPleiotropyInvasionTrajectory.T1+YesPleiotropyInvasionTrajectory.T(end)];
Y1 = [YesPleiotropyInvasionTrajectory.Y; YesPleiotropyInvasionTrajectory.Y1];
T = T1(T1>YesPleiotropyInvasionTrajectory.T(end)-4*YesPleiotropyInvasionTrajectory.params.T);
Y = Y1(T1>YesPleiotropyInvasionTrajectory.T(end)-4*YesPleiotropyInvasionTrajectory.params.T,:);

resident = Y(:,2)+ Y(:,4) + 0.5*(Y(:,5)+Y(:,6)) + Y(:,8) + Y(:,10);
mutant = Y(:,3)+ Y(:,7) + 0.5*(Y(:,5)+Y(:,6)) + Y(:,9) + Y(:,11);
mutantfraction = mutant./(mutant+resident);

% Obtain mutant growth rate -- no pleiotropy
[rate,prefactor] = MutantGrowthRate(YesPleiotropyInvasionTrajectory.T1,YesPleiotropyInvasionTrajectory.Y1);
%plot mutant fraction
plot(T,mutantfraction,'-k','LineWidth',1.5); hold on;
temp_time = [(-3:.1:-.1)';YesPleiotropyInvasionTrajectory.T1;YesPleiotropyInvasionTrajectory.T1(end)+(.1:.1:3)']+YesPleiotropyInvasionTrajectory.T(end);
plot(temp_time,prefactor*exp((temp_time-YesPleiotropyInvasionTrajectory.T(end))*rate),'-','LineWidth',2,'Color',.5*[1 1 1]);
set(gca,'TickLabelInterpreter','latex','FontSize',16,'YTick',0:.1:1,'Box','off');
xlabel('Time (hr)','FontSize',22,'Interpreter','latex');
ylabel('Mutant fraction','FontSize',22,'Interpreter','latex');
ylim([0 1]);
xline(YesPleiotropyInvasionTrajectory.T(end),'LineWidth',1,'LineStyle','--','Color','k');
%annotation('textarrow','String','mutant added','Interpreter','latex',...
%    'FontSize',16,'X',[YesPleiotropyInvasionTrajectory.T(end)/T(end) YesPleiotropyInvasionTrajectory.T(end)/T(end)+.1],...
%    'Y',[.7 .7]);
legend('Mutant fraction','Exponential Fit','Interpreter','Latex','Box','off','FontSize',18,'Location','best');

% Heatmaps of mutant growth rate -- Yes pleiotropy;
nexttile(6);
z1 = unique(YesPleiotropyDataTable.z1);
z2 = unique(YesPleiotropyDataTable.z2);
M = unstack(YesPleiotropyDataTable,'growth_rate','z2');
growthratematrix = table2array(M(:,2:end));
[Z1,Z2] = meshgrid(z1,z2);
imagesc(z1,z2,growthratematrix');
axis xy;
axis equal
hold on;
contour(Z1,Z2,growthratematrix',[0 0]);
set(gca,'FontSize',16,'TickLabelInterpreter','latex','XTick',0:.1:1,'YTick',0:.1:1);
scatter(YesPleiotropyESS(1),YesPleiotropyESS(2),'filled','MarkerFaceColor','k');
xlabel('$\phi^1$','FontSize',22,'Interpreter','latex');
ylabel('$\phi^2$','FontSize',22,'Interpreter','latex');
colormap('hot');
c = colorbar;
c.Label.String = 'Growth rate (hr$^{-1}$)';
c.Label.Interpreter = 'latex';
c.Label.Rotation = -90;
c.TickLabelInterpreter = 'latex';
c.FontSize = 16;
c.Label.FontSize = 22;


