% This script generates figure S3 of the manuscript which creates the
% 1. Sampled strategies for multispecies dynamics
% 2. Mean trait values over time
% 3. Variance of trait values over time
% for both, pleiotropic and non-pleiotropic cases. 


%Date: 02/23/2026
%Author: Tapan Goel


clear all;
close all;
%% Generate figure
addpath('Utils/');
f = figure('Position',[200,100,1904,935]);
t = tiledlayout(2,3);


%% No Pleiotropy case %%%%%%%%%%%%
coESS = [0 1];
if ~isfile('Data/NoPleiotropy/MultiSpeciesDynamics_Period=3.00,n=100,z1=0.000,z2=1.000.mat')
    [coESS, OrderedStrategies, PopFractions, T, TraitMean, TraitStdev] = MultiSpeciesTrajectories(coESS, 3, 100, 2000, Pleiotropy= false, Save= true);
end
load('Data/NoPleiotropy/MultiSpeciesDynamics_Period=3.00,n=100,z1=0.000,z2=1.000.mat');
% Scatter plot of strategies -- no pleiotropy
nexttile(1);
scatter(OrderedStrategies(:,1),OrderedStrategies(:,2),40,'r','Marker','x');
hold on;
scatter(coESS(1),coESS(2),80,0*[1 1 1],'filled','Marker','o');
pbaspect([1 1 1]);
set(gca,'FontSize',14,'TickLabelInterpreter','latex','XTick',0:.2:1,'YTick',0:.2:1,'Box','on');
xlabel('$\phi^1$','Interpreter','latex','FontSize',22);
ylabel('$\phi^2$','Interpreter','latex','FontSize',22);

%Mean trait trajectory -- no pleiotropy
nexttile(2);
plot(T(1:30:end), TraitMean(1:30:end,1),'LineWidth',2,'Color','r'); hold on;
plot(T(1:30:end), TraitMean(1:30:end,2),'LineWidth',2,'Color','r','LineStyle','--');
yline(coESS(1),'--');
yline(coESS(2),'--');
set(gca,'FontSize',14,'TickLabelInterpreter','latex','YLim',[-.01 1.01],'box','off');
legend('$\bar{\phi^1}$','$\bar{\phi^2}$','Interpreter','latex','FontWeight','Bold','FontSize',16,'location','best','box','off');
xlabel('Time (hr)','Interpreter','latex','FontSize',22);
ylabel('Mean trait value','Interpreter','latex','FontSize',22);

%Standard deviation trait trajectory -- no pleiotropy
nexttile(3);
plot(T(1:30:end),TraitStdev(1:30:end,1),'LineWidth',2,'Color','r'); hold on;
plot(T(1:30:end),TraitStdev(1:30:end,2),'LineWidth',2,'Color','r','LineStyle','--');
set(gca,'FontSize',14,'TickLabelInterpreter','latex','box','off');
legend('$\sigma_{\phi^1}$','$\sigma_{\phi^2}$','Interpreter','latex','FontWeight','Bold','FontSize',16,'location','best','box','off');
xlabel('Time (hr)','Interpreter','latex','FontSize',22);
ylabel('Standard deviation of trait value','Interpreter','latex','FontSize',22);

clear coESS OrderedStrategies PopFractions T TraitStdev TraitMean;

%% Yes Pleiotropy case %%%%%%%%%%%%
coESS = [0.061 0.881];
if ~isfile('Data/YesPleiotropy/MultiSpeciesDynamics_Period=3.00,n=100,z1=0.061,z2=0.881.mat')
    [coESS, OrderedStrategies, PopFractions, T, TraitMean, TraitStdev] = MultiSpeciesTrajectories(coESS, 3, 100, 2000, Pleiotropy= false, Save= true);
end
load('Data/YesPleiotropy/MultiSpeciesDynamics_Period=3.00,n=100,z1=0.061,z2=0.881.mat');

%Scatter plot of strategies -- yes pleiotropy
nexttile(4);
scatter(OrderedStrategies(:,1),OrderedStrategies(:,2),40,'r','Marker','x');
hold on;
scatter(coESS(1),coESS(2),80,0*[1 1 1],'filled','Marker','o');
pbaspect([1 1 1]);
set(gca,'FontSize',14,'TickLabelInterpreter','latex','XTick',0:.2:1,'YTick',0:.2:1,'Box','on');
xlabel('$\phi^1$','Interpreter','latex','FontSize',22);
ylabel('$\phi^2$','Interpreter','latex','FontSize',22);

%Mean trait trajectory -- yes pleiotropy
nexttile(5);
plot(T(1:30:end), TraitMean(1:30:end,1),'LineWidth',2,'Color','r'); hold on;
plot(T(1:30:end), TraitMean(1:30:end,2),'LineWidth',2,'Color','r','LineStyle','--');
yline(coESS(1),'--');
yline(coESS(2),'--');
set(gca,'FontSize',14,'TickLabelInterpreter','latex','YLim',[-.01 1.01],'box','off');
legend('$\bar{\phi^1}$','$\bar{\phi^2}$','Interpreter','latex','FontWeight','Bold','FontSize',16,'location','best','box','off');
xlabel('Time (hr)','Interpreter','latex','FontSize',22);
ylabel('Mean trait value','Interpreter','latex','FontSize',22);

%Standard deviation trait trajectory -- yes pleiotropy
nexttile(6);
plot(T(1:30:end),TraitStdev(1:30:end,1),'LineWidth',2,'Color','r'); hold on;
plot(T(1:30:end),TraitStdev(1:30:end,2),'LineWidth',2,'Color','r','LineStyle','--');
set(gca,'FontSize',14,'TickLabelInterpreter','latex','box','off');
legend('$\sigma_{\phi^1}$','$\sigma_{\phi^2}$','Interpreter','latex','FontWeight','Bold','FontSize',16,'location','best','box','off');
xlabel('Time (hr)','Interpreter','latex','FontSize',22);
ylabel('Standard deviation of trait value','Interpreter','latex','FontSize',22);
