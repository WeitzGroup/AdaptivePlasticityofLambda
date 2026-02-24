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
[coESS, OrderedStrategies, PopFractions, Time, TraitMean, TraitStdev] = MultiSpeciesTrajectories([0 1], 3, 100, 2000, Pleiotropy= false);

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
plot(Time(1:30:end), TraitMean(1:30:end,1),'LineWidth',2,'Color','r'); hold on;
plot(Time(1:30:end), TraitMean(1:30:end,2),'LineWidth',2,'Color','r','LineStyle','--');
yline(coESS(1),'--');
yline(coESS(2),'--');
set(gca,'FontSize',14,'TickLabelInterpreter','latex','YLim',[-.01 1.01],'box','off');
legend('$\bar{\phi^1}$','$\bar{\phi^2}$','Interpreter','latex','FontWeight','Bold','FontSize',16,'location','best','box','off');
xlabel('Time (hr)','Interpreter','latex','FontSize',22);
ylabel('Mean trait value','Interpreter','latex','FontSize',22);

%Standard deviation trait trajectory -- no pleiotropy
nexttile(3);
plot(Time(1:30:end),TraitStdev(1:30:end,1),'LineWidth',2,'Color','r'); hold on;
plot(Time(1:30:end),TraitStdev(1:30:end,2),'LineWidth',2,'Color','r','LineStyle','--');
set(gca,'FontSize',14,'TickLabelInterpreter','latex','box','off');
legend('$\sigma_{\phi^1}$','$\sigma_{\phi^2}$','Interpreter','latex','FontWeight','Bold','FontSize',16,'location','best','box','off');
xlabel('Time (hr)','Interpreter','latex','FontSize',22);
ylabel('Standard deviation of trait value','Interpreter','latex','FontSize',22);


%% Yes Pleiotropy case %%%%%%%%%%%%
[coESS, OrderedStrategies, PopFractions, Time, TraitMean, TraitStdev] = MultiSpeciesTrajectories([0.06 .88], 3, 100, 2000, Pleiotropy= true);

% Scatter plot of strategies -- yes pleiotropy
nexttile(1);
scatter(OrderedStrategies(:,1),OrderedStrategies(:,2),40,'r','Marker','x');
hold on;
scatter(coESS(1),coESS(2),80,0*[1 1 1],'filled','Marker','o');
pbaspect([1 1 1]);
set(gca,'FontSize',14,'TickLabelInterpreter','latex','XTick',0:.2:1,'YTick',0:.2:1,'Box','on');
xlabel('$\phi^1$','Interpreter','latex','FontSize',22);
ylabel('$\phi^2$','Interpreter','latex','FontSize',22);

%Mean trait trajectory -- yes pleiotropy
nexttile(2);
plot(Time(1:30:end), TraitMean(1:30:end,1),'LineWidth',2,'Color','r'); hold on;
plot(Time(1:30:end), TraitMean(1:30:end,2),'LineWidth',2,'Color','r','LineStyle','--');
yline(coESS(1),'--');
yline(coESS(2),'--');
set(gca,'FontSize',14,'TickLabelInterpreter','latex','YLim',[-.01 1.01],'box','off');
legend('$\bar{\phi^1}$','$\bar{\phi^2}$','Interpreter','latex','FontWeight','Bold','FontSize',16,'location','best','box','off');
xlabel('Time (hr)','Interpreter','latex','FontSize',22);
ylabel('Mean trait value','Interpreter','latex','FontSize',22);

%Standard deviation trait trajectory -- yes pleiotropy
nexttile(3);
plot(Time(1:30:end),TraitStdev(1:30:end,1),'LineWidth',2,'Color','r'); hold on;
plot(Time(1:30:end),TraitStdev(1:30:end,2),'LineWidth',2,'Color','r','LineStyle','--');
set(gca,'FontSize',14,'TickLabelInterpreter','latex','box','off');
legend('$\sigma_{\phi^1}$','$\sigma_{\phi^2}$','Interpreter','latex','FontWeight','Bold','FontSize',16,'location','best','box','off');
xlabel('Time (hr)','Interpreter','latex','FontSize',22);
ylabel('Standard deviation of trait value','Interpreter','latex','FontSize',22);
