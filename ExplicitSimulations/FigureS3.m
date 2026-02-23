% This script generates figure S3 of the manuscript which creates the
% 1. Sampled strategies for multispecies dynamics
% 2. Mean trait values over time
% 3. Variance of trait values over time
% for both, pleiotropic and non-pleiotropic cases. 

% The data for these plots are contained in
% 'Data/MultispeciesDynamics_NoPleiotropy_n=100,Period=3,z1=0.00,z2=1.00.mat'
% and 'Data/MultispeciesDynamics_YesPleiotropy_n=100,Period=3,z1=0.06,z2=0.88.mat'

%Date: July 16, 2025
%Author: Tapan Goel


clear all;
close all;
%% Generate figure
addpath('Utils\');
f = figure('Position',[200,100,1904,935]);
t = tiledlayout(2,3);


%% No Pleiotropy case %%%%%%%%%%%%

load('Data/MultispeciesDynamics_NoPleiotropy_n=100,Period=3,z1=0.00,z2=1.00.mat');

% Scatter plot of strategies -- no pleiotropy
nexttile(1);
scatter(params.z(:,1),params.z(:,2),40,'r','Marker','x');
hold on;
scatter(params.z(1,1),params.z(1,2),80,0*[1 1 1],'filled','Marker','o');
pbaspect([1 1 1]);
set(gca,'FontSize',14,'TickLabelInterpreter','latex','XTick',0:.2:1,'YTick',0:.2:1,'Box','on');
xlabel('$\phi^1$','Interpreter','latex','FontSize',22);
ylabel('$\phi^2$','Interpreter','latex','FontSize',22);

%Mean trait trajectory -- no pleiotropy
nexttile(2);
plot(T(1:30:end),mean_phi(1:30:end,1),'LineWidth',2,'Color','r'); hold on;
plot(T(1:30:end), mean_phi(1:30:end,2),'LineWidth',2,'Color','r','LineStyle','--');
yline(strategies_ordered(end,1),'--');
yline(strategies_ordered(end,2),'--');
set(gca,'FontSize',14,'TickLabelInterpreter','latex','YLim',[-.01 1.01],'box','off');
legend('$\bar{\phi^1}$','$\bar{\phi^2}$','Interpreter','latex','FontWeight','Bold','FontSize',16,'location','best','box','off');
xlabel('Time (hr)','Interpreter','latex','FontSize',22);
ylabel('Mean trait value','Interpreter','latex','FontSize',22);

%Standard deviation trait trajectory -- no pleiotropy
nexttile(3);
std_phi = sqrt(popfractions*(strategies_ordered.^2)-mean_phi.^2);
plot(T(1:30:end),std_phi(1:30:end,1),'LineWidth',2,'Color','r'); hold on;
plot(T(1:30:end), std_phi(1:30:end,2),'LineWidth',2,'Color','r','LineStyle','--');
%yline(strategies_ordered(1,1),'--');
%yline(strategies_ordered(1,2),'--');
set(gca,'FontSize',14,'TickLabelInterpreter','latex','box','off');
legend('$\sigma_{\phi^1}$','$\sigma_{\phi^2}$','Interpreter','latex','FontWeight','Bold','FontSize',16,'location','best','box','off');
xlabel('Time (hr)','Interpreter','latex','FontSize',22);
ylabel('Standard deviation of trait value','Interpreter','latex','FontSize',22);


% clear data that is not in use anymore
clear E1 iter1 L mean_phi params popfractions S std_phi strategies_ordered T V;

%% Yes Pleiotropy case %%%%%%%%%%%%
load('Data/MultispeciesDynamics_YesPleiotropy_n=100,Period=3,z1=0.06,z2=0.88.mat');

% Scatter plot of strategies -- yes pleiotropy
nexttile(4);
scatter(params.z(:,1),params.z(:,2),40,'k','Marker','x');
hold on;
scatter(params.z(1,1),params.z(1,2),80,0*[.5 .5 .5],'filled','Marker','o');
pbaspect([1 1 1]);
set(gca,'FontSize',14,'TickLabelInterpreter','latex','XTick',0:.2:1,'YTick',0:.2:1,'box','on');
xlabel('$\phi^1$','Interpreter','latex','FontSize',22);
ylabel('$\phi^2$','Interpreter','latex','FontSize',22);

%Mean trait trajectory -- yes pleiotropy
nexttile(5);
plot(T(1:30:end),mean_phi(1:30:end,1),'LineWidth',2,'Color','r'); hold on;
plot(T(1:30:end), mean_phi(1:30:end,2),'LineWidth',2,'Color','r','LineStyle','--');
yline(strategies_ordered(end,1),'--');
yline(strategies_ordered(end,2),'--');
set(gca,'FontSize',14,'TickLabelInterpreter','latex','YLim',[-.01 1.01],'box','off');
legend('$\bar{\phi^1}$','$\bar{\phi^2}$','Interpreter','latex','FontWeight','Bold','FontSize',16,'location','best','box','off');
xlabel('Time (hr)','Interpreter','latex','FontSize',22);
ylabel('Mean trait value','Interpreter','latex','FontSize',22);

%Standard deviation trait trajectory -- yes pleiotropy
nexttile(6);
std_phi = sqrt(popfractions*(strategies_ordered.^2)-mean_phi.^2);
plot(T(1:30:end),std_phi(1:30:end,1),'LineWidth',2,'Color','r'); hold on;
plot(T(1:30:end), std_phi(1:30:end,2),'LineWidth',2,'Color','r','LineStyle','--');
%yline(strategies_ordered(1,1),'--');
%yline(strategies_ordered(1,2),'--');
set(gca,'FontSize',14,'TickLabelInterpreter','latex','Box','off');
legend('$\sigma_{\phi^1}$','$\sigma_{\phi^2}$','Interpreter','latex','FontWeight','Bold','FontSize',16,'location','best','box','off');
xlabel('Time (hr)','Interpreter','latex','FontSize',22);
ylabel('Standard deviation of trait value','Interpreter','latex','FontSize',22);
