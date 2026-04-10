%% This script generates figure 1 of the manuscript.
% Data for the figures courtesy of Ido Golding - data from Zeng et al. 2010 (DOI: 10.1016/j.cell.2010.03.034)

% Author: Joshua Weitz
% Date: 07/15/2025
% Modified: 04/10/2026, Tapan Goel, combined scripts to plot figures 1a and
% 1b


%% Data courtesy of Ido Golding - data from Zeng et al. 2010 (DOI: 10.1016/j.cell.2010.03.034)
MOI = [1 2 3 4 5];  % row % Multiplicity of infection
length = [0.7 0.9 1.1 1.3 1.5]; % col % Normalized cell length

plys = 1- [
    0.5532    0.7298    0.8304    0.8441    0.8967
    0.5241    0.7150    0.7246    0.7367    0.7595
    0.3701    0.5206    0.6190    0.7604    0.7730
    0.3260    0.4464    0.5988    0.5933    0.6814
    0.2143    0.2497    0.3953    0.4992    0.6631
]; % Probability of lysogeny. Entry (i,j) is the probability of lysogeny at MOI i with normalized cell length length(j).

%% Generate Figure 1a
colors = [0 0 0;...
          170 68 0;...
          212 85 0;...
          255 153 85;...
          255 204 170;...
          255 255 255];

colors = colors/255; %% Colors for markers for different MOIs

clf;
% automatically create postscript whenever
% figure is drawn

tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

set(gcf,'DefaultLineMarkerSize',10);
set(gcf,'DefaultAxesLineWidth',2);
set(gcf,'PaperPositionMode','auto');

tmpm = ['o','o','o','o','o'];
for i=1:5
  tmph=errorbar(MOI(i),mean(plys(i,:)),std(plys(i,:))/sqrt(size(plys(i,:),2))); %% Divide the std by square root of number of samples.
  set(tmph,'marker',tmpm(i));
  set(tmph,'markersize',10);
  set(tmph,'color',colors(i,:));
  set(tmph,'markerfacecolor',colors(i,:));
  hold on
end

set(gca,'fontsize',20);
ylim([0 1]);
xlim([0.5 5.5]);
% fix up tickmarks
set(gca,'xtick',[1:1:5]);
set(gca,'ytick',[0:0.2:1]);

xlabel('Viral MOI, $k$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Probability of lysogeny, $\phi^k$','fontsize',20,'verticalalignment','bottom','interpreter','latex');

exportgraphics(gcf,'Figure1a.pdf');

clear tmp*

%% Generate Figure 1b
clf;

tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

set(gcf,'DefaultLineMarkerSize',10);
set(gcf,'DefaultAxesLineWidth',2);
set(gcf,'PaperPositionMode','auto');

tmpm = ['o','o','o','o','o'];
for i=1:5
  tmph=plot(MOI(i)./length,plys(i,:),tmpm(i));
  set(tmph,'markersize',10);
  set(tmph,'color',colors(i,:));
  set(tmph,'markerfacecolor',colors(i,:));
  hold on
end

set(gca,'fontsize',20);
ylim([0 1]);
set(gca,'xtick',[0:1:8]);
% legend
for i=1:5,
  tmps{i}=sprintf('$k = %d$',MOI(i));
end
tmplh = legend(tmps,'Location','NorthWest');
set(tmplh,'interpreter','latex');
legend('boxoff');

xlabel('Viral concentration, $k/\mathrm{vol}$ (1/$\mu$m)','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Probability of lysogeny, $\phi^{k}$','fontsize',20,'verticalalignment','bottom','interpreter','latex');

exportgraphics(gcf,'Figure1b.pdf');

clear tmp*
