%% plot Conv(Gamma) and D under different number of samples
% dependencies: Yalmip
% 2024-02-21
%  By Jiayi Ding & Shuai Lu

%% Initialize
clc;clear;close all;
warning off;
yalmip('clear');

%% Import data
filepath = [cd '\data\data.mat'];
temp = load(filepath);
data = temp.data;

%% new figure
h_fig = figure();          % gcf: get current figure
h_axis = gca;              % gca: get current axis

%% set position & color
% position, color,
left = 10; bottom = 0; width = 20; height = 18;
% units:inches|centimeters|normalized|points|{pixels}|characters
set(h_fig, 'Units','centimeters', 'position', [left, bottom, width, height], 'color', 'w');

%% Remove the blank edge
set(gca,'LooseInset',get(gca,'TightInset'));

%% num = 20
num_sample = 20;
data.P = data.dataset.P_agg(1:num_sample,:);
data.price = data.dataset.price(1:num_sample,:);
[k,av] = convhull(data.P);

% % plot conv(Gamma) 
legend_1_1 = plot(data.P(k, 1), data.P(k,2), LineWidth=1, Color=[0 0 0]); hold on;
legend_1_2 = fill(data.P(k, 1), data.P(k,2), [0 0 0]);
set(legend_1_2,'edgealpha',0,'facealpha',0.2);
hold on;

% % plot D
D = sdpvar(2, 1);
Constraints_D = [];
Constraints_D = Constraints_D + ( ...
    data.price(k, 1:2)*D >= ...
    sum(data.price(k, 1:2).*data.dataset.P_agg(k, 1:2),2));
v = vertex(Constraints_D,D);
[k,av] = convhull(v');
legend_1_3 = plot(v(1,k),v(2,k), 'LineWidth',5,'LineStyle',':','Color',[0 0 0]);
hold on;

%% num = 50
num_sample=50;
data.P = data.dataset.P_agg(1:num_sample,:);
data.price = data.dataset.price(1:num_sample,:);
[k,av] = convhull(data.P);

% % plot conv(Gamma) 
legend_2_1=plot(data.P(k, 1), data.P(k,2), LineWidth=1, Color=[0 0 1]); hold on;
legend_2_2 = fill(data.P(k, 1), data.P(k,2), [0 0.4 1]);
hold on;
set(legend_2_2,'edgealpha',0,'facealpha',0.2);

% % plot D
D = sdpvar(2, 1);
Constraints_D = [];
Constraints_D = Constraints_D + (...
    data.price(k, 1:2)*D >= ...
    sum(data.price(k, 1:2).*data.dataset.P_agg(k, 1:2),2));
v = vertex(Constraints_D,D);
[k,av] = convhull(v');
legend_2_3 = plot(v(1,k),v(2,k), 'LineWidth',3,'LineStyle','-.','Color',[0 0 1]);
hold on;

%% num = 200
num_sample=200;
data.P = data.dataset.P_agg(1:num_sample,:);
data.price = data.dataset.price(1:num_sample,:);
[k,av] = convhull(data.P);

% % plot conv(Gamma) 
legend_3_1 = plot(data.P(k, 1), data.P(k,2), LineWidth=1, Color=[1 0 0]); hold on;
legend_3_2 = fill(data.P(k, 1), data.P(k,2), [1 0 0]);
hold on;
set(legend_3_2,'edgealpha',0,'facealpha',0.2);

% % plot D
D = sdpvar(2, 1);
Constraints_D = [];
Constraints_D = Constraints_D + (...
    data.price(k, 1:2)*D >= ...
    sum(data.price(k, 1:2).*data.dataset.P_agg(k, 1:2),2));
v = vertex(Constraints_D,D);
[k,av] = convhull(v');
legend_3_3 = plot(v(1,k),v(2,k), 'LineWidth',1,'LineStyle','--','Color',[1 0 0]);
hold on;

%% set asix
axis([-100 700 -200 900]);
set(h_axis, 'XTick', 0:200:600);                        % or xticks();
set(h_axis, 'YTick', -200:200:800);                     % or xticks();

%% lengend
legend1 = legend([legend_1_2,legend_1_3,legend_2_2,legend_2_3,legend_3_2,legend_3_3], ...
    '$\ \rm{Conv}(\Gamma) (|\rm{K}|=20)$','$\ \partial \Pi (|\rm{K}|=20)$', ...
    '$\ \rm{Conv}(\Gamma) (|\rm{K}|=50)$','$\ \partial \Pi (|\rm{K}|=50)$', ...
    '$\ \rm{Conv}(\Gamma) (|\rm{K}|=200)$','$\ \partial \Pi (|\rm{K}|=200)$', ...
    'Orientation','horizontal','NumColumns',2, ...
    'Interpreter','latex');

% Create xlabel
xlabel('$P^1$ (kW)', 'Interpreter','latex');
% Create ylabel
ylabel('$P^2$ (kW)', 'Interpreter','latex');

hold on; grid on;

%% Font
set(h_axis, 'FontName', 'Times New Roman', 'FontSize', 28);
set(legend1, 'FontSize', 24, 'Location', 'north', 'Color','none', 'EdgeColor','none');

%% style
set(h_axis,'Color','none');

%% save
filepath = [cd '\results\'];
str = [filepath 'Fig_1'];
saveas(gca, str, 'fig');
saveas(gca, str, 'svg');
saveas(gca, str, 'emf');
