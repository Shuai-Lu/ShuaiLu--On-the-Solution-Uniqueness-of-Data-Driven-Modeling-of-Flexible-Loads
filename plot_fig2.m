%% plot Ωphy under different numbers of virtual battery
%  2024-02-21
%  By Jiayi Ding & Shuai Lu

%% Initialize
clc;clear;close all;
warning off;
yalmip('clear');

% %% Declaring Global Variables
% global model;

model.set1 = [];
model.set1.var = [];
model.set1.cons = [];

model.set2 = [];
model.set2.var = [];
model.set2.cons = [];

%% Import data
filepath = [cd '\data\identification_1_50.mat'];
temp = load(filepath);
model.settings = temp.model.settings;
model.set1.para = temp.model.results;
model.set1.data = temp.model.para;

filepath = [cd '\data\identification_2_50.mat'];
temp = load(filepath);
model.set2.para = temp.model.results;
model.set2.data = temp.model.para;

%%
h_fig = figure();          % gcf: get current figure
h_axis = gca;              % gca: get current axis

%% set position & color
% position, color,
left = 1; bottom = 1; width = 20; height = 18;
% units:inches|centimeters|normalized|points|{pixels}|characters
set(h_fig, 'Units','centimeters', 'position', [left, bottom, width, height], 'color', 'w');
%% Remove the blank edge
set(gca,'LooseInset',get(gca,'TightInset'));

%% plot Ω based on data
model.set1.data.P_agg = model.set1.data.P_agg(1:model.settings.num_sample,:);
model.set1.data.price = model.set1.data.price(1:model.settings.num_sample,:);
[k,av] = convhull(model.set1.data.P_agg);

% % plot vertices
legend_1_1 = scatter(model.set1.data.P_agg(:,1), model.set1.data.P_agg(:,2),20,[0.9290 0.6940 0.1250]); % Vertice(F)
set(legend_1_1,'MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerEdgeColor',[0.9290 0.6940 0.1250]);
hold on;

% % plot price vector
for i = 1:model.settings.num_sample
    quiver_price = quiver(model.set1.data.P_agg(i,1), model.set1.data.P_agg(i,2), ...
        model.set1.data.price(i, 1)/norm(model.set1.data.price(i, :), 2),...
        model.set1.data.price(i, 2)/norm(model.set1.data.price(i, :), 2), 40, 'b');
    quiver_price.MaxHeadSize = 10;
    hold on;
end
legend_1_2 = quiver_price;

% % plot conv(Gamma)
legend_1_3 = plot(model.set1.data.P_agg(k, 1), model.set1.data.P_agg(k,2), LineWidth=2, Color='r'); hold on;
legend_1_4 = fill(model.set1.data.P_agg(k, 1), model.set1.data.P_agg(k,2),'r'); hold on; % Conv(F)
set(legend_1_4,'edgealpha',0,'facealpha',0.2);

% %  plot D
model.set1.var.D = sdpvar(2, 1);
model.set1.Cons_D = [];
model.set1.Cons_D = model.set1.Cons_D + (...
    model.set1.data.price(k, 1:2)*model.set1.var.D >=...
    sum(model.set1.data.price(k, :).*model.set1.data.P_agg(k, :),2));
v = vertex(model.set1.Cons_D,model.set1.var.D);
[k,av] = convhull(v');
legend_1_5 = plot(v(1,k),v(2,k),'Color', [0 0 0],'LineWidth',0.5); %
legend_1_6 = fill(v(1,k),v(2,k),'k'); hold on; % Conv(F)
set(legend_1_6,'edgealpha',0,'facealpha',0.1);

%% plot Ωphy when N = 1
% % Define var
model.set1.var.P_agg = sdpvar(model.set1.data.num_period, 1, 'full');
model.set1.var.Load_adj = sdpvar(model.set1.data.num_period, model.set1.data.num_adjload, 'full');
model.set1.var.P_vb = sdpvar(model.set1.data.num_period, model.set1.data.num_vb, 'full');
% % Define cons
model.set1.Cons_theta = [];
% Adjustable load power constraints
model.set1.Cons_theta = model.set1.Cons_theta + (...
    model.set1.para.load_adj_low <= model.set1.var.Load_adj <= ...
    model.set1.para.load_adj_up );
% vb constraints
model.set1.Cons_theta = model.set1.Cons_theta + (...
    model.set1.para.p_vb_low <= model.set1.var.P_vb <= ...
    model.set1.para.p_vb_up );
for i=1 : model.set1.data.num_vb
    model.set1.Cons_theta = model.set1.Cons_theta + (...
        model.set1.para.e_vb_min(:,i) <= ...
        model.set1.data.vb.Gamma(:,:,i)*model.set1.var.P_vb(:,i) + ...
        model.set1.para.e_vb_0(i) <= ...
        model.set1.para.e_vb_max(:,i) );
end
% aggregate power
model.set1.Cons_theta = model.set1.Cons_theta + (...
    model.set1.var.P_agg == ...
    sum(model.set1.para.load_fixed,2) + ...
    sum(model.set1.var.Load_adj,2) + ...
    sum(model.set1.var.P_vb,2) );
% plot Ωphy
v = vertex(model.set1.Cons_theta, model.set1.var.P_agg);
[k,av] = convhull(v');
legend_2 = plot(v(1,k),v(2,k),'LineStyle','--', ...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5);
hold on;

%% plot Ωphy when N = 2
% % Define var
model.set2.var.P_agg = sdpvar(model.set2.data.num_period, 1, 'full');
model.set2.var.Load_adj = sdpvar(model.set2.data.num_period, model.set2.data.num_adjload, 'full');
model.set2.var.P_vb = sdpvar(model.set2.data.num_period, model.set2.data.num_vb, 'full');
% % Define cons
model.set2.Cons_theta = [];
% Adjustable load power constraints
model.set2.Cons_theta = model.set2.Cons_theta + (...
    model.set2.para.load_adj_low <= model.set2.var.Load_adj <= ...
    model.set2.para.load_adj_up );
% vb constraints
model.set2.Cons_theta = model.set2.Cons_theta + (...
    model.set2.para.p_vb_low <= model.set2.var.P_vb <= ...
    model.set2.para.p_vb_up );
for i=1 : model.set2.data.num_vb
    model.set2.Cons_theta = model.set2.Cons_theta + (...
        model.set2.para.e_vb_min(:,i) <= ...
        model.set2.data.vb.Gamma(:,:,i)*model.set2.var.P_vb(:,i) + ...
        model.set2.para.e_vb_0(i) <= ...
        model.set2.para.e_vb_max(:,i) );
end
% aggregate power
model.set2.Cons_theta = model.set2.Cons_theta + (...
    model.set2.var.P_agg == ...
    sum(model.set2.para.load_fixed,2) + ...
    sum(model.set2.var.Load_adj,2) + ...
    sum(model.set2.var.P_vb,2) );
% % calculate the vertices of the full space
[KKTConstraints, details] = kkt(model.set2.Cons_theta, sum(model.set2.var.P_agg(:)), model.set2.var.P_agg);
id_var = [];
for k = 1 : size(details.f, 1)
    if isa(details.f(k), 'sdpvar')
        id_var = [id_var k];
    end
end
id_nonvar = setdiff(1:size(details.f,1), id_var);
[V,nr,nre] = lcon2vert(full(details.A), details.b, ...
    full(details.E(id_nonvar,:)), details.f(id_nonvar), 1e-8);
% % project on the aggregate vars
V_agg = V * details.E(id_var,:)';
% % get the vertices of the aggregate vars
k = convhull(V_agg);
% % plot Ωphy
legend_3= plot(V_agg(k,1),V_agg(k,2), ...
    'LineStyle','--','LineWidth',1.5);   % Ωphy

%% set asix
axis([-100 700 -200 900]);
set(h_axis, 'XTick', 0:200:600);                        % or xticks();
set(h_axis, 'YTick', -200:200:800);                     % or xticks();

%% lengend
mylegend = legend([legend_1_1, legend_1_2, legend_1_4, legend_1_6, legend_2, legend_3],...
    {'Vertice($\Gamma$)',...
    'Price vector',...
    '$\rm{Conv}(\Gamma)$',...
    'D',...
    '$\partial \Omega_{phy}$ ($\rm{N_{vb}}=1$)',...
    '$\partial \Omega_{phy}$ ($\rm{N_{vb}}=2$)'},...
    'Orientation','horizontal','NumColumns',2, ...
    'Interpreter','latex');

% Create xlabel
xlabel('$P^1$ (kW)', 'Interpreter','latex');
% Create ylabel
ylabel('$P^2$ (kW)', 'Interpreter','latex');
%
hold on; grid on;box on;

%% Font
set(h_axis, 'FontName', 'Times New Roman', 'FontSize', 28);
set(mylegend, 'FontSize', 24, 'Location', 'north', 'Color','none', 'EdgeColor','none');

%% style
set(h_axis,'Color','none');

%% save
filepath = [cd '\results\'];
str = [filepath 'Fig_2'];
saveas(gca, str, 'fig');
saveas(gca, str, 'svg');
saveas(gca, str, 'emf');