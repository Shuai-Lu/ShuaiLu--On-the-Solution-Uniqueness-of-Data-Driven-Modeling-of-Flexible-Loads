%% generate (price, power) data
%  2024-02-21
%  By Jiayi Ding & Shuai Lu

%% Initialize
clc;clear;close all;
warning off;
yalmip('clear');

%% Display current time
StartTime = datestr(now,31); % Start time
fprintf('%s\n', ['- Start time: ' StartTime]);

%% Declaring Global Variables
% global data model CodeTime;
CodeTime.Overall = tic; % start stopwatch timer
update_random = 0;

%% Param Settings
% Initialize data structure
data = [];
% Set the number of periods
data.para.num_period = 2;
% Set the delta period
data.para.delta_period = 1;
% Set the number of fixed loads
data.para.num_fixedload = 1; 
% Set the number of adjustable loads
data.para.num_adjload = 1;   
% Set the number of virtual battery (vb)
data.para.num_vb = 4;        
% Set the number of samples
data.para.num_sample = 200;  
% Set the fixed load values
data.para.Load_fix = [120 150];
% Set the lower limit for adjustable loads
data.para.Load_adj_low = 10;
% Set the upper limit for adjustable loads
data.para.Load_adj_up = 70;
% Set the lower power limits for virtual battery
data.para.P_vb_low = [-95;-75;-60;-20];
% Set the upper power limits for virtual battery
data.para.P_vb_up = [40;140;120;40];
% Set the lower energy limits for virtual battery
data.para.E_vb_low = [40;80;120;240];
% Set the upper energy limits for virtual battery
data.para.E_vb_up = [180;250;300;450];
% Set the reliability coefficients for virtual battery
data.para.sigma = [0.95; 0.90; 0.85; 0.88];
% Set the initial energy level for virtual battery
data.para.E_vb_0 = data.para.E_vb_low;

% Calculate Gamma Matrix
data.para.Gamma = zeros(data.para.num_period, data.para.num_period, data.para.num_vb);
for n=1:data.para.num_vb
    for i=1:data.para.num_period
        for j=1:data.para.num_period
            if i<j
                data.para.Gamma(i,j,n)=0;
            else
                data.para.Gamma(i,j,n)=data.para.sigma(n)^(i-j);
            end
        end
    end
end

% Randomly generated price signals
filename = 'random_price_benchmark.mat'; % set to 'random_price_benchmark.mat' to reproduce the results in paper
filepath = [cd '\data\' filename];
if ~ exist(filepath) || update_random
    random_price = rand(data.para.num_sample, data.para.num_period)*2-1; % [-1, 1]
    data.para.price = random_price;
    save(filepath, "random_price");
else
    temp = load(filepath);
    data.para.price = temp.random_price;
end


%% Calculate the optimal power under different prices
%% Initialize the model
model = [];
model.var = [];
model.cons = [];
model.obj = [];

%% Define vars
% Define adjustable load variable (num_adjload*num_period)
model.var.Load_adj = sdpvar(data.para.num_adjload,data.para.num_period,'full');
% Define the vb power variable (num_vb*num_period)
model.var.P_vb = sdpvar(data.para.num_vb,data.para.num_period,'full');
% Define the vb storage capacity variable (num_vb*num_period)
model.var.E_vb = sdpvar(data.para.num_vb,data.para.num_period,'full');
% Define the aggregate power variable (1*num_period)
model.var.P_agg = sdpvar(1,data.para.num_period,'full');


%% Define cons
% % load_adj
model.cons = model.cons + (...
    (data.para.Load_adj_low <= model.var.Load_adj <= ...
    data.para.Load_adj_up): ...
    'Adjustable load power constraints' );
% % virtual battery
model.cons = model.cons + (...
    (data.para.P_vb_low*ones(1,data.para.num_period) <= ...
    model.var.P_vb <= ...
    data.para.P_vb_up*ones(1,data.para.num_period)): ...
    'vb power constraints' );
for i=1:data.para.num_vb
    model.cons = model.cons + (...
        (data.para.E_vb_low(i) <= ...
        data.para.Gamma(:,:,i)*model.var.P_vb(i,:)'*data.para.delta_period + ...
        data.para.E_vb_0(i) <= ...
        data.para.E_vb_up(i)): ...
        'vb storage energy constraint' );
end
% % Aggregate power
model.cons = model.cons + (...
    (model.var.P_agg == sum(data.para.Load_fix,1) + ...
    sum(model.var.Load_adj,1)+sum(model.var.P_vb,1)): ...
    'Aggregate power' );

%% Solve the model
for i = 1 : data.para.num_sample
    if rem(i,10) == 0
        fprintf('%s%d\n', ' - Solve the model ', i);
    end
    model.obj = data.para.price(i,:) * model.var.P_agg'; % obj
    model.ops = sdpsettings('solver', 'gurobi','verbose',0);  % 用gurobi求解器
    model.sol = optimize(model.cons, model.obj, model.ops);
    if model.sol.problem == 0  
        data.dataset.P_agg(i,:) = value(model.var.P_agg);
    else
        disp('Hmm, something went wrong!');
        model.sol.info
        yalmiperror(model.sol.problem)
    end
end

% %
data.dataset.price = data.para.price;

%% Save the data
filepath = [cd '\data\data.mat'];
save(filepath, "data");

%% Plot the results
h_fig = figure();          % gcf: get current figure
h_axis = gca;              % gca: get current axis
% set position & color
left = 1; bottom = 1; width = 20; height = 18;
% units:inches|centimeters|normalized|points|{pixels}|characters
set(h_fig, 'Units','centimeters', 'position', [left, bottom, width, height], 'color', 'w');
% Remove the blank edge
set(gca,'LooseInset',get(gca,'TightInset'));

% % get the P_agg data
data.para.P_temp = data.dataset.P_agg(1:data.para.num_sample,:);
[k,av] = convhull(data.para.P_temp);

% % plot vertices of conv(Gamma)
legend1 = scatter(data.para.P_temp(:,1), data.para.P_temp(:,2),20,[0.9290 0.6940 0.1250]); % Vertice(F)
set(legend1 ,'MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerEdgeColor',[0.9290 0.6940 0.1250]);
hold on;

% % plot price vector
for i=1:data.para.num_sample
    quiver_price = quiver(data.para.P_temp(i,1), data.para.P_temp(i,2), ... 
        data.para.price(i, 1)/norm(data.para.price(i, :), 2), ...
        data.para.price(i, 2)/norm(data.para.price(i, :), 2), 40, 'b');
    quiver_price.MaxHeadSize = 10;
    hold on;
end
legend2 = quiver_price;

% % plot conv(Gamma) 
plot(data.para.P_temp(k, 1), data.para.P_temp(k,2), LineWidth=2, Color='r'); hold on;
legend3 = fill(data.para.P_temp(k, 1), data.para.P_temp(k,2),'r'); % Conv(F)
hold on;
set(legend3,'edgealpha',0,'facealpha',0.2);

% % plot D
D = sdpvar(2, 1);
Constraints_D = [];
Constraints_D = Constraints_D + (...
    data.para.price(k, 1:2)*D >= ...
    sum(data.para.price(k, 1:2).*data.dataset.P_agg(k, 1:2),2));
v = vertex(Constraints_D,D);
[k,av] = convhull(v');
plot(v(1,k),v(2,k),'Color', [0 0 0],'LineWidth',0.5);
legend4 = fill(v(1,k),v(2,k),'k'); %
hold on;
set(legend4,'edgealpha',0,'facealpha',0.1);

% set asix
axis([-100 700 -100 800]);
set(h_axis, 'XTick', 0:200:800);                        % or xticks();
set(h_axis, 'YTick', -200:200:1000);                        % or xticks();

% lengend
mylegend = legend([legend1, legend2, legend3, legend4],...
    {'Vertice($\Gamma$)',...
    'Price vector',...
    '$\rm{Conv}(\Gamma)$',...
    'D'},...
    'Orientation','horizontal','NumColumns',2, ...
    'Interpreter','latex');

% Create xlabel
xlabel('$P^1$ (kW)', 'Interpreter','latex');
% Create ylabel
ylabel('$P^2$ (kW)', 'Interpreter','latex');
hold on; grid on;box on;

% Font
set(h_axis, 'FontName', 'Times New Roman', 'FontSize', 28);
set(mylegend, 'FontSize', 24, 'Location', 'north', 'Color','none', 'EdgeColor','none');

% style
set(h_axis,'Color','none');

% save
filepath = [cd '\results\'];
saveas(gca, [filepath 'Samples.fig']);
saveas(gca, [filepath 'Samples.svg']);
saveas(gca, [filepath 'Samples.emf']);
hold off;

%%
CodeTime.Overall = toc(CodeTime.Overall);
EndTime = datestr(now,31);
fprintf('%s\n', ['- End time: ' EndTime]);
fprintf('%s%.2f%s\n', '- Total time: ', CodeTime.Overall, '(s)');




