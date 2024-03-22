%% identify the Physics-Based Model
%  dependencies: Yalmip, Gurobi
%  2024-02-21
%  By Jiayi Ding & Shuai Lu

%% Initialize
clc; clear; close all;
warning off;
yalmip('clear');

%% set num_vb
num_vb = 1; % set to 1 for Ωphy^1 and 2 for Ωphy^2

%% load data
filepath = [cd '\data\data.mat'];
temp = load(filepath);
data = temp.data;

%% Display current time
StartTime = datestr(now, 31); % Start time
fprintf('%s\n', ['- Start time: ' StartTime]);

%% Declaring Global Variables
% global model CodeTime;
CodeTime.Overall = tic; % start stopwatch timer
update_random = 0;

%% Param Settings
model.settings.num_sample = 50; % this parameter can be adjusted
model.settings.big_M = 1e3;
model.para.num_fixedload = 1;
model.para.num_adjload = 1;
model.para.num_vb = num_vb;
model.para.price = data.dataset.price(1:model.settings.num_sample, :);
model.para.P_agg = data.dataset.P_agg(1:model.settings.num_sample, :);
model.para.num_period = data.para.num_period;

% Paramteres of virtual battery
model.para.vb.sigma(:, 1) = [1; 0.95]; % the length need to be equal to num_vb
model.para.vb.Gamma = zeros(data.para.num_period, data.para.num_period, model.para.num_vb);
for n = 1:model.para.num_vb
    for i = 1:data.para.num_period
        for j = 1:i
            model.para.vb.Gamma(i, j, n) = model.para.vb.sigma(n) ^ (i - j);
        end
    end
end

fprintf('%s\n', '- Model the forward problem ...');
%% get the concise form based on the kkt function
model.var = [];
model.cons = [];
% % Define vars
% parameters
id_P_fixed = 1:model.para.num_fixedload;
id_P_adjload_LB = id_P_fixed(end) + 1:id_P_fixed(end) + model.para.num_adjload;
id_P_adjload_UB = id_P_adjload_LB(end) + 1:id_P_adjload_LB(end) + model.para.num_adjload;
id_P_vb_LB = id_P_adjload_UB(end) + 1:id_P_adjload_UB(end) + model.para.num_vb;
id_P_vb_UB = id_P_vb_LB(end) + 1:id_P_vb_LB(end) + model.para.num_vb;
id_E_vb_min = id_P_vb_UB(end) + 1:id_P_vb_UB(end) + model.para.num_vb;
id_E_vb_max = id_E_vb_min(end) + 1:id_E_vb_min(end) + model.para.num_vb;
id_E_vb_initial = id_E_vb_max(end) + 1:id_E_vb_max(end) + model.para.num_vb;
% % define vars
model.var.parameter = sdpvar(data.para.num_period, id_E_vb_initial(end), 'full');
% price
model.var.price = rand(data.para.num_period, 1);
% fixed load
model.var.P_fixed = sdpvar(data.para.num_period, 1, 'full');
% adjustable load
model.var.P_adjload_LB = sdpvar(data.para.num_period, model.para.num_adjload, 'full');
model.var.P_adjload_UB = sdpvar(data.para.num_period, model.para.num_adjload, 'full');
% virtual battery
model.var.P_vb_LB = sdpvar(data.para.num_period, model.para.num_vb, 'full');
model.var.P_vb_UB = sdpvar(data.para.num_period, model.para.num_vb, 'full');
model.var.E_vb_min = sdpvar(data.para.num_period, model.para.num_vb, 'full');
model.var.E_vb_max = sdpvar(data.para.num_period, model.para.num_vb, 'full');
model.var.E_vb_initial = sdpvar(data.para.num_period, model.para.num_vb, 'full');
% power
model.var.P_adjload = sdpvar(data.para.num_period, model.para.num_adjload, 'full');
model.var.P_vb = sdpvar(data.para.num_period, model.para.num_vb, 'full');
model.var.E_vb = sdpvar(data.para.num_period, model.para.num_vb, 'full');
model.var.P_agg = sdpvar(data.para.num_period, 1, 'full');
% intermediate variable
model.var.temp_E_vb = sdpvar(data.para.num_period, model.para.num_vb, 'full');

% % link parameters and vars
model.cons = model.cons + ( ...
    model.var.parameter(:, id_P_fixed) == model.var.P_fixed);
model.cons = model.cons + ( ...
    model.var.parameter(:, id_P_adjload_LB) == model.var.P_adjload_LB);
model.cons = model.cons + ( ...
    model.var.parameter(:, id_P_adjload_UB) == model.var.P_adjload_UB);
model.cons = model.cons + ( ...
    model.var.parameter(:, id_P_vb_LB) == model.var.P_vb_LB);
model.cons = model.cons + ( ...
    model.var.parameter(:, id_P_vb_UB) == model.var.P_vb_UB);
model.cons = model.cons + ( ...
    model.var.parameter(:, id_E_vb_min) == model.var.E_vb_min);
model.cons = model.cons + ( ...
    model.var.parameter(:, id_E_vb_max) == model.var.E_vb_max);
model.cons = model.cons + ( ...
    model.var.parameter(:, id_E_vb_initial) == model.var.E_vb_initial);

% % Define cons
% aggregate power
model.cons = model.cons + ( ...
    model.var.P_agg(:, 1) == model.var.P_fixed(:, 1) + ...
    sum(model.var.P_adjload(:, :), 2) + ...
    sum(model.var.P_vb(:, :), 2));

% power of adjustable load
model.cons = model.cons + ( ...
    model.var.P_adjload_LB(:, :) <= model.var.P_adjload(:, :) <= ...
    model.var.P_adjload_UB(:, :));

% power & energy of vb
model.cons = model.cons + ( ...
    model.var.P_vb_LB(:, :) <= model.var.P_vb(:, :) <= ...
    model.var.P_vb_UB(:, :));

for n = 1:model.para.num_vb
    model.var.temp_E_vb(:, n) = ...
        model.para.vb.Gamma(:, :, n) * model.var.P_vb(:, n) * ...
        data.para.delta_period + ...
        model.var.E_vb_initial(:, n);
end

model.cons = model.cons + ( ...
    model.var.E_vb_min <= model.var.temp_E_vb(:, :) <= ...
    model.var.E_vb_max);

% objective
model.obj = model.var.price' * model.var.P_agg(:);

% % derive the concise form based on kkt function
fprintf('%s\n', '- Derive the concise form ...');
model.ops = sdpsettings('kkt.dualbounds', 0, 'verbose', 0);
[model.cons_kkt, model.details] = ...
    kkt(model.cons, model.obj, model.var.parameter, model.ops);

%% model the inverse optimization problem
fprintf('%s\n', '- Model the invesre optimization problem ...');
model.inverseproblem = [];
num_dual_c = size(model.details.c, 1);
num_dual_eq = size(model.details.f, 1);
num_dual_ineq = size(model.details.b, 1);

% % parameters
model.inverseproblem.var.parameters = sdpvar(data.para.num_period, id_E_vb_initial(end), 'full');

% % primal vars
model.inverseproblem.var.x = sdpvar(num_dual_c, model.settings.num_sample, 'full');
model.inverseproblem.var.P_agg = sdpvar(data.para.num_period, model.settings.num_sample, 'full');
model.inverseproblem.var.delta = sdpvar(data.para.num_period, model.settings.num_sample, 'full');

% % dual vars
model.inverseproblem.var.c = sdpvar(num_dual_c, model.settings.num_sample, 'full');
model.inverseproblem.var.x = sdpvar(num_dual_c, model.settings.num_sample, 'full');
model.inverseproblem.var.dual_eq = sdpvar(num_dual_eq, model.settings.num_sample, 'full');
model.inverseproblem.var.dual_ineq = sdpvar(num_dual_ineq, model.settings.num_sample, 'full');
model.inverseproblem.var.dual_ineq_binary = binvar(num_dual_ineq, 1, model.settings.num_sample, 'full');

% % link paramters
id_f_sdpvar = [];
for i = 1:size(model.details.f, 1)
    if isa(model.details.f(i), 'sdpvar')
        id_f_sdpvar = [id_f_sdpvar; i];
    end

end

model.details.f(id_f_sdpvar) = reshape(model.inverseproblem.var.parameters, [], 1);

% % cons
model.inverseproblem.cons = [];

% % parameter
model.inverseproblem.cons = model.inverseproblem.cons + ( ...
    ones(data.para.num_period - 1, 1) * model.inverseproblem.var.parameters(1, id_E_vb_initial) == ...
    model.inverseproblem.var.parameters(2:end, id_E_vb_initial));

% % primal cons
% Loop through all samples
k = 1:model.settings.num_sample;
% Add inequality constraint: A * x <= b for each sample
model.inverseproblem.cons = model.inverseproblem.cons + ( ...
    model.details.A * model.inverseproblem.var.x(:, k) <= model.details.b * ones(1, model.settings.num_sample));
% Add equality constraint: E * x == f for each sample
model.inverseproblem.cons = model.inverseproblem.cons + ( ...
    model.details.E * model.inverseproblem.var.x(:, k) == model.details.f * ones(1, model.settings.num_sample));
% Add constraint to relate the aggregate variable P_agg with x for each sample
model.inverseproblem.cons = model.inverseproblem.cons + ( ...
    model.inverseproblem.var.P_agg(:, k) == ...
    model.inverseproblem.var.x(end - model.para.num_period + 1:end, k));

% % dual cons
% set the dual variables for each sample
for k = 1:model.settings.num_sample
    model.inverseproblem.var.dual_c(:, k) = [model.details.c(1:end - 2); model.para.price(k, :)'];
end

% Loop through all samples
k = 1:model.settings.num_sample;
% Add complementary slackness constraint
model.inverseproblem.cons = model.inverseproblem.cons + ( ...
    model.inverseproblem.var.dual_c(:, k)' + ...
    model.inverseproblem.var.dual_ineq(:, k)' * model.details.A + ...
    model.inverseproblem.var.dual_eq(:, k)' * model.details.E == 0);
% Big-M constraints for inequality constraints
model.inverseproblem.cons = model.inverseproblem.cons + ( ...
    0 <= (model.details.b * ones(1, model.settings.num_sample) - model.details.A * model.inverseproblem.var.x(:, k)) <= ...
    model.settings.big_M * (1 - model.inverseproblem.var.dual_ineq_binary(:, k)));
% Bound constraints for equality dual variables
model.inverseproblem.cons = model.inverseproblem.cons + ( ...
    - model.settings.big_M <= model.inverseproblem.var.dual_eq(:, k) <= model.settings.big_M);
% Non-negativity and bounding constraints for inequality dual variables
model.inverseproblem.cons = model.inverseproblem.cons + ( ...
    0 <= model.inverseproblem.var.dual_ineq(:, k) <= ...
    model.settings.big_M * model.inverseproblem.var.dual_ineq_binary(:, k));

% % obj
model.inverseproblem.cons = model.inverseproblem.cons + ( ...
    model.inverseproblem.var.delta >= ...
    model.para.P_agg(1:model.settings.num_sample, :)' - model.inverseproblem.var.P_agg);
model.inverseproblem.cons = model.inverseproblem.cons + ( ...
    model.inverseproblem.var.delta >= ...
    - model.para.P_agg(1:model.settings.num_sample, :)' + model.inverseproblem.var.P_agg);
model.inverseproblem.obj = sum(model.inverseproblem.var.delta(:));

% model.inverseproblem.obj = sum(model.inverseproblem.var.delta(:).^2); % quadratic objective is time consuming

% % solve
fprintf('%s\n', '- Solve the inverse optimization problem ...');
ops = sdpsettings('solver', 'gurobi', 'verbose', 2);
model.inverseproblem.sol = optimize(model.inverseproblem.cons, model.inverseproblem.obj, ops);

%% get the results
fprintf('%s\n', '- Obtain the results ...');
%% Solve
if model.inverseproblem.sol.problem == 0 % Extract and display value
    model.inverseproblem.var = myFun_GetValue(model.inverseproblem.var);
    model.results.load_fixed = model.inverseproblem.var.parameters(:, id_P_fixed); % model.var.P_adjload_LB  = 10
    model.results.load_adj_low = model.inverseproblem.var.parameters(:, id_P_adjload_LB); % model.var.P_adjload_LB  = 10;
    model.results.load_adj_up = model.inverseproblem.var.parameters(:, id_P_adjload_UB); % model.var.P_adjload_UB = 50;
    model.results.p_vb_low = model.inverseproblem.var.parameters(:, id_P_vb_LB); % model.var.P_vb_LB = -80;
    model.results.p_vb_up = model.inverseproblem.var.parameters(:, id_P_vb_UB); % model.var.P_vb_UB = 40;
    model.results.e_vb_max = model.inverseproblem.var.parameters(:, id_E_vb_max); % model.var.E_vb_max = 180;
    model.results.e_vb_min = model.inverseproblem.var.parameters(:, id_E_vb_min); % model.var.E_vb_min = 40;
    model.results.e_vb_0 = model.inverseproblem.var.parameters(1, id_E_vb_initial); % model.var.E_vb_initial = 40;
    model.results.loss = value(model.inverseproblem.obj);
else
    disp('Hmm, something went wrong!');
    model.inverseproblem.sol.info
    yalmiperror(model.inverseproblem.sol.problem)
end

%% save
fprintf('%s\n', '- Save the data ...');

%% Save
filepath = [cd '\data\identification_' num2str(model.para.num_vb), '_', ...
                num2str(model.settings.num_sample) '.mat'];
save(filepath, 'model');

%% End
CodeTime.Overall = toc(CodeTime.Overall);
EndTime = datestr(now, 31);
fprintf('%s\n', ['- End time: ' EndTime]);
fprintf('%s%.2f%s\n', '- Total time: ', CodeTime.Overall, '(s)');
