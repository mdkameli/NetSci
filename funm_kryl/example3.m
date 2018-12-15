% non-Hermitian case
% Compute expm(A)*b using the restarted Krylov method FUNM_KRYL
% please check index.htm for further details

randn('state',0);
rand('state',0);
A = rand(500,500);
[A,dummy] = qr(A);
A = .8*A;
b = randn(500,1);
exact = expm(A)*b;
b = b/norm(exact); exact = exact/norm(exact);


% Use expm at the core of the restart algorithm
% which becomes slower with each cycle

param.function = @expm;       % other choices: 'expBA', 'expCF', ...
param.restart_length = 5;
param.max_restarts = 50;
param.hermitian = 0;          % set 0 if A is not Hermitian
param.V_full = 1;             % set 1 if you need Krylov basis
param.H_full = 1;             % if using rational functions you can set this 0
param.exact = exact;          % if not known set to []
param.bound = 1;              % returns upper and lower bounds (after some cycles)
param.stopping_accuracy = 1e-16;  % stopping accuracy
param.inner_product = @inner_product;
param.thick = [];             % thick-restart function  
param.min_decay = 0.95;       % we desire linear error reduction of rate < .95 
param.waitbar = 1;            % show waitbar 
param.reorth_number = 0;      % #reorthogonalizations
param = param_init(param);    % check and correct param structure

[f,out1] = funm_kryl(A,b,param);


% Use best rational approximation to exp on unit disk at the core
% of the restart algorithm -> constant work per cycle

param.function = parfrac('cfud',8,@exp);
param.H_full = 0;
param = param_init(param);

[f,out2] = funm_kryl(A,b,param);


% ---------------------------------
% Plot results

report(param,out1);
set(gcf,'Name','using expm');
report(param,out2);
set(gcf,'Name','using rational approximation');
