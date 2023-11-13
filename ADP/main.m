clc
clear
close all
%%Index setting
% bus idx
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
% branch idx
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
% gen idx
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
% cost idx
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

mpc = loadcase('case33mg');
%mpc = loadcase('case69');
mpc         = runopf(mpc);

id_bus      = mpc.bus(:,BUS_I);
id_gen      = mpc.gen(:,GEN_BUS);
Nbus        = size(mpc.bus,1);
Ngen        = numel(id_gen);
Nbranch     = size(mpc.branch,1);

entries_pf{1} = 1:Nbus;                        % vmag
entries_pf{2} = (Nbus+1):2*Nbus;               % vang
entries_pf{3} = (2*Nbus+1):(2*Nbus+Ngen);      % Pg
entries_pf{4} = (2*Nbus+Ngen+1):2*(Nbus+Ngen); % Qg 

baseMVA     = mpc.baseMVA;              % baseMVA
cost_param  = mpc.gencost(:,5:end);     % objective coefficients
vmax        = mpc.bus(:,VMAX);                                            
vmin        = mpc.bus(:,VMIN);
Pgmin       = mpc.gen(:,PMIN)/baseMVA;
Qgmin       = mpc.gen(:,QMIN)/baseMVA;
Pgmax       = mpc.gen(:,PMAX)/baseMVA;
Qgmax       = mpc.gen(:,QMAX)/baseMVA;
Fmax        = mpc.branch(:,RATE_A)/baseMVA;

Pd          = mpc.bus(:,PD)/baseMVA;   
Qd          = mpc.bus(:,QD)/baseMVA;
Cg          = sparse(id_gen,1:Ngen,ones(Ngen,1),Nbus,Ngen);

%% lower & upper bounds
lbx         = [vmin;-inf*ones(Nbus,1);Pgmin;Qgmin];       
ubx         = [vmax; inf*ones(Nbus,1);Pgmax;Qgmax];
%% initial state x0
vang0       = mpc.bus(:,VA)/180*pi;
vmag0       = mpc.bus(:,VM);
U0          = vmag0.*cos(vang0);
W0          = vmag0.*sin(vang0);
Pg0         = mpc.gen(:,PG)/baseMVA;
Qg0         = mpc.gen(:,QG)/baseMVA;
x0      = vertcat(vmag0, vang0, Pg0, Qg0);
%% equality & inequality constraints - current balance constraints
% create Ybus Yf Yt
[Ybus, Yf, Yt] = makeYbus(mpc);
Gbus           = real(Ybus);
Bbus           = imag(Ybus);
Gf             = real(Yf);
Bf             = imag(Yf);
Gt             = real(Yt);
Bt             = imag(Yt);
from_bus       = mpc.branch(:, F_BUS);                           %% list of "from" buses
to_bus         = mpc.branch(:, T_BUS);  
Cf             = sparse(1:Nbranch,from_bus,ones(Nbranch,1),Nbranch,Nbus);
Ct             = sparse(1:Nbranch,to_bus,ones(Nbranch,1),Nbranch,Nbus);
% slack
slack_bus_entries =  find(mpc.bus(:,BUS_TYPE) == REF);
idx_ref = entries_pf{2}(slack_bus_entries);
eq_ref = @(x)x(idx_ref);
% power flow equation
eq_pf          = @(x)create_local_power_flow_equation_pol(x(entries_pf{1}),x(entries_pf{2}),...
    x(entries_pf{3}),x(entries_pf{4}),Gbus,Bbus,Pd,Qd,Cg);
test_pf        = eq_pf(x0);
Npf            = numel(eq_pf(x0));
% upper & lower bound for voltage magnitude
ineq_voltage   = @(x) x(entries_pf{1});
% line Limit
ineq_line = [];
idx_limit = [];
Nlimit   = 0;
g   = @(x)vertcat(eq_pf(x),eq_ref(x));

% lbx = vertcat(vmin,)
% ubx = vertcat(vmax,)
lbg = vertcat(vmin, -inf*ones(Nlimit,1), zeros(Npf+1,1));
ubg = vertcat(vmax, zeros(Npf+Nlimit+1,1));
%% solver options
import casadi.*
% tolerance
tol        = 1e-8;
options.ipopt.tol             = tol;
options.ipopt.constr_viol_tol = tol;
options.ipopt.compl_inf_tol   = tol;
options.ipopt.acceptable_tol  = tol;
options.ipopt.acceptable_constr_viol_tol = tol;
options.ipopt.print_level = 5;
% options.ipopt.grad_f = fgrad;
options.print_time        = 5;
options.ipopt.max_iter    = 1000;

Nx         =  numel(x0);
x_SX       =   SX.sym('x',Nx,1);
constraint = g(x_SX);
i = 1;
Points = zeros(8,2);
for c1 = -1:1
    for c2 = -1:1
        if c1== 0 && c2 ==0
            continue;
        else
            f   = @(x) c1*x(entries{3}(1))+c2*x(entries{4}(1));
            nlp = struct('x',x_SX,'f',objective,'g',constraint);
            S   = nlpsol('solver','ipopt', nlp,options);
            sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
                    'lbx', lbx, 'ubx', ubx);
            xopt= full(sol.x);
            Points(i,:) = [xopt(2*Nbus+1),xopt(2*Nbus+Ngen)];
            i = i+1;
        end
    end
end