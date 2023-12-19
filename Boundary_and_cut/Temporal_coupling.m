%% Only consider 6 time periods
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

mpc = ext2int(loadcase('case33_modified'));
mpc = ext2int(mpc);
% load('Pd2_test.mat');
% load("Qd2_test.mat");
%% parameters
baseMVA     = mpc.baseMVA;
baseKV      = mpc.bus(1,BASE_KV);
Umax        = mpc.bus(:,VMAX).^2;                                            
Umin        = mpc.bus(:,VMIN).^2;
Pgmin       = mpc.gen(:,PMIN)/baseMVA;  
Qgmin       = mpc.gen(:,QMIN)/baseMVA; 
Pgmax       = mpc.gen(:,PMAX)/baseMVA; 
Qgmax       = mpc.gen(:,QMAX)/baseMVA; 
Pramp       = [1;1;1];
Pd          = mpc.bus(:,PD)/baseMVA; Pd = repmat(Pd,1,6);%Pd=[Pd Pd2];
Qd          = mpc.bus(:,QD)/baseMVA; Qd = repmat(Qd,1,6);%Qd=[Qd Qd2];

id_gen      = mpc.gen(:,GEN_BUS);
id_slack    =  find(mpc.bus(:,BUS_TYPE) == REF);
id_gen_slack  = find(id_gen == id_slack);
id_gen_nslack = find(id_gen ~= id_slack);
Ngen_nslack = numel(id_gen_nslack);
Nbus        = size(mpc.bus,1);
Ngen        = numel(id_gen);
Nbranch     = size(mpc.branch,1);

branch_r   = mpc.branch(:,BR_R);
branch_x   = mpc.branch(:,BR_X);

from_bus       = mpc.branch(:, F_BUS);                         
to_bus         = mpc.branch(:, T_BUS);  
Cf             = sparse(1:Nbranch,from_bus,ones(Nbranch,1),Nbranch,Nbus);
Ct             = sparse(1:Nbranch,to_bus,ones(Nbranch,1),Nbranch,Nbus);
C              = Cf - Ct;
Cg             = sparse(id_gen,1:Ngen,ones(Ngen,1),Nbus,Ngen);
Cg_nslack      = Cg(:,id_gen_nslack);
Cg_slack       = Cg(:,id_gen_slack);
%% Problem formulation
% Time period
T = 6;
% beta_V
A_v = [eye(Nbus*T) zeros(Nbus*T,2*T*Ngen+2*T*Nbranch);
       -eye(Nbus*T) zeros(Nbus*T,2*T*Ngen+2*T*Nbranch)];
% beta_Generator at non-slack bus 
A_gen_ns = [zeros(Ngen*T,Nbus*T) eye(Ngen*T) zeros(Ngen*T,Ngen*T+2*T*Nbranch);
            zeros(Ngen*T,Nbus*T) -eye(Ngen*T) zeros(Ngen*T,Ngen*T+2*T*Nbranch);
            zeros(Ngen*T,Nbus*T+Ngen*T) eye(Ngen*T) zeros(Ngen*T,2*T*Nbranch);
            zeros(Ngen*T,Nbus*T+Ngen*T) -eye(Ngen*T) zeros(Ngen*T,2*T*Nbranch);];% bounds for P Q
% temporal coupling constraints
C_diff = zeros(1,Ngen+1);
C_diff(1,1) = -1;
C_diff(1,end) = 1;
C_tc = zeros((T-1)*Ngen,Ngen*T);% temporal coupling
for i = 1:(T-1)*Ngen
    C_tc(i,i:i+Ngen) = C_diff;
end
C_tc = vertcat(C_tc,-C_tc);
A_tc = [zeros(2*(T-1)*Ngen,Nbus*T),C_tc,zeros(size(C_tc)),zeros(2*(T-1)*Ngen,2*T*Nbranch)];% no ramp for Q
% beta_voltage_constraints for injection power
C_comb = C;
BR_comb = diag(branch_r);
BX_comb = diag(branch_x);
for i = 2:T
    C_comb = blkdiag(C_comb,C);
    BR_comb = blkdiag(BR_comb,diag(branch_r));
    BX_comb = blkdiag(BX_comb,diag(branch_x));
end
A_inj = [C_comb zeros(Nbranch*T,2*T*Ngen) -2*BR_comb -2*BX_comb;
         -C_comb zeros(Nbranch*T,2*T*Ngen) 2*BR_comb  2*BX_comb];
% beta_power flow balance
Cg_comb = Cg;
C_comb2 = C';
for i = 2:T
    Cg_comb = blkdiag(Cg_comb,Cg);
    C_comb2 = blkdiag(C_comb2,C');
end
A_eq = [zeros(Nbus*T) Cg_comb zeros(Nbus*T, Ngen*T) -C_comb2 zeros(Nbus*T, Nbranch*T);
        zeros(Nbus*T) -Cg_comb zeros(Nbus*T, Ngen*T) C_comb2 zeros(Nbus*T, Nbranch*T);
        zeros(Nbus*T) zeros(Nbus*T, Ngen*T) Cg_comb zeros(T*Nbus, Nbranch*T) -C_comb2;
        zeros(Nbus*T) zeros(Nbus*T, Ngen*T) -Cg_comb zeros(T*Nbus, Nbranch*T) C_comb2;];

% Integrate all matrices
A = vertcat(A_v,A_gen_ns,A_tc,A_inj,A_eq);

% lower and upper bounds on U
b0_U_u = repmat(Umax,T,1);
b0_U_l = -repmat(Umin,T,1);
% lower and upper bounds on P/Q
b0_P_u = repmat(Pgmax,T,1);
b0_P_l = -repmat(Pgmin,T,1);
b0_Q_u = repmat(Qgmax,T,1);
b0_Q_l = -repmat(Qgmin,T,1);
% ramp rate
b0_ramp = repmat(Pramp,2*(T-1),1);
% power flow equation
tol_cons = 0;
b0_inj = tol_cons*ones(2*Nbranch*T,1);
Pd = Pd(:);
Qd = Qd(:);
b0_pf = vertcat(Pd+tol_cons,-Pd+tol_cons,Qd+tol_cons,-Qd+tol_cons);

b0 = vertcat(b0_U_u,b0_U_l,b0_P_u,b0_P_l,b0_Q_u,b0_Q_l,b0_ramp,b0_inj,b0_pf);
%% condense model with umbrella constraints
[A_u, b_u] = E_UCI(A, b0);
% parameters for projection varaibles
idx_pqs = linspace(1, 1+Ngen*(T-1),T);
B = A_u(:, [Nbus*T+idx_pqs Nbus*T+Ngen*T+idx_pqs]);
% B_p = A_u(:, Nbus*T+idx_pqs);
% B_q = A_u(:, Nbus*T+Ngen*T+idx_pqs);
remainingIndices = setdiff(1:size(A_u,2), [Nbus*T+idx_pqs Nbus*T+Ngen*T+idx_pqs]);
A_de = A_u(:, remainingIndices);


% initilization 
D0 = [ eye(T) zeros(T);
     -eye(T) zeros(T);
      zeros(T) eye(T);
      zeros(T) -eye(T);];
v = [Pgmax(id_gen_slack)*ones(T,1);
     -Pgmin(id_gen_slack)*ones(T,1);
     Qgmax(id_gen_slack)*ones(T,1);
     -Qgmin(id_gen_slack)*ones(T,1)];

%M = 10*sum(Pd);
M = 80;
P_half = 0.5*(2*sum(Pd)-sum(Pgmin(id_gen_nslack))-sum(Pgmax(id_gen_nslack)));
Q_half = 0.5*(2*sum(Qd)-sum(Qgmin(id_gen_nslack))-sum(Qgmax(id_gen_nslack)));
%z_0 = [P_half;Q_half];
z_0 = [0.2*ones(T,1);zeros(T,1)];
% tolerance
tol = 1e-3;

% res = 20;
% z_p = linspace(-0.2,0.8,res);
% z_q = linspace(-2,1.5,res);
% mesh = zeros(res);
% for i = 1:res
%     for j =1:res
%         z_B = [z_p(i);z_q(j)];
%         [mesh(i,j),~] = Boundart_check(b0, beta, gamma, z_B);
%     end
% end


K = 1; K_2 = 0;
while K ~= 0
    [K, z_s] = Boundary_search(b_u, D0, v, B, A_de, M);

    lambda_l = 0;
    lambda_u = 1;
    lambda = 0.5*(lambda_l + lambda_u);
    if K <= 1e-16
        break;
    else
        while ~(lambda_u - lambda_l <= tol && K_2 > 0)
            lambda = 0.5*(lambda_l + lambda_u);
            z_B = lambda*z_s + (1-lambda)*z_0;
            [K_2, h_s] = Boundart_check(b_u, A_de, B, z_B);
            if K_2 == 0
                lambda_l = lambda;
            else 
                lambda_u = lambda;
            end
        end
        D0 = [D0;-(h_s'*B)];
        v = [v;-(h_s'*b_u)];
    end
end

% 设置 x 轴和 y 轴的范围
x = linspace(-0.2, 0.7, 400);
y = linspace(-1.2, 1.2, 400);
[X, Y] = meshgrid(x, y);

% 初始化绘图
fig=figure; box on; hold all; set(fig, 'Position', [100, 100, 850, 650]);

% 设置坐标轴范围和网格线间隔
xlim([-0.2, 0.7]);
ylim([-1.2, 1.2]);
xticks(-0.2:0.1:0.7);
yticks(-1.2:0.5:1.2);
grid on;

% 绘制每个不等式定义的线
for i = 1:size(T,1)
    if T(i,1) == 0
        line(xlim, [v(i)/T(i,2) v(i)/T(i,2)], 'Color', 'r');
    elseif T(i,2) == 0
        line([v(i)/T(i,1) v(i)/T(i,1)], ylim, 'Color', 'r');
    else
        plot(x, (v(i) - T(i,1)*x)/T(i,2), 'r');
    end
end

% 检查每个点是否满足所有不等式
inside = all((A * [X(:), Y(:)]' <= b)');

% 绘制可行区域
fill(X(inside), Y(inside), 'g', 'FaceAlpha', 0.3);

% 设置坐标轴标签和标题
xlabel('x');
ylabel('y');
title('Linear Inequality Constraints');
hold off;