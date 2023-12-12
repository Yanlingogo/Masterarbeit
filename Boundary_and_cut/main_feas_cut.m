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
%% parameters
baseMVA     = mpc.baseMVA;
baseKV      = mpc.bus(1,BASE_KV);
Umax        = mpc.bus(:,VMAX).^2;                                            
Umin        = mpc.bus(:,VMIN).^2;
Pgmin       = mpc.gen(:,PMIN)/baseMVA;  
Qgmin       = mpc.gen(:,QMIN)/baseMVA; 
Pgmax       = mpc.gen(:,PMAX)/baseMVA; 
Qgmax       = mpc.gen(:,QMAX)/baseMVA; 
Pd          = mpc.bus(:,PD)/baseMVA;
Qd          = mpc.bus(:,QD)/baseMVA;

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
Cg_nslack      = sparse(id_gen_nslack,1:Ngen_nslack,ones(Ngen_nslack,1),Nbus,Ngen_nslack);
Cg_slack       = sparse(id_gen_slack, 1, ones(1),Nbus,1);
%% Problem formulation
beta = [eye(Nbus) zeros(Nbus,2*Ngen_nslack+2*Nbranch);
        -eye(Nbus) zeros(Nbus,2*Ngen_nslack+2*Nbranch);
        zeros(Ngen_nslack,Nbus) eye(Ngen_nslack) zeros(Ngen_nslack,Ngen_nslack+2*Nbranch);
        zeros(Ngen_nslack,Nbus) -eye(Ngen_nslack) zeros(Ngen_nslack,Ngen_nslack+2*Nbranch);
        zeros(Ngen_nslack,Nbus+Ngen_nslack) eye(Ngen_nslack) zeros(Ngen_nslack,2*Nbranch);
        zeros(Ngen_nslack,Nbus+Ngen_nslack) -eye(Ngen_nslack) zeros(Ngen_nslack,2*Nbranch);
        Cf-Ct zeros(Nbranch,2*Ngen_nslack) diag(-2*branch_r) diag(-2*branch_x);
        -(Cf-Ct) zeros(Nbranch,2*Ngen_nslack) diag(2*branch_r) diag(2*branch_x);% branch flow equations
        zeros(Nbus) Cg_nslack zeros(Nbus, Ngen_nslack) -Cf'+Ct' zeros(Nbus, Nbranch);
        zeros(Nbus) -Cg_nslack zeros(Nbus, Ngen_nslack) Cf'-Ct' zeros(Nbus, Nbranch);
        zeros(Nbus) zeros(Nbus, Ngen_nslack) Cg_nslack zeros(Nbus, Nbranch) -Cf'+Ct';
        zeros(Nbus) zeros(Nbus, Ngen_nslack) -Cg_nslack zeros(Nbus, Nbranch) Cf'-Ct';];

gamma = [zeros(2*Nbus+4*Ngen_nslack+2*Nbranch,2);
         Cg_slack zeros(size(Cg_slack));
         -Cg_slack zeros(size(Cg_slack));
         zeros(size(Cg_slack)) Cg_slack;
         zeros(size(Cg_slack)) -Cg_slack];

b0 = [Umax;-Umin;Pgmax(id_gen_nslack);-Pgmin(id_gen_nslack);...
      Qgmax(id_gen_nslack);-Qgmin(id_gen_nslack);zeros(2*Nbranch,1);...
      Pd;-Pd;Qd;-Qd];
% initilization 
T = [ 1 0;
     -1 0;
      0 1;
      0 -1;];
v = [sum(Pd)-sum(Pgmin(id_gen_nslack));
     -(sum(Pd)-sum(Pgmax(id_gen_nslack)));
     sum(Qd)-sum(Qgmin(id_gen_nslack));
     -(sum(Qd)-sum(Qgmax(id_gen_nslack)))];
%M = 10*sum(Pd);
M = 10;
P_half = 0.5*(2*sum(Pd)-sum(Pgmin(id_gen_nslack))-sum(Pgmax(id_gen_nslack)));
Q_half = 0.5*(2*sum(Qd)-sum(Qgmin(id_gen_nslack))-sum(Qgmax(id_gen_nslack)));
%z_0 = [P_half;Q_half];
z_0 = [0.2;0];
% tolerance
tol = 1e-3;

K = 1;
while K >= 0
    [K, z_s] = Boundary_search(b0, T, v, gamma, beta, M);
    if K <= 0
        break;
    else
        lambda_l = 0;
        lambda_u = 1;
        lambda = 0.5*(lambda_l + lambda_u);
        while lambda_u - lambda_l >= tol
            lambda = 0.5*(lambda_l + lambda_u);
            z_B = lambda*z_s + (1-lambda)*z_0;
            [K_2, h_s] = Boundart_check(b0, beta, gamma, z_B);
            if K_2 <= 1e-16
                lambda_l = lambda;
            else 
                lambda_u = lambda;
            end
        end
        T = [T;-(h_s'*gamma)];
        v = [v;-(h_s'*b0)];
    end
end



% 使用 polytope 函数（需要 Multi-Parametric Toolbox）
P = Polyhedron('A', T, 'b', v);

% 计算顶点
V = P.V;

% 绘制多边形
fill(V(:,1), V(:,2), 'b');
alpha(0.3); % 使填充半透明
xlabel('P');
ylabel('Q');
grid on;