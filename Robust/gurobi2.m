import casadi.*

% 定义决策变量
x = SX.sym('x', 2);

% 定义目标函数的二次和线性部分
H = [2, 1;    % 二次目标函数系数矩阵
     1, 2];
g = [-1; -1]; % 线性目标函数系数

% 约束的系数
A = [1, 1; 
     1, -1; 
     1, 0;
     0, 1];
lba = [1; -inf; 0; 0];
uba = [2; 1; 1; inf];

% 创建QP问题结构
qp = struct('h', H, 'g', g, 'a', A);

% 设置Gurobi选项
opts = struct;
opts.print_time = true;
opts.verbose = true;
opts.verbose_init = true;

% 创建求解器
solver = qpsol('solver', 'gurobi', qp, opts);

% 调用求解器
sol = solver('lba', lba, 'uba', uba);

% 输出结果
disp('Solution:');
disp(full(sol.x));
