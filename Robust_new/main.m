mpc = loadcase('case9');

opts = mpoption;
opts.opf.violation   = 1e-12;
opts.mips.costtol    = 1e-12;
opts.mips.gradtol    = 1e-12;
opts.mips.comptol    = 1e-12;
opts.opf.ignore_angle_lim = true;
opts.out.all = 0; 

mpc = runopf(mpc,opts);

id_load = find((mpc.bus(:,PD).^2+mpc.bus(:,QD).^2) ~= 0);
pl0 = mpc.bus(id_load,3)/mpc.baseMVA;
ql0 = mpc.bus(id_load,4)/mpc.baseMVA;
Sigmap0 = diag(pl0.^2);
Sigmaq0 = diag(ql0.^2);
gamma_0 = 0;
mpc.uncertainty.Sigmap = Sigmap0; 
mpc.uncertainty.Sigmaq = Sigmaq0;
mpc.uncertainty.gamma0 = gamma_0; 

mpc.gen(:,end+1) = [0 0.5 0.5]; % 分配slack bus 不平衡量给PV节点的比例
mpc.gen(:,10) = -mpc.gen(:,9)*0.2; % PV节点PMIN
mpc.bus(:,9) = deg2rad(mpc.bus(:,9)); % 角度转弧度

idx_slack = find(mpc.bus(:,2)==3);
plot_rng = [-2 3 -1 1]; %画图范围
resolution = 100; % 画图分辨率
contour_plot(mpc,idx_slack,plot_rng,resolution);