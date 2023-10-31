function [mpc] = runpf_cvxr(mpc)
    %% Index setting
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

    Nbus = size(mpc.bus,1);
    Ngen = size(mpc.gen,1);
    Npq = sum(mpc.bus(:,BUS_TYPE)==1);
    bus = mpc.bus;
    gen = mpc.gen;

    bus_type = bus(:,BUS_TYPE);
    slack_bus = bus_type == 3;
    idx_nslack = find(bus_type ~= 3);
    id_gen = gen(:,GEN_BUS);
    idx_pq = find(bus(:,BUS_TYPE) == 1);

    gencost = mpc.gencost(:,COST:end);

    Cg = sparse(id_gen, (1:Ngen), ones(1, Ngen), Nbus, Ngen);
    Cl = sparse(idx_pq, (1:Npq), ones(1, Npq), Nbus, Npq);

    gen_status = gen(:,GEN_STATUS);% no load status
    pg0 = gen(:,GEN_STATUS).*gen(:,PG); 
    ppq0 = bus(idx_pq,PD);
    qpq0 = bus(idx_pq,QD);
    pinj0 = Cg*pg0 -Cl*ppq0; qinj0 = -Cl*qpq0;
    Sinj = (pinj0 + 1i * qinj0)/mpc.baseMVA; % in p.u.
    
    [Y,Yf,Yt] = makeYbus(mpc);

    vm = bus(:,VM);
    vm(id_gen) = gen(:,VG);
    va = deg2rad(bus(:,VA)); % Angle value to radian system
    va(slack_bus) = 0;
    delta = 0;
    alpha = gen(:,end);

    max_iter = 30;
    nf = zeros(max_iter,1); ndx = zeros(max_iter,1);

    for iter = 1:max_iter
        v_cpx = vm.*cos(va)+1i*(vm.*sin(va));
        S_bal = v_cpx.*conj(Y*v_cpx) - Sinj-Cg*alpha*delta;
        f=[real(S_bal);imag(S_bal(idx_pq))];
        J1 = real(diag(v_cpx) * conj(Y .* (ones(Nbus, 1) * (1i * v_cpx)')) + diag(1i * v_cpx) * diag(conj(Y * v_cpx)));
        J2 = real(diag(v_cpx) * conj(Y .* (ones(Nbus, 1) * (v_cpx./vm)')) + diag(v_cpx ./ vm) * diag(conj(Y * v_cpx)));
        J3 = imag(diag(v_cpx) * conj(Y .* (ones(Nbus, 1) * (1i * v_cpx)')) + diag(1i * v_cpx) * diag(conj(Y * v_cpx)));
        J4 = imag(diag(v_cpx) * conj(Y .* (ones(Nbus, 1) * (v_cpx./vm)')) + diag(v_cpx ./ vm) * diag(conj(Y * v_cpx)));
        J  = [J1(:,idx_nslack) J2(:,idx_pq) -Cg*alpha; J3(idx_pq,idx_nslack) J4(idx_pq,idx_pq) zeros(length(idx_pq), 1)];
        dx = -J\f;

        va(idx_nslack) = va(idx_nslack)+dx(1:Nbus-1);
        vm(idx_pq) = vm(idx_pq)+dx(Nbus:end-1);
        delta = delta+dx(end);

        nf(iter) = norm(f); ndx(iter) = norm(dx);
        if isnan(nf(iter)) || ((nf(iter)<1e-3)&&(ndx(iter)<1e-3))
            break; 
        end 
    end

    mpc.bus(:,VM) = vm;
    mpc.bus(:,VA) = va;
    mpc.delta = delta;


    v_cpx=vm.*cos(va)+1i*vm.*sin(va);
    S_inj=v_cpx.*conj(Y*v_cpx)+1i*Cl*(qpq0/mpc.baseMVA);
    qg_inj = Cg' * diag(1 ./ sum(Cg, 2)) * imag(S_inj)*mpc.baseMVA;
    mpc.gen(:,QG) = qg_inj(id_gen);
    mpc.cost = gencost(:,1)'*(pg0+alpha.*delta).^2 + gencost(:,2)'*(pg0+alpha.*delta) + sum(gencost(:,3));

    if nf(end)>1e-8
        mpc.delta = NaN;
    end
end

