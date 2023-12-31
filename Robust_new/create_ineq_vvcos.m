function [vvcos_sin] = create_ineq_vvcos(g_u_vvcos,g_l_vvcos,g_u_vvsin,g_l_vvsin,Dvvu_uu,Dvvu_ll,Dvvu_ul,Dvvu_lu,Dvvl_uu,Dvvl_ll,Dvvl_ul,Dvvl_lu,Dcuu,Dcul,Dclu,Dcll,Dsuu,Dsul,Dslu,Dsll,v_u,v_l,Phi_u,Phi_l, Psi_vvcos0,Psi_vvsin0,Phi0,v0,idx_fr,idx_to,onoff_pq)
    %<=0
    % cos1
    g_gcos11 = -g_u_vvcos+(Psi_vvcos0+cos(Phi0).*Dvvu_uu+v0(idx_fr).*v0(idx_to).*Dcuu+1/4*(Dvvu_uu+Dcuu).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u);
    g_gcos12 = -g_u_vvcos+(Psi_vvcos0+cos(Phi0).*Dvvu_uu+v0(idx_fr).*v0(idx_to).*Dcul+1/4*(Dvvu_uu+Dcul).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l);
    g_gcos13 = -g_u_vvcos+(Psi_vvcos0+cos(Phi0).*Dvvu_ul+v0(idx_fr).*v0(idx_to).*Dcuu+1/4*(Dvvu_ul+Dcuu).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u);
    g_gcos14 = -g_u_vvcos+(Psi_vvcos0+cos(Phi0).*Dvvu_ul+v0(idx_fr).*v0(idx_to).*Dcul+1/4*(Dvvu_ul+Dcul).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l);
    g_gcos15 = -g_u_vvcos+(Psi_vvcos0+cos(Phi0).*Dvvu_lu+v0(idx_fr).*v0(idx_to).*Dcuu+1/4*(Dvvu_lu+Dcuu).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u);
    g_gcos16 = -g_u_vvcos+(Psi_vvcos0+cos(Phi0).*Dvvu_lu+v0(idx_fr).*v0(idx_to).*Dcul+1/4*(Dvvu_lu+Dcul).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l);
    g_gcos17 = -g_u_vvcos+(Psi_vvcos0+cos(Phi0).*Dvvu_ll+v0(idx_fr).*v0(idx_to).*Dcuu+1/4*(Dvvu_ll+Dcuu).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u);
    g_gcos18 = -g_u_vvcos+(Psi_vvcos0+cos(Phi0).*Dvvu_ll+v0(idx_fr).*v0(idx_to).*Dcul+1/4*(Dvvu_ll+Dcul).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l);
    cos1 = vertcat(g_gcos11,g_gcos12,g_gcos13,g_gcos14,g_gcos15,g_gcos16,g_gcos17,g_gcos18);
    % cos2
    g_gcos21 = g_l_vvcos-(Psi_vvcos0+cos(Phi0).*Dvvl_uu+v0(idx_fr).*v0(idx_to).*Dclu-1/4*(Dvvl_uu-Dclu).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u);
    g_gcos22 = g_l_vvcos-(Psi_vvcos0+cos(Phi0).*Dvvl_uu+v0(idx_fr).*v0(idx_to).*Dcll-1/4*(Dvvl_uu-Dcll).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l);
    g_gcos23 = g_l_vvcos-(Psi_vvcos0+cos(Phi0).*Dvvl_ul+v0(idx_fr).*v0(idx_to).*Dclu-1/4*(Dvvl_ul-Dclu).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u);
    g_gcos24 = g_l_vvcos-(Psi_vvcos0+cos(Phi0).*Dvvl_ul+v0(idx_fr).*v0(idx_to).*Dcll-1/4*(Dvvl_ul-Dcll).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l);
    g_gcos25 = g_l_vvcos-(Psi_vvcos0+cos(Phi0).*Dvvl_lu+v0(idx_fr).*v0(idx_to).*Dclu-1/4*(Dvvl_lu-Dclu).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u);
    g_gcos26 = g_l_vvcos-(Psi_vvcos0+cos(Phi0).*Dvvl_lu+v0(idx_fr).*v0(idx_to).*Dcll-1/4*(Dvvl_lu-Dcll).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l);
    g_gcos27 = g_l_vvcos-(Psi_vvcos0+cos(Phi0).*Dvvl_ll+v0(idx_fr).*v0(idx_to).*Dclu-1/4*(Dvvl_ll-Dclu).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u);
    g_gcos28 = g_l_vvcos-(Psi_vvcos0+cos(Phi0).*Dvvl_ll+v0(idx_fr).*v0(idx_to).*Dcll-1/4*(Dvvl_ll-Dcll).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0)+v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l);
    cos2 = vertcat(g_gcos21,g_gcos22,g_gcos23,g_gcos24,g_gcos25,g_gcos26,g_gcos27,g_gcos28);
    %sin1
    g_gsin11 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvu_uu+v0(idx_fr).*v0(idx_to).*Dsuu+1/4*(Dvvu_uu+Dsuu).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin12 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvu_uu+v0(idx_fr).*v0(idx_to).*Dsul+1/4*(Dvvu_uu+Dsul).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    g_gsin13 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvu_ul+v0(idx_fr).*v0(idx_to).*Dsuu+1/4*(Dvvu_ul+Dsuu).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin14 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvu_ul+v0(idx_fr).*v0(idx_to).*Dsul+1/4*(Dvvu_ul+Dsul).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    g_gsin15 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvu_lu+v0(idx_fr).*v0(idx_to).*Dsuu+1/4*(Dvvu_lu+Dsuu).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin16 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvu_lu+v0(idx_fr).*v0(idx_to).*Dsul+1/4*(Dvvu_lu+Dsul).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    g_gsin17 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvu_ll+v0(idx_fr).*v0(idx_to).*Dsuu+1/4*(Dvvu_ll+Dsuu).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin18 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvu_ll+v0(idx_fr).*v0(idx_to).*Dsul+1/4*(Dvvu_ll+Dsul).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    sin1 = vertcat(g_gsin11,g_gsin12,g_gsin13,g_gsin14,g_gsin15,g_gsin16,g_gsin17,g_gsin18);
    %sin2
    g_gsin21 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvl_uu+v0(idx_fr).*v0(idx_to).*Dsuu+1/4*(Dvvl_uu+Dsuu).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin22 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvl_uu+v0(idx_fr).*v0(idx_to).*Dsul+1/4*(Dvvl_uu+Dsul).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    g_gsin23 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvl_ul+v0(idx_fr).*v0(idx_to).*Dsuu+1/4*(Dvvl_ul+Dsuu).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin24 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvl_ul+v0(idx_fr).*v0(idx_to).*Dsul+1/4*(Dvvl_ul+Dsul).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    g_gsin25 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvl_lu+v0(idx_fr).*v0(idx_to).*Dsuu+1/4*(Dvvl_lu+Dsuu).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin26 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvl_lu+v0(idx_fr).*v0(idx_to).*Dsul+1/4*(Dvvl_lu+Dsul).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    g_gsin27 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvl_ll+v0(idx_fr).*v0(idx_to).*Dsuu+1/4*(Dvvl_ll+Dsuu).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin28 = -g_u_vvsin+(Psi_vvsin0+sin(Phi0).*Dvvl_ll+v0(idx_fr).*v0(idx_to).*Dsul+1/4*(Dvvl_ll+Dsul).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    sin2 = vertcat(g_gsin21,g_gsin22,g_gsin23,g_gsin24,g_gsin25,g_gsin26,g_gsin27,g_gsin28);
    %sin3
    g_gsin31 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvl_uu+v0(idx_fr).*v0(idx_to).*Dslu-1/4*(Dvvl_uu-Dslu).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin32 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvl_uu+v0(idx_fr).*v0(idx_to).*Dsll-1/4*(Dvvl_uu-Dsll).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    g_gsin33 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvl_ul+v0(idx_fr).*v0(idx_to).*Dslu-1/4*(Dvvl_ul-Dslu).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin34 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvl_ul+v0(idx_fr).*v0(idx_to).*Dsll-1/4*(Dvvl_ul-Dsll).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    g_gsin35 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvl_lu+v0(idx_fr).*v0(idx_to).*Dslu-1/4*(Dvvl_lu-Dslu).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin36 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvl_lu+v0(idx_fr).*v0(idx_to).*Dsll-1/4*(Dvvl_lu-Dsll).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    g_gsin37 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvl_ll+v0(idx_fr).*v0(idx_to).*Dslu-1/4*(Dvvl_ll-Dslu).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin38 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvl_ll+v0(idx_fr).*v0(idx_to).*Dsll-1/4*(Dvvl_ll-Dsll).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    sin3 = vertcat(g_gsin31,g_gsin32,g_gsin33,g_gsin34,g_gsin35,g_gsin36,g_gsin37,g_gsin38);
    %sin4
    g_gsin41 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvu_uu+v0(idx_fr).*v0(idx_to).*Dslu-1/4*(Dvvu_uu-Dslu).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin42 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvu_uu+v0(idx_fr).*v0(idx_to).*Dsll-1/4*(Dvvu_uu-Dsll).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    g_gsin43 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvu_ul+v0(idx_fr).*v0(idx_to).*Dslu-1/4*(Dvvu_ul-Dslu).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin44 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvu_ul+v0(idx_fr).*v0(idx_to).*Dsll-1/4*(Dvvu_ul-Dsll).^2-onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    g_gsin45 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvu_lu+v0(idx_fr).*v0(idx_to).*Dslu-1/4*(Dvvu_lu-Dslu).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin46 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvu_lu+v0(idx_fr).*v0(idx_to).*Dsll-1/4*(Dvvu_lu-Dsll).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    g_gsin47 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvu_ll+v0(idx_fr).*v0(idx_to).*Dslu-1/4*(Dvvu_ll-Dslu).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u);
    g_gsin48 = g_l_vvsin-(Psi_vvsin0+sin(Phi0).*Dvvu_ll+v0(idx_fr).*v0(idx_to).*Dsll-1/4*(Dvvu_ll-Dsll).^2-onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0)-v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0)-v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l);
    sin4 = vertcat(g_gsin41,g_gsin42,g_gsin43,g_gsin44,g_gsin45,g_gsin46,g_gsin47,g_gsin48);

    vvcos_sin = vertcat(cos1,cos2,sin1,sin2,sin3,sin4);
    %48*Nbranch
end

