function [residues] = create_ineq_residues(g_u_vvcos,g_l_vvcos,g_u_vvsin,g_l_vvsin,g_u_vv,g_l_vv,gamma_opt,K_plus,K_minus,xi0,xi_gamma,Phi_max,Phi_min,vmax,vmin,idx_pq)
    
    Rmax = K_plus*[g_u_vvcos;g_u_vvsin;g_u_vv]+...
        K_minus*[g_l_vvcos;g_l_vvsin;g_l_vv]+xi0+gamma_opt*xi_gamma-[Phi_max;vmax(idx_pq)];
    Rmin = -K_plus*[g_l_vvcos;g_l_vvsin;g_l_vv]-...
        K_minus*[g_u_vvcos;g_u_vvsin;g_u_vv]-xi0+gamma_opt*xi_gamma+[Phi_min;vmin(idx_pq)];

    residues = vertcat(Rmax,Rmin);
    %2Nbranch+2Npq
end

