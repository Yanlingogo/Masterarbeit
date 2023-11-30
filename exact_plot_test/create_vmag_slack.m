function [v_slack] = create_vmag_slack(v,slack_v,vmax,vmin)
    ineq1 =  v + slack_v - vmax;
    ineq2 = -v + slack_v + vmin;

    v_slack = vertcat(ineq1,ineq2);

    %2Nbus
end

