%% iv
%
%  Function to create IV curve data based on one-diode model parameters of
%  a PV module.
%
%  Params:
%      n     - ideality factor
%      i_o   - diodes saturation current
%      i_pv  - photocurrent
%      rs    - series resistance
%      rsh   - shunt resistance
%      N     - number of PV cells in the PV module
%      range - range can be a vector or the max value (in the latter case,
%              range will be 0-max.
%
%  Returns:
%      v_out - voltage values of the iv curve
%      i_out - current values of the iv curve
%

function [v_out,i_out] = iv(n, i_o, i_pv, rs, rsh, N, range)
    k = 1.38064852E-23;
    q = 1.60217662E-19;
    T = 298.15;
    vt = k*T/q;

    if(isscalar(range))
        v_range = 0:range/1000:range;
    else
        v_range = range;
    end
    
    f = @(i) i_pv - i_o * (exp((v_range + i*rs)/(n*N*vt)) - 1) - (v_range + i*rs)/rsh - i;
    
    i_num = fsolve(f, i_pv*ones(size(v_range)));

    i_out = i_num(i_num >= 0);
    v_out = v_range(1:length(i_out));
end