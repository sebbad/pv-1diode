%% Villalva 2009
%  TBased on "Comprehensive approach to modeling and
%  simulation of photovoltaic arrays" Villalva et al. 2009.
%  
%  params:
%    pv_data   - remarkable points of the pv device
%    rs_range  - max. value for rs iteration
%    n         - prescribed n value
%    tol       - stop condition for early iteration break
%
%  returns:
%    intern    - structure with internal parameters of interest
%    n         - estimated ideality factor
%    i_0       - estimated value for i_0
%    i_pv      - estimated value for i_pv
%    rs        - estimated value for rs
%    rsh       - estimated value for rsh
%

function [intern, n, i_0, i_pv, rs, rsh] = villalva2009(pv_data, rs_range, n, tol)

%% Set parameters
v_oc = pv_data(1);
i_sc = pv_data(2);
v_mp = pv_data(3);
i_mp = pv_data(4);
N = pv_data(5);

k = 1.38064852e-23;
q = 1.60217662e-19;
T = 298.15;

vt = k*T/q;

intern = struct;

%% Calculate i0 and initial rsh
i_0 = i_sc/(exp(v_oc/(n*N*vt)) - 1);
rsh = v_mp/(i_sc-i_mp) - (v_oc-v_mp)/i_mp;

ipv_array = [];
rsh_array = [];
e_array = [];
e = 1;


%% Iterate through rs
for r = 1:length(rs_range)

    rs = rs_range(r);

    % Calculate i_pv based on Eq.10
    i_pv = (rsh+rs)/rsh * i_sc;
    ipv_array = [ipv_array i_pv];

    % Calculate new rsh based on Eq.9
    rsh = v_mp*(v_mp+i_mp*rs) / (v_mp*i_pv - v_mp*i_0*exp((v_mp+i_mp*rs)/(n*N*vt)) + v_mp*i_0 - v_mp*i_mp);
    rsh_array = [rsh_array rsh];

    % Predict the iv-curve with the estimated values and calculate p_mp
    try
    [v,i] = iv(n,i_0,i_pv,rs,rsh,N,v_oc);
    p = v.*i;

    % Calculate error of the p_mp estimation
    e = abs(max(p) - v_mp*i_mp);
    e_array = [e_array e];
    catch
    end
    
    % Break if the error is small enough
    if (e < tol)
        intern.tol = e;
        break;
    end
end

%  Check if we found a set with e < tol
if (r == length(rs_range))
    disp('Fail');
    
    intern.tol = min(e_array);
    [~, idx] = min(e_array);
    
    rs = rs_range(idx);
    rsh = rsh_array(idx);
    i_pv = ipv_array(idx);    
end