%% Cubas 2014
%  Based on "On the analytical approach for modeling
%  photovoltaic systems behavior", Cubas et al. 2014
%  
%  params:
%    pv_data   - remarkable points of the pv device
%    rs_guess  - initial guess for rs value
%
%  returns:
%    intern    - structure with internal parameters of interest
%    n         - estimated ideality factor
%    i_0       - estimated value for i_0
%    i_pv      - estimated value for i_pv
%    rs        - estimated value for rs
%    rsh       - estimated value for rsh
%

function [intern, n, i_0, i_pv, rs, rsh] = cubas2014(pv_data, rs_guess)

addpath('../');

%% Set initial constants
n = 1.3;

v_oc = pv_data(1);
i_sc = pv_data(2);
v_mp = pv_data(3);
i_mp = pv_data(4);
N = pv_data(5);

T = 298.15;
k = 1.38064852E-23;
q = 1.60217662E-19;

Vt = k*T/q;

intern = struct;


%% Rs as a function of a (the ideality factor)
%  Equation 16 in the article.
%  Needs to be solved numerically.
eq16 = @(r_s) (n*N*Vt*v_mp*(2*i_mp-i_sc))/((v_mp*i_sc+v_oc*(i_mp-i_sc))*(v_mp-i_mp*r_s)-n*N*Vt*(v_mp*i_sc-v_oc*i_mp))-exp((v_mp+i_mp*r_s-v_oc)/(n*N*Vt));

%  Initial guess of rs.
rs = fsolve(eq16, rs_guess);


%% Once we have rs, the other parameters (rsh, io, ipv) can be determined explicitly.
rsh = ((v_mp-i_mp*rs)*(v_mp-rs*(i_sc-i_mp)-n*N*Vt))/((v_mp-i_mp*rs)*(i_sc-i_mp)-n*N*Vt*i_mp);    % Eq. 17
i_0 = ((rsh+rs)*i_sc-v_oc)/(rsh*exp(v_oc/(n*N*Vt)));                                             % Eq. 9
i_pv = (rsh+rs)/rsh * i_sc;                                                                      % Eq. 7

end