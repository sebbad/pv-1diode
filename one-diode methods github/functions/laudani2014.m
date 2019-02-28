%% Laudani2014
%  Based on "Identification of the on-diode model for
%  photovoltaic modules from datasheet values", Laudani et al. 2014.
%  
%  params:
%    pv_data   - remarkable points of the pv device
%    rs_guess  - initial guess for rs value
%    n_range   - max. n value to estimate rs for
%
%  returns:
%    intern    - structure with internal parameters of interest
%    n         - estimated ideality factor
%    i_0       - estimated value for i_0
%    i_pv      - estimated value for i_pv
%    rs        - estimated value for rs
%    rsh       - estimated value for rsh
%

function [intern, n, i_0, i_pv, rs, rsh] = laudani2014(pv_data, rs_guess, n_range)

addpath('../');

%% Set initial constants
v_oc = pv_data(1);
i_sc = pv_data(2);
v_mp = pv_data(3);
i_mp = pv_data(4);
N = pv_data(5);

k = 1.38064852E-23;             % boltzmann constant
q = 1.60217662E-19;             % electron charge
T = 298.15;                     % temperature (in K)
vt = k*T/q;                     % thermal voltage

intern = struct;

%% Estimation of the max. ideality factor (n)
%  Defined symbols to simplify equations
P1 = v_mp*i_mp;                
P2 = (v_oc-v_mp)*i_mp;
P3 = (v_oc-v_mp)*(i_sc-i_mp);
P4 = v_mp*(i_sc-i_mp);

%  Eq. 29: rs = f(n)
f1 = @(n,r) (P2-P1)*exp((i_sc*r)/(n*N*vt)) + (P1-P4)*exp((v_oc)/(n*N*vt)) + ((P1-P3)*(i_mp*r-v_mp)/(n*N*vt) + (P4-P2))*exp((v_mp+i_mp*r)/(n*N*vt));
%  Eq. 26: rs_max = f(n) (Based on an earlier work of the same author)
rs_max = @(x) v_mp/i_mp + ((N*x*vt)/i_mp)*(1+lambertw(-1,-exp((v_oc-x*N*vt-2*v_mp)/(N*x*vt))));

rs_vec = [];
rs_max_vec = [];

for i=1:length(n_range)
    f = @(r) f1(n_range(i),r);                       % reduce f1 to f(r)
    rs_temp = fsolve(f,rs_guess);                    % numerically solve f
    rs_vec = [rs_vec rs_temp];                       % add the new rs to the array
    rs_max_vec = [rs_max_vec rs_max(n_range(i))];    % calculate the max. rs for n
end

intern.rs_vec = rs_vec;
intern.rs_max_vec = rs_max_vec;

%  According to the article, n_max is the intersection between the rs and
%  rs_max curves.
diff = abs(rs_vec-rs_max_vec);           % calculate the difference for each value pair
[~, min_diff_idx] = min(diff);           % find the closest value to the intersection

n_max = n_range(min_diff_idx);

intern.nmax = n_max;


%% Estimation of the parameters

%  Follow the heuristic rule n = 0.9 * n_max
n = 0.9 * n_max;

%  Determine rs the same way as before (Eq. 29)
f = @(r) f1(n,r);
rs = fsolve(f,rs_guess);

denom = ((v_mp+rs*i_mp-v_oc)*exp((i_sc*rs)/(n*N*vt)) + (v_oc-rs*i_sc)*exp((v_mp+i_mp*rs)/(n*N*vt)) + (rs*i_sc-rs*i_mp-v_mp)*exp((v_oc)/(n*N*vt)));

%  Determine rsh (Eq. 22)
gsh = (exp((v_oc)/(n*N*vt))*(i_mp-i_sc) + exp((v_mp+i_mp*rs)/(n*N*vt))*i_sc - exp((i_sc*rs)/(n*N*vt))*i_mp)/denom;
rsh = 1/gsh;

%  Determine i_0 (Eq. 23)
i_0 = (v_oc*(i_sc-i_mp) - v_mp*i_sc)/denom;

%  Determine i_pv (Eq. 24)
i_pv = (i_sc*v_oc*(exp((v_mp+i_mp*rs)/(n*N*vt)) - 1) + i_sc*v_mp*(1 - exp((v_oc)/(n*N*vt))) + i_mp*v_oc*(1 - exp((i_sc*rs)/(n*N*vt))))/denom;

end

