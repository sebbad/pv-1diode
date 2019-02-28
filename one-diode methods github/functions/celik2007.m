%% Celik 2007
%  Based on "Modelling and experimental verification of the
%  operating current of mono-crystalline photovoltaic modules using four-
%  and five-parameter models" Celik et al. 2007
%  
%  params:
%    pv_data   - remarkable points of the pv device
%    pv_curve  - experimental curve to extract rs0 and rsh0 from
%    r_avg     - number of samples to take into account for rs0 and rsh0
%    g_ratio   - ratio between the scaled and reference irradiance
%    condition
%
%  returns:
%    intern    - structure with internal parameters of interest
%    n         - estimated ideality factor
%    i_0       - estimated value for i_0
%    i_pv      - estimated value for i_pv
%    rs        - estimated value for rs
%    rsh       - estimated value for rsh
%

function [intern, n, i_0, i_pv, rs, rsh] = celik2007(pv_data, pv_curve, r_avg, g_ratio)

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


%% Calculate rs0 and rsh0 from the slopes of the provided iv_curve
rs0 = rs_slope(pv_curve, r_avg);
rsh0 = rsh_slope(pv_curve, r_avg);

intern.rs0 = rs0;
intern.rsh0 = rsh0;


%% Analytically estimate the five parameters of the one-diode model

rsh = rsh0;

%  Calculate n for reference conditions first
n = (v_mp+i_mp*rs0-v_oc)/(vt*(log(i_sc-v_mp/rsh-i_mp)-log(i_sc-v_oc/rsh)+(i_mp/(i_sc-v_oc/rsh)))) / N;

%  Scale i_sc and v_oc for irradiance condition
i_sc = i_sc * g_ratio;
v_oc = v_oc + n*N*vt * log(g_ratio);

i_0 = (i_sc-v_oc/rsh)*exp(-v_oc/(n*N*vt));

rs = rs0 - (n*N*vt/i_0 * exp(-v_oc/(n*N*vt)));

i_pv = i_sc*(1 + rs/rsh) + i_0*(exp((i_sc*rs)/(n*N*vt))-1);
