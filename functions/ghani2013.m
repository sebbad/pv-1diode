%% Ghani 2013
%  Based on "Numerical calculation of series and shunt
%  resistance and diode quality factor of a photovoltaic cell using the
%  Lambert W-function" Ghani et al. 2013
%  
%  params:
%    pv_data   - remarkable points of the pv device
%    rs_guess  - initial guess for rs value
%    rsh_guess - initial guess for rsh value
%    tol_r     - stop condition for r loop
%    tol_n     - stop condition for n loop
%
%  returns:
%    intern    - structure with internal parameters of interest
%    n         - estimated ideality factor
%    i_0       - estimated value for i_0
%    i_pv      - estimated value for i_pv
%    rs        - estimated value for rs
%    rsh       - estimated value for rsh
%

function [intern, n, i_0, i_pv, rs, rsh] = ghani2013(pv_data, rs_guess, rsh_guess, tol_r, tol_n)

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


%% Initial assumptions
i_pv = i_sc;
i_0 = @(n) i_sc/(exp(v_oc/(n*N*vt)) - 1);


%% Set up the equations
%  The method looks to find solutions that make z and w zero. To do so, a
%  number of partial derivitives are needed.

syms rs rsh n;

theta1 = (rs*i_0(n)*rsh)/(n*N*vt*(rs+rsh)) * exp((rsh*(rs*i_pv+rs*i_0(n)+v_mp))/(n*N*vt*(rs+rsh)));

z = (rsh*(i_0(n)+i_pv))/(rs+rsh) - v_mp/(rs+rsh) - (n*N*vt/rs)*lambertw(theta1) - i_mp;

theta2 = (i_0(n)*rsh)/(n*N*vt) * exp((rsh*(-i_mp+i_pv+i_0(n)))/(n*N*vt));

w = rsh*i_pv + rsh*i_0(n) - 2*i_mp*rs - 2*i_mp*rsh - n*N*vt*lambertw(theta2) ...
    + i_mp*rsh*(lambertw(theta2)/(1+lambertw(theta2)));

dth1_drs_sym = diff(theta1, rs);
dth1_drs = matlabFunction(dth1_drs_sym);

dz_drs_sym = diff(z, rs);
dz_drs = matlabFunction(dz_drs_sym);

dth1_drsh_sym = diff(theta1, rsh); 
dth1_drsh = matlabFunction(dth1_drsh_sym);

dz_drsh_sym = diff(z, rsh);
dz_drsh = matlabFunction(dz_drsh_sym);

dw_drs_sym = diff(w, rs);
dw_drs = matlabFunction(dw_drs_sym);

dth2_drsh_sym = diff(theta2, rsh);
dth2_drsh = matlabFunction(dth2_drsh_sym);

dw_drsh_sym = diff(w, rsh);
dw_drsh = matlabFunction(dw_drsh_sym);

z = matlabFunction(z);
w = matlabFunction(w);


%% Start with the iteration process

%  Just some starting values for the FF error and n
n = 1 - 0.001;

n_array = [];
rs_e_array = [];
rsh_e_array = [];
eps_e_array = [];

%  Iterate over different n's as long as the FF error is not good enough
for iter=1:20000  % stop condition if tolerance is not reached
    
    n = n + 0.001;
    n_array = [n_array n];

    clear rs rsh;

    % Starting values for the Newton-Raphson method
    rs = rs_guess;
    rsh = rsh_guess;
    eps = 1;
    max_iter = 1000;
    
    eps_i_array = [];
    rs_i_array = [];
    rsh_i_array = [];    

    while eps > tol_r && max_iter ~= 0

        % Correct rs and rsh
        rs = rs - (z(n,rs,rsh)*dw_drsh(n,rsh) - w(n,rs,rsh)*dz_drsh(n,rs,rsh))/(dz_drs(n,rs,rsh)*dw_drsh(n,rsh) - dz_drsh(n,rs,rsh)*dw_drs());
        rs_i_array = [rs_i_array rs];
        
        rsh = rsh - (w(n,rs,rsh)*dz_drs(n,rs,rsh) - z(n,rs,rsh)*dw_drs())/(dz_drs(n,rs,rsh)*dw_drsh(n,rsh) - dz_drsh(n,rs,rsh)*dw_drs());
        rsh_i_array = [rsh_i_array rsh];

        % Calculate the error
        eps = max([abs(z(n,rs,rsh)), abs(w(n,rs,rsh))]);
        eps_i_array = [eps_i_array eps];
        
        max_iter = max_iter - 1;

    end
    
    intern.tol_r = eps;
    
    if(max_iter == 0)
        intern.tol_r = min(eps_i_array);
        [~, idx] = min(eps_i_array);
        rs = rs_i_array(idx);
        rsh = rsh_i_array(idx);
    end
    
    
    rs_e_array = [rs_e_array rs];
    rsh_e_array = [rsh_e_array rsh];

    % Once we found good estimates for rs and rsh, we want to calculate the
    % FF error. First we implement the iv curve function (Eq.5)
    iv = @(v) -v/(rs+rsh) - (lambertw((rs*i_0(n)*rsh)/(n*N*vt*(rs+rsh)) * exp(rsh*(rs*i_pv+rs*i_0(n)+v)/(n*N*vt*(rs+rsh)))))*n*N*vt/rs + (rsh*(i_0(n)+i_pv))/(rs+rsh);

    % Calculate v,i and p for a fine v range
    v = linspace(0,v_oc,1000);
    i = iv(v);
    p = v.*i;
    
    % Find the VOC and ISC of our model
    [~,idx] = min(abs(i));
    v_oc_m = v(idx);
    i_sc_m = iv(0);    

    % Calculate the FF of our model and compare it to the datasheet FF
    ff_m = max(p)/(v_oc_m*i_sc_m);
    ff_a = (v_mp*i_mp)/(v_oc*i_sc);

    eps_ff = abs(ff_a - ff_m);
    eps_e_array = [eps_e_array eps_ff];
      
    if(eps_ff < tol_n)
        break;
    end
end

intern.tol_n = eps_ff;

if(eps_ff > tol_n)
    intern.tol_n = min(eps_e_array);
    [~, idx] = min(eps_e_array);
    n = n_array(idx);
    rs = rs_e_array(idx);
    rsh = rsh_e_array(idx);
end

i_0 = i_0(n);

end