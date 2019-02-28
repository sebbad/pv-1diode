%% Hussein 2017
%  Based on "A simple approach to extract the unknown
%  parameters of PV modules" Hussein 2017. 
%  
%  params:
%    pv_data   - remarkable points of the pv device
%    rs_guess  - initial guess for rs value
%    n_start   - start value for ideality factor
%    eval_data - experimental data to compare the model error with
%    tol       - stop condition for rs search
%
%  returns:
%    intern    - structure with internal parameters of interest
%    n         - estimated ideality factor
%    i_0       - estimated value for i_0
%    i_pv      - estimated value for i_pv
%    rs        - estimated value for rs
%    rsh       - estimated value for rsh
%

function [intern, n, i_0, i_pv, rs, rsh] = hussein2017(pv_data, rs_guess, n_start, eval_data, tol)

addpath('../');

%% Set initial constants
v_oc = pv_data(1);
i_sc = pv_data(2);
v_mp = pv_data(3);
i_mp = pv_data(4);
N = pv_data(5);

n = n_start;            % intial value for n
rs = rs_guess;          % initial value for rs

iter = 10000;           % maximum number of iterations for n
rms_err = 1e6;          % start with a large rms error
rms_n = 1e6;

T = 298.15;             % temperature (in K)
k = 1.38064852E-23;     % boltzmann constant
q = 1.60217662E-19;     % electron charge

vt = k*T/q;             % thermal voltage

intern = struct;


%% Numerically evaluate r for a given n

%  Equation 17 from the article
f1 = @(n,r) v_mp*(exp((v_mp+i_mp*r)/(n*N*vt)))*(i_sc*v_oc-i_sc*v_mp-i_mp*v_oc) ...
           - i_sc*v_mp*n*N*vt*(exp(v_oc/(n*N*vt))-1 - exp((v_mp+i_mp*r)/(n*N*vt))-1) ...
           + i_mp*v_oc*n*N*vt*(exp((i_sc*r)/(n*N*vt))-1 - exp((v_mp+i_mp*r)/(n*N*vt))-1) ...
           + 2*i_mp*v_mp*n*N*vt*(exp(v_oc/(n*N*vt))-1 - exp((i_sc*r)/(n*N*vt))-1);

%  Equation 18 from the article
f2 = @(n,r) (v_oc*v_mp*(exp((v_mp+i_mp*r)/(n*N*vt))))/(n*N*vt)*(i_sc-i_mp) ...
            + i_sc*(exp((i_sc*r)/(n*N*vt)))*(v_oc-2*v_mp) ...
            + (exp((v_mp+i_mp*r)/(n*N*vt)))*(i_sc*v_mp - i_mp*v_oc - (i_sc*(v_mp)^2)/(n*N*vt));

%  Equation 14 from the article
eq14 = @(n,r) (v_oc*(exp((v_mp+i_mp*r)/(n*N*vt))-1 - exp((i_sc*r)/(n*N*vt))-1) - v_mp*(exp(v_oc/(n*N*vt))-1 - exp((i_sc*r)/(n*N*vt))-1)) ...
              / (i_mp*(exp(v_oc/(n*N*vt))-1 - exp((i_sc*r)/(n*N*vt))-1) - i_sc*(exp(v_oc/(n*N*vt))-1 - exp((v_mp+i_mp*r)/(n*N*vt))-1)) - r;

%  Equation 12 from the article
eq12 = @(n,r,rsh) (i_sc*(1 + r/rsh) * (exp(v_oc/(n*N*vt))-1) - v_oc/rsh) / (exp(v_oc/(n*N*vt))-1 - exp((i_sc*r)/(n*N*vt))-1);

%  Equation 5 from the article
eq5 = @(n,ipv,rsh) ipv/(exp(v_oc/(n*N*vt))-1) - v_oc/((exp(v_oc/(n*N*vt))-1)*rsh);


%  As long as we haven't tried the maximum number of iterations (or found
%  the min)
while iter

    x = f1(n,rs);
    y = f2(n,rs);

    % Check if f1/f2 is smaller than the tolerance
    while (abs(x/y) > tol)
        rs = rs - x/y;    % here there is an error in the article: the article states '+', but it should be '-'

        x = f1(n,rs);
        y = f2(n,rs);
    end
    
    intern.tol = abs(x/y);

    % When we found rs, we can calculate the remaining parameters
    rsh = eq14(n,rs);
    i_pv = eq12(n,rs,rsh);
    i_0 = eq5(n,i_pv,rsh);

    try
        % Estimate the rms error based on the obtained parameters
        [~, i] = iv(n, i_0, i_pv, rs, rsh, N, eval_data(:,1));
        err = i - eval_data(:,2);
        rms_n = rms(err);
    catch
    end

    % If the obtained rms error is smaller than the previous iteration, we
    % continue ...
    if(rms_n <= rms_err)
            rms_err = rms_n;
            n = n+1e-3; % change n, and try again
            iter = iter-1;
    else
        break; % otherwise we stop with the current iteration
    end
    
    intern.rms = rms_err;

end
end
