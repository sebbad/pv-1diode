%% Brano 2010
%  Based on "An improved five-parameter model for photovoltaic modules" Brano et al. 2010.
%  
%  params:
%    pv_data   - remarkable points of the pv device
%    pv_curve  - experimental curve to extract rs0 and rsh0 from
%    r_avg     - number of samples to take into account for rs0 and rsh0
%    n_guess   - initial guess for ideality factor
%
%  returns:
%    intern    - structure with internal parameters of interest
%    n         - estimated ideality factor
%    i_0       - estimated value for i_0
%    i_pv      - estimated value for i_pv
%    rs        - estimated value for rs
%    rsh       - estimated value for rsh
%

function [intern, n, i_0, i_pv, rs, rsh] = brano2010(pv_data, pv_curve, r_avg, n_guess)

%% Set parameters
v_oc = pv_data(1);
i_sc = pv_data(2);
v_mp = pv_data(3);
i_mp = pv_data(4);

T = 298.15;

intern = struct;


%% Calculate rs0 and rsh0 from the slopes of the provided iv_curve
rs0 = rs_slope(pv_curve, r_avg);
rsh0 = rsh_slope(pv_curve, r_avg);

intern.rs0 = rs0;
intern.rsh0 = rsh0;


%% Prepare equations to be used

%  To solve the equations in their general form, we need the variables
%  declared as symbolic variables.
syms n il i0 rs rsh

%  Equation 15 from the article: Solve for i0 and create a matlab function
eq15 = @(n,i0,il,rs,rsh) (il - (v_mp+i_mp*rs)/(rsh) - i_mp) / (exp((v_mp+i_mp*rs)/(n*T)) - 1) - i0;
eq15_sym = solve(eq15,i0);
eq15_i0 = matlabFunction(eq15_sym);

%  Equation 11 from the article: Solve for il and create a matlab function
eq11 = @(n,i0,il,rsh) i_sc - il + i0*(exp((i_sc*rs)/(n*T)) - 1) + (i_sc*rs)/rsh;
eq11_sym = solve(eq11,il);
eq11_il = matlabFunction(eq11_sym);

%  Equation 12 from the article: Solve for rsh and create a matlab function
eq12 = @(n,i0,rs,rsh) ((i0/(n*T))*exp((i_sc*rs)/(n*T)) + 1/rsh) / (1 + rs*((i0/(n*T))*exp((i_sc*rs)/(n*T)) + 1/rsh )) - 1/rsh0;
eq12_sym = solve(eq12,rsh);
eq12_rsh = matlabFunction(eq12_sym);

%  Equation 13 from the article: Solve for n and create a matlab function
eq13 = @(n,i0,il,rsh) il - i0*(exp(v_oc/(n*T)) - 1) - v_oc/rsh;
eq13_sym = solve(eq13,n);
eq13_n = matlabFunction(eq13_sym);

%  Equation 14 from the article: Solve for rs and create a matlab function
eq14 = @(n,i0,rs,rsh) ((i0/(n*T))*exp(v_oc/(n*T)) + 1/rsh) / (1 + rs*((i0/(n*T))*exp(v_oc/(n*T)) + 1/rsh)) - 1/rs0;
eq14_sym = solve(eq14,rs);
eq14_rs = matlabFunction(eq14_sym);


%% Use the equations to determine the parameters based on a double try-and-error procedure

clear n i0 il rs rsh;

%  Set initial values for the parameters
il = i_sc;
rsh = rsh0;
rs = rs0;
n = n_guess;  % not the same n as usual. Here n' = n*k/q (could be updated for conformity)


%  An array to save the error development of rs
diff_rs = [];

%  Loop for adjustments of rs (Could be based on improvement instead)
for j=1:20

    % Array to save the error development of n
    diff_n = [];

    % Loop for adjustments of n
    for i=1:100

        % Calculate i0, il and rsh based on the current rs and n
        i0 = eq15_i0(il,n,rs,rsh);
        il = eq11_il(i0,n,rs,rsh);
        rsh = eq12_rsh(i0,n,rs);

        % Calculate n and compare it to our current estimate
        n1 = eq13_n(i0,il,rsh);
        diff_n = [diff_n; abs(1-n/n1)];

        % Update n for next iteration
        n = n1;

    end

    % Calculate rs and compare it to our current estimate
    rs1 = eq14_rs(i0,n,rsh);
    diff_rs = [diff_rs; abs(1-rs/rs1)];

    % Update rs for the next iteration
    rs = rs1;

end

i_0 = i0;
i_pv = il;

end