%% Example runs of the one-diode-model methods included
%
%  The script uses data provided in Ghani et al. 2013 for a 1 cell PV
%  device. For each method, the five model parameters are extracted given
%  a number of assumptions required by the individual methods.
%
%  Once the parameters are extracted, an IV curve is predicted using the
%  'iv' auxilliary function provided. This predicted IV curve is then
%  compared to the experimental data provided in 'ghani_data.mat'.
%

clearvars

%  Add path to the respective model functions
addpath('../functions/');
addpath('../functions/aux/');

%  Data for evaluation
load('ghani_data.mat');

data = ghani_experimental;


%% Set initial constants
%  Remarkable points of the targeted PV module (Ghani2013, Table B1).
v_mp = 0.4245;                  % mpp voltage
i_mp = 0.51;                    % mpp current
v_oc = 0.5969;                  % open circuit voltage
i_sc = 0.5768;                  % short circuit current
N = 1;                          % number of cells in series

k = 1.38064852E-23;             % boltzmann constant
q = 1.60217662E-19;             % electron charge

%  Aggregate data to pass to method implementations
pv_data = [v_oc;i_sc;v_mp;i_mp;N];


%% Celik2007

r_averages = 5;

[celik_intern, celik_n, celik_i_0, celik_i_pv, celik_rs, celik_rsh] = celik2007(pv_data, data, r_averages, 1);
try
    [celik_v, celik_i] = iv(celik_n, celik_i_0, celik_i_pv, celik_rs, celik_rsh, N, data(:,1));
    celik_err = celik_i - data(1:length(celik_i),2);
    celik_mae = mean(abs(celik_err));
    celik_mape = mean(abs(celik_err./ghani_experimental(1:length(celik_i),2)));
    celik_rms = rms(celik_err);
catch
end


%% Villalva2009

n_guess = 1.5;
rs_range = linspace(0,1,1000);
tol = 1e-7;

[villalva_intern, villalva_n, villalva_i_0, villalva_i_pv, villalva_rs, villalva_rsh] = villalva2009(pv_data, rs_range, n_guess, tol);
try
    [villalva_v, villalva_i] = iv(villalva_n, villalva_i_0, villalva_i_pv, villalva_rs, villalva_rsh, N, data(:,1));
    villalva_err = villalva_i - data(1:length(villalva_i),2);
    villalva_mae = mean(abs(villalva_err));
    villalva_mape = mean(abs(villalva_err./data(1:length(villalva_i),2)));
    villalva_rms = rms(villalva_err);
catch
end


%% Brano2010

r_averages = 5;
n_guess = 1 * (k*N)/q;  % TODO: adjustment should be included in the function instead

[brano_intern, brano_n, brano_i_0, brano_i_pv, brano_rs, brano_rsh] = brano2010(pv_data, data, r_averages, n_guess);
brano_n = brano_n*q/(k*N);
try
    [brano_v, brano_i] = iv(brano_n, brano_i_0, brano_i_pv, brano_rs, brano_rsh, N, data(:,1));
    brano_err = brano_i - data(1:length(brano_i),2);
    brano_mae = mean(abs(brano_err));
    brano_mape = mean(abs(brano_err./data(1:length(brano_i),2)));
    brano_rms = rms(brano_err);
catch
end


%% Ghani2013

rs_guess = 0.1;
rsh_guess = 10;
tol_r = 1e-7;
tol_n = 1e-3;

[ghani_intern, ghani_n, ghani_i_0, ghani_i_pv, ghani_rs, ghani_rsh] = ghani2013(pv_data, rs_guess, rsh_guess, tol_r, tol_n);
try
    [ghani_v, ghani_i] = iv(ghani_n, ghani_i_0, ghani_i_pv, ghani_rs, ghani_rsh, N, data(:,1));
    ghani_err = ghani_i - data(1:length(ghani_i),2);
    ghani_mae = mean(abs(ghani_err));
    ghani_mape = mean(abs(ghani_err./data(1:length(ghani_i),2)));
    ghani_rms = rms(ghani_err);
catch
end


%% Cubas2014

rs_guess = 1;

[cubas_intern, cubas_n, cubas_i_0, cubas_i_pv, cubas_rs, cubas_rsh] = cubas2014(pv_data, rs_guess);
try
    [cubas_v, cubas_i] = iv(cubas_n, cubas_i_0, cubas_i_pv, cubas_rs, cubas_rsh, N, data(:,1));
    cubas_err = cubas_i - data(1:length(cubas_i),2);
    cubas_mae = mean(abs(cubas_err));
    cubas_mape = mean(abs(cubas_err./data(1:length(cubas_i),2)));
    cubas_rms = rms(cubas_err);
catch
end


%% Laudani2014

rs_guess = 0.1;
n_range = linspace(0.4,2,1000);

[laudani_intern, laudani_n, laudani_i_0, laudani_i_pv, laudani_rs, laudani_rsh] = laudani2014(pv_data, rs_guess, n_range);
try
    [laudani_v, laudani_i] = iv(laudani_n, laudani_i_0, laudani_i_pv, laudani_rs, laudani_rsh, N, data(:,1));
    laudani_err = laudani_i - data(1:length(laudani_i),2);
    laudani_mae = mean(abs(laudani_err));
    laudani_mape = mean(abs(laudani_err./data(1:length(laudani_i),2)));
    laudani_rms = rms(laudani_err);
catch
end


%% Hussein2017

rs_guess = 0.1;
n_start = 1;
tol = 1e-7;

[hussein_intern, hussein_n, hussein_i_0, hussein_i_pv, hussein_rs, hussein_rsh] = hussein2017(pv_data, rs_guess, n_start, data, tol);
try
    [hussein_v, hussein_i] = iv(hussein_n, hussein_i_0, hussein_i_pv, hussein_rs, hussein_rsh, N, data(:,1));
    hussein_err = hussein_i - data(1:length(hussein_i),2);
    hussein_mae = mean(abs(hussein_err));
    hussein_mape = mean(abs(hussein_err./data(1:length(hussein_i),2)));
    hussein_rms = rms(hussein_err);
catch
end

