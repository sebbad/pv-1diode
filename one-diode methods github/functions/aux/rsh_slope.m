%% rsh_slope
%
%  Function to calculate rsh0 from the dI/dV slope at short circuit
%  condition. The function expects a matrix of IV data points (v,i) with
%  the first entry in the matrix corresponding to the short circuit
%  condition. The function, moreover, can calculate multiple slopes and average them.
%
%  Params:
%    iv_data    - an iv curve in matrix format containing voltage values in
%                 column 1 and current values in column 2.
%    num_points - the number of slopes to average
%
%  Returns:
%    rsh0        - the calculated rsh0 value.
%

function rsh0 = rsh_slope(iv_data, num_points)

buf = (iv_data(2:(num_points+1),2) - iv_data(1,2)) ./ (iv_data(2:(num_points+1),1) - iv_data(1,1));

buf_mean = mean(buf);

rsh0 = -1/buf_mean;

end