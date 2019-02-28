%% rs_slope
%  
%  Function to calculate rs0 from the dI/dV slope at open circuit
%  condition. The function expects a matrix of IV data points (v,i) with the last
%  entry in the matrix corresponding to the open circuit condition. The
%  function, moreover, can calculate multiple slopes and average them.
%
%  Params:
%    iv_data    - an iv curve in matrix format containing voltage values in
%                 column 1 and current values in column 2.
%    num_points - the number of slopes to average
%
%  Returns:
%    rs0        - the calculated rs0 value.
%


function rs0 = rs_slope(iv_data, num_points)

buf = (iv_data(end-num_points:end-1,2) - iv_data(end,2)) ./ (iv_data(end-num_points:end-1,1) - iv_data(end,1));

buf_mean = mean(buf);

rs0 = -1/buf_mean;

end