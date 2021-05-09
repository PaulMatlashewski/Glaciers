function xout = Jacobian_test_single(func,z,parameters,delta)
% xout = Jacobian_test_single(func,z,parameters,delta)
%Input:
%func: function handle for function that outputs two arguments, the first
%is the function to differentiate (must be n-by-1 column vector), and the
%second is the analytical Jacobian to be tested (must be n-by-n matrix)
%z: value of input argument to func at which Jacobian is to be tested, must be n-by-1 column vector
%parameters = parameter structure to pass to func and Jacobian; set to empty array if not applicable
%delta: step length in divided difference (optional argument, default is %sqrt(eps)
%Output:
%diff: evaluation of exact Jacobian formula
%diff_num: divided difference approximation
%error: fractional error for every entry in Jacobian, may have NaN's where jacobian is zero
%
%Christian Schoof, January 2014
%corrected bug that overwrites user-specified delta, July 2016
%cosmetic changes, July 2020

if nargin < 4, delta = sqrt(eps); end

if isempty(parameters)
    [~, xout.diff] = func(z);
else
    [~, xout.diff] = func(z,parameters);
end

xout.diff_num = zeros(size(xout.diff));
xout.error = zeros(size(xout.diff));
for p = 1:length(z)
    delta_z = [zeros(p-1,1); 1; zeros(length(z)-p,1)];
    if isempty(parameters)
        [func_plus, ~] = func(z+delta*delta_z);
        [func_0, ~] = func(z);
         xout.diff_num(:,p) = delta^(-1)*(func_plus-func_0);
    else
        [func_plus, ~] = func(z+delta*delta_z,parameters);
        [func_0, ~] = func(z,parameters);
        xout.diff_num(:,p) = delta^(-1)*(func_plus-func_0);
    end
    %compute divided difference approximation to Jacobian
    xout.error(:,p) = abs(xout.diff(:,p)-xout.diff_num(:,p))./(abs(xout.diff_num(:,p)) + abs(xout.diff(:,p)));
    %fractional error (may have NaN's where Jacobian is zero)
    if length(z) < 40
        xout.norm_error = norm(full(xout.diff)-full(xout.diff_num))/norm(full(xout.diff));
    end
    %fractional error in norm, compute only for sufficiently small matrices
end
