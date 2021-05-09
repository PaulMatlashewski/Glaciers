function [xout, error_flag, faux] = Newton_single(f,x,parameters,srchparams)
%[xout, error_flag, faux] = Newton_single(f,x,parameters,srchparams)
%Newton: Newton's method for finding zeros of function f. Works the same as
%Newton_v2 but calls only a single function f that outputs the function
%vector and its derivative matrix at the same time; useful if it is cheaper
%to compute function and its derivatives simultaneously
%
%   Input arguments:
%   f:      function handle of function whose roots are to be found, along
%           with its derivative with respect to x
%           If Newton_single requires one or two output arguments
%           f must take the form
%               [fout, Dfout] = f(x,parmaeters)
%           If Newton_single requires three output arguments, f must take the form
%               [fout, Dfout, faux] = f(x,parameters)
%           The Newton iteration will solve for fout=0; faux can be
%           additional outputs computed within f that may be used elsewhere
%   x:      initial guess; must have correct dimensions for an input
%           argument of f
%   parameters: optional parameter structure to be passed to f, df. If
%           parameters is not passed to Newton_v2, an empty array is passed
%           to f, df instead. If parameters os passed to Newton_v2 but has
%           a field 'noparameters' set to 'true', no parameters structure
%           is passed to f, df.
%   srchparams: optional structure containing search parameters for Newton
%               iteration. This can contain the fields:
%               itmax:  maximum number of iterations
%               toldelta:   error bound on x; last Newton step must be <=
%                           0.5*toldelta
%               tolF:    error bound on f = 0; 0.5*norm(f,2) <=
%                           tolF to terminate iteration
%               verbose:    0 produces no runtime output
%                           1 displays step size and current value of
%                           function f
%                           2 also displays current iterate
%
%   The algorithm used is Newton's method as a default; intend to upgrade
%   to line search algorithm
%
%   Output:
%   xout:       root of f
%   error_flag: returns logical one if convergence to prescribed tolerance
%               has not been achieved
%   faux:       optional second output of function f, if auxiliary
%               calculations within f generate additional output
%
%Christian Schoof, January 2014

%NOTE: Check against Numerical Recipes for more efficient version; see
%below for comments regarding line search.

    %set up Newton iteration parameters
    if nargin == 4
        if isfield(srchparams,'itmax')
            itmax = srchparams.itmax;
        else
            itmax = 25;
        end
        if isfield(srchparams,'tolF')
            tolF = srchparams.tolF;
        else
            tolF = sqrt(eps);
        end
        if isfield(srchparams,'toldelta')
            toldelta = srchparams.toldelta;
        else
            toldelta = sqrt(eps);
        end
        if isfield(srchparams,'verbose')
            verbose = srchparams.verbose;
        else
            verbose = 0;
        end
    else
        itmax = 25;
        tolF = sqrt(eps);
        toldelta = sqrt(eps);
        verbose = 0;
    end

    %Deal with missing input parameters
    if nargin < 3 || isempty(parameters)
        parameters.noparameters = true;
    end
    
    %set error_flag to false if required
    if nargout > 1, error_flag = false; end
    
    %set up first Newton iteration step
    if nargout < 3
        if isfield(parameters,'noparameters') && parameters.noparameters
            [F, DF] = f(x);
        else
            [F, DF] = f(x,parameters);
        end
    else
        if isfield(parameters,'noparameters') && parameters.noparameters
            [F, DF, faux] = f(x);
        else
            [F, DF, faux] = f(x,parameters);
        end
    end
    if verbose, disp(strcat('Initial norm of f =',num2str(norm(F,2)))), end
    Dx = toldelta + eps;
    iter = 1;
    while iter < itmax && (norm(F) >= tolF || norm(Dx) >= toldelta)
        Dx = -DF\F;
        x = x + Dx;
        if nargout < 3
            if isfield(parameters,'noparameters') && parameters.noparameters
                [F, DF] = f(x);
            else
                [F, DF] = f(x,parameters);
            end
        else
            if isfield(parameters,'noparameters') && parameters.noparameters
                [F, DF, faux] = f(x);
            else
                [F, DF, faux] = f(x,parameters);
            end
        end
        if verbose, disp(strcat('Updated norm of f =',num2str(norm(F,2)))), disp(strcat('Size of step Dx =',num2str(norm(Dx,2)))), end
        iter = iter + 1;
        if verbose == 2, disp('Current iterate ='), disp(x), end
    end
    %output warning, set error_flag if convergence not achieved
    if iter >= itmax, warning(strcat('Newton did not converge, current norm of f = ',num2str(norm(F,2)))), if nargout == 2, error_flag = true; end, end
    xout = x;
end