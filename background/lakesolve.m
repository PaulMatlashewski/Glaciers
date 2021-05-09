function fsolution = lakesolve(parameters,nsteps,nskip,test)
%fsolution = lakesolve(parameters,nsteps,nskip)
%   Detailed explanation goes here

%create parameter structure if necessary
if nargin == 0 || ~isfield(parameters,'npts')
    parameters.npts = 800;
end
if ~isfield(parameters,'dt')
    parameters.dt = 1e-2;
end
if ~isfield(parameters,'dx')
    parameters.dx = .1;
end
if ~isfield(parameters,'ops')
   parameters.ops = grid(parameters.npts,parameters.dx);
end
if ~isfield(parameters,'epsilon')
    parameters.epsilon = 1;
end
if ~isfield(parameters,'deltainv')
    parameters.deltainv = 1e2;
end
if ~isfield(parameters,'gamma')
    parameters.gamma = 1;
end
if ~isfield(parameters,'nu')
    parameters.nu = 5e-2;
end
if ~isfield(parameters,'r')
    parameters.r = .9;
end
if ~isfield(parameters,'slope')
    parameters.slope = 1;
end
if ~isfield(parameters,'a')
    parameters.a = 1;
end
if ~isfield(parameters,'b')
    parameters.b = 1;
end
if ~isfield(parameters,'C')
    parameters.C = 1;
end
if ~isfield(parameters,'gamma')
    parameters.gamma = 1;
end
if ~isfield(parameters,'alpha')
    parameters.alpha = 1;
end
if ~isfield(parameters,'kappa0')
    parameters.kappa0 = 1;
end
if ~isfield(parameters,'beta')
    parameters.beta = 2.5e-1;
end
if ~isfield(parameters,'mu')
    parameters.mu = 1e-6;
end
if ~isfield(parameters,'ckappa')
    parameters.ckappa = 3;
end
if ~isfield(parameters,'akappa')
    parameters.akappa = 0;
end
if ~isfield(parameters,'bkappa')
    parameters.bkappa = 0;
end
if ~isfield(parameters,'fluxadjust')
    parameters.fluxadjust = false;
end
if ~isfield(parameters,'voladjust')
    parameters.voladjust = false;
end

if nargin < 2 || isempty(nsteps)
    nsteps = 1000;
end

if nargin < 3 || isempty(nskip)
    nskip = 48*4;
end

%number of points
npts = parameters.npts;

%time step size
dt = parameters.dt;

%solver settings
srchparams.toldelta = 1e-5;
srchparams.tolnorm = 1e-5;


if nargin == 4  && test
    parameters.hprev = rand(npts,1);
    parameters.Nprev = rand(npts,1);
    vin = 2*rand(3*npts,1)-1;
    if parameters.fluxadjust
        vin = [vin; rand(1)];
    end
    if parameters.voladjust
        vin = [vin; rand(1)];
    end
    fsolution = Jacobian_test_single(@lakeevolve,vin,parameters);
    fsolution.parameters = parameters;
    return
else
    if isfield(parameters,'hprev') && isfield(parameters,'Nprev')
        vin = [ones(npts,1);parameters.hprev;parameters.Nprev];
    else
        vin = [ones(npts,1); 1+1e-6*rand(npts,1); 1+1e-2*rand(npts,1)];
    end
    if parameters.fluxadjust
        if isfield(parameters,'aprev')
            vin = [vin; parameters.aprev];
        else
            vin = [vin;1];
        end
    end
    if parameters.voladjust
        if isfield(parameters,'mprev')
            vin = [vin; parameters.mprev];
        else
            vin = [vin;1];
        end
    end
end


%initialize
ii = 1;
jj = 0;
t = 0;
fsolution.parameters = parameters;
fsolution.u = zeros(npts,nsteps);
fsolution.h = zeros(npts,nsteps);
fsolution.N = zeros(npts,nsteps);
if parameters.fluxadjust
    fsolution.a = zeros(nsteps,1);
end
if parameters.voladjust
    fsolution.m = zeros(nsteps,1);
end
fsolution.t = zeros(nsteps,1);

%initial condition
parameters.hprev = vin(npts+1:2*npts);
parameters.Nprev = vin(2*npts+1:3*npts);

% run through time stepping
while ii < nsteps
    vout = Newton_single(@lakeevolve,vin,parameters,srchparams);
    if ii==1 && jj == 0
        fsolution.u(:,ii) = vout(1:npts);
        fsolution.h(:,ii) = vin(npts+1:2*npts);
        fsolution.N(:,ii) = vin(2*npts+1:3*npts);
        if parameters.fluxadjust
            fsolution.a(ii) = vin(3*npts+1);
            if parameters.voladjust
                fsolution.m(ii) = vin(3*npts+2);
            end
        elseif parameters.voladjust
            fsolution.m(ii) = vin(3*npts+1);
        end
        fsolution.t(ii) = t;
    end
    t = t+dt;
    parameters.hprev = vout(npts+1:2*npts);
    parameters.Nprev = vout(2*npts+1:3*npts);
    vin = vout;
    if jj < nskip
        jj = jj+1;
    else
        jj = 0;
        ii = ii+1
        fsolution.u(:,ii) = vout(1:npts);
        fsolution.h(:,ii) = vout(npts+1:2*npts);
        fsolution.N(:,ii) = vout(2*npts+1:3*npts);
        fsolution.t(ii) = t;
    end
end

end

