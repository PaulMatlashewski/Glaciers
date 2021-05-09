function [fout, Dfout] = lakeevolve(v,parameters)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%simulation mode
if isfield(parameters,'fluxadjust')
    fluxadjust = parameters.fluxadjust;
else
    fluxadjust = false;
end
if isfield(parameters,'voladjust')
    voladjust = parameters.voladjust;
else
    voladjust = false;
end

if voladjust
    q0 = kappafun(1,1,1,parameters);
end

%grid information
ops = parameters.ops;
npts = parameters.npts;

%operators
grad = ops.grad;
div = ops.div;
mean = ops.mean;
grad_staggered = ops.grad_staggered;
div_staggered = ops.div_staggered;
mean_staggered = ops.mean_staggered;

%time step
dt = parameters.dt;
hprev = parameters.hprev;
Nprev = parameters.Nprev;

%other parameters
epsilon = parameters.epsilon;
deltainv = parameters.deltainv;
nu = parameters.nu;
r = parameters.r;
slope = parameters.slope;

%unpack variables
u = v(1:npts);              %u on staggered grid
h = v(npts+1:2*npts);       %h on original grid
N = v(2*npts+1:3*npts);     %N on original grid
if fluxadjust
    a = v(3*npts+1);
end
if voladjust
    if ~fluxadjust
        m = v(3*npts+1);
    else
        m = v(3*npts+2);
    end
end

%Velocity solve
sigma = 2*h.*grad_staggered*u;
meanN = mean*N;
meanh = mean*h;
gradh = grad*h;
if nargout > 1
    [tau, dtaudu, dtaudmeanN] = taufun(u,meanN,parameters);
    [hw, dhwdN] = hwfun(N,parameters);
else
    tau =  taufun(u,meanN,parameters);
    hw = hwfun(N,parameters);
end
gradhw = grad*hw;
taud =  meanh.*(slope - deltainv*gradh - nu*gradhw);
fout1 = 2*epsilon*div_staggered*sigma - tau + taud;
if nargout > 1
    dsigmadu = spdiags(2*h,0,npts,npts)*grad_staggered;
    dsigmadh = spdiags(sigma./h,0,npts,npts);
    dtauddh =  spdiags(slope - deltainv*gradh - nu*gradhw,0,npts,npts)*mean...
        + spdiags(meanh,0,npts,npts)*(-deltainv*grad);
    dtauddN = spdiags(meanh,0,npts,npts)*(-nu*grad)*spdiags(dhwdN,0,npts,npts);
    dfout1du = 2*epsilon*div_staggered*dsigmadu - spdiags(dtaudu,0,npts,npts);
    dfout1dh = 2*epsilon*div_staggered*dsigmadh + dtauddh;
    dfout1dN = - spdiags(dtaudmeanN,0,npts,npts)*mean + dtauddN;
end
if fluxadjust
    dfout1da = sparse(npts,1);
end
if voladjust
    dfout1dm = sparse(npts,1);
end
%Thickness solve
if fluxadjust
    fout2 = [deltainv*(h - hprev)/dt + div*(u.*meanh) - a; sum(u.*meanh)/npts-1];
else
    fout2 = deltainv*(h - hprev)/dt + div*(u.*meanh);
end
if nargout > 1
    if fluxadjust
        dfout2dh = [deltainv*speye(npts)/dt + div*spdiags(u,0,npts,npts)*mean; sparse(ones(npts,1),(1:npts).',u/npts,1,npts)];
        dfout2du = [div*spdiags(meanh,0,npts,npts); sparse(ones(npts,1),(1:npts).',meanh/npts,1,npts)];
        dfout2dN = sparse(npts+1,npts);
        dfout2da = -[sparse((1:npts).',ones(npts,1),ones(1,npts),npts,1);sparse(1,1)];
    else
        dfout2dh = deltainv*speye(npts)/dt + div*spdiags(u,0,npts,npts)*mean;
        dfout2du = div*spdiags(meanh,0,npts,npts);
        dfout2dN = sparse(npts,npts);
    end
end
if voladjust
    if fluxadjust
        dfout2dm = sparse(npts+1,1);
    else
        dfout2dm = sparse(npts,1);
    end
end

%water pressure solve
hwprev = hwfun(Nprev,parameters);
gradN = grad*N;
Psi = slope - r*deltainv*gradh + gradN;
if nargout > 1
    [kappa,dkappadmeanN,dkappadPsi,dkappadu] = kappafun(meanN,Psi,u,parameters);
else
    kappa = kappafun(meanN,Psi,u,parameters);
end
q = kappa.*Psi;
if voladjust
    fout3 = [(hw-hwprev)/dt + div*q- m;  sum(q)/npts-q0];
else
    fout3 = (hw-hwprev)/dt + div*q;
end
if nargout > 1
    if voladjust
        dPsidh = -r*deltainv*grad;
        dPsidN = grad;
        dkappadN = spdiags(dkappadmeanN,0,npts,npts)*mean + spdiags(dkappadPsi,0,npts,npts)*dPsidN;
        dkappadh =  spdiags(dkappadPsi,0,npts,npts)*dPsidh;
        dqdN = spdiags(kappa,0,npts,npts)*dPsidN + spdiags(Psi,0,npts,npts)*dkappadN;
        dqdh = spdiags(kappa,0,npts,npts)*dPsidh + spdiags(Psi,0,npts,npts)*dkappadh;
        dqdu = spdiags(Psi,0,npts,npts)*spdiags(dkappadu,0,npts,npts);
        dqmeandN = ones(1,npts)*dqdN/npts;
        dqmeandh = ones(1,npts)*dqdh/npts;
        dqmeandu = ones(1,npts)*dqdu/npts;
        dfout3dN = [spdiags(dhwdN/dt,0,npts,npts) + div*dqdN; sparse(ones(npts,1),(1:npts).',dqmeandN,1,npts)];
        dfout3dh = [div*dqdh; sparse(ones(npts,1),(1:npts).',dqmeandh,1,npts)];
        dfout3du = [div*dqdu; sparse(ones(npts,1),(1:npts).',dqmeandu,1,npts)];    
        dfout3dm = -[sparse((1:npts).',ones(npts,1),ones(1,npts),npts,1);sparse(1,1)];
    else
        dPsidh = -r*deltainv*grad;
        dPsidN = grad;
        dkappadN = spdiags(dkappadmeanN,0,npts,npts)*mean + spdiags(dkappadPsi,0,npts,npts)*dPsidN;
        dkappadh =  spdiags(dkappadPsi,0,npts,npts)*dPsidh;
        dqdN = spdiags(kappa,0,npts,npts)*dPsidN + spdiags(Psi,0,npts,npts)*dkappadN;
        dqdh = spdiags(kappa,0,npts,npts)*dPsidh + spdiags(Psi,0,npts,npts)*dkappadh;
        dqdu = spdiags(Psi,0,npts,npts)*spdiags(dkappadu,0,npts,npts);
        dfout3dN = spdiags(dhwdN/dt,0,npts,npts) + div*dqdN;
        dfout3dh = div*dqdh;
        dfout3du = div*dqdu;
    end
end
if fluxadjust
    if voladjust
        dfout3da = sparse(npts+1,1);
    else
        dfout3da = sparse(npts,1);
    end
end

%assemlbe output  
fout = [fout1; fout2; fout3];
if fluxadjust
    if voladjust
        %disp('fout1')
        %size(dfout1du), size(dfout1dh), size(dfout1dN), size(dfout1da), size(dfout1dm)
        %disp('fout2')
        %size(dfout2du), size(dfout2dh), size(dfout2dN), size(dfout2da), size(dfout2dm)
        %disp('fout3')
        %size(dfout3du), size(dfout3dh), size(dfout3dN), size(dfout3da), size(dfout3dm)
        Dfout = [dfout1du, dfout1dh, dfout1dN, dfout1da, dfout1dm;...
            dfout2du, dfout2dh, dfout2dN, dfout2da, dfout2dm;
            dfout3du, dfout3dh, dfout3dN, dfout3da, dfout3dm];    
    else
        Dfout = [dfout1du, dfout1dh, dfout1dN, dfout1da;...
            dfout2du, dfout2dh, dfout2dN, dfout2da;
            dfout3du, dfout3dh, dfout3dN, dfout3da]; 
    end
elseif voladjust
        Dfout = [dfout1du, dfout1dh, dfout1dN, dfout1dm;...
            dfout2du, dfout2dh, dfout2dN, dfout2dm;
            dfout3du, dfout3dh, dfout3dN, dfout3dm]; 
else
    Dfout = [dfout1du, dfout1dh, dfout1dN;...
        dfout2du, dfout2dh, dfout2dN;
        dfout3du, dfout3dh, dfout3dN];  
end

end

function [fout, Dfu, DfN] = taufun(u,N,parameters)
    a = parameters.a;
    b = parameters.b;
    C = parameters.C;
    fout = C*abs(u).^(a-1).*u.*max(N,0).^b;
    if nargout > 1
        Dfu = a*fout./u;
        Dfu(u==0) = 0;      %deal with indeterminate case
        DfN = b*fout./N;
        DfN(N==0) = 0;      %deal with indeterminate case
    end
end

function [fout, DfN] = hwfun(N,parameters)
    alpha = parameters.alpha;
    gamma = parameters.gamma;
    nu = parameters.nu;
    fout = exp(-alpha*N);
    fout(N<0) = gamma*(1 -N(N<0)/nu);
    if nargout > 1
        DfN = -alpha*gamma*exp(-alpha*N);
        DfN(N<0) = -gamma/nu;
    end
end

function [fout, DfN, Dfpsi, Dfu] = kappafun(N,psi,u,parameters)
    beta = parameters.beta;
    kappa0 = parameters.kappa0;
    akappa = parameters.akappa;
    bkappa = parameters.bkappa;
    ckappa = parameters.ckappa;
    mu = parameters.mu;
    fout = kappa0*exp(-beta*max(N,0)).*abs(psi).^akappa.*abs(u).^bkappa + max(0,-N).^ckappa/mu;
    fout(fout > 1e8) = 1e8;
    if nargout > 1
        DfN = -beta*kappa0*exp(-beta*N).*abs(psi).^akappa.*abs(u).^bkappa;
        DfN(N < 0) = -ckappa*(-N(N<0)).^(ckappa-1)/mu;
        Dfpsi = akappa*kappa0*exp(-beta*max(N,0)).*abs(psi).^(akappa-1).*abs(u).^bkappa.*sign(psi);
        Dfu =  bkappa*kappa0*exp(-beta*max(N,0)).*abs(psi).^akappa.*abs(u).^(bkappa-1).*sign(u);
        DfN(fout > 1e8) = 0;
        Dfpsi(fout > 1e8) = 0;
        Dfu(fout > 1e8) = 0;
    end
end
