function ops = grid(npts,dx)
%ops = grid(npts,dx)
%   Divided difference differential and averaging operators for a
%   one-dimensional periodic staggered grid with constant grid spacing. We assume
%   npts grid points on each grid and define a staggered grid as notionally
%   shifted by 1/2 grid point from index i to index i+1/2. If a vector v is
%   defined on the original grid, then v_i = v(i), whereas a vector
%   v_stagger defined on the staggered grid has v_{i+1/2} = v_stagger(i)
%   The following operators are defined:
%       grad:   if f = grad*v with v defined on the original grid, then
%               f_{i+1/2} = (v_{i+1}-v_i)/dx is defined on the staggered
%               grid
%       div:    if m = div*q with q defined on the staggered grid and m on
%               the original grid, then m_i = (q_{i+1/2}-q_{i-1/2})/dx
%       mean:   if f = mean*h with h defined on the original grid, then
%               f_{i+1/2} = (h_i + h_{i+1})/2 is defined on the staggered
%               grid
%       grad_staggered:   if f = grad*v with v defined on the staggered grid,
%               then f_{i+1} = (v_{i+3/2} - v_{i+1/2})/dx is defined on the
%               original gid
%       div_staggered:    if m = div*q with q defined on the original grid,
%               then m_{i+1/2} = (q_{i+1}-q_i)/dx is defined on the
%               staggered grid. div_stagger replicates grad.
%       mean_staggered:   if f = mean*u with u defined on the staggered grid,
%               then f_i = (u_{i+1/2}-u_{i-1/2})/2 is defined on the
%               original grid
%With a periodic domain with npts grid points, we identify i with i-npts
%for i > npts, and with i+npts for i<1
%
%Input arguments
%       npts:   number of grid points
%       dx:     grid spacing
%Output argument
%       ops:    structure containing the operators grad, div, mean,
%               grad_staggered, div_staggered, mean_staggered as npts-by-npts
%               sparse matrices

ii =1:npts;
jjplus = mod((ii),npts)+1;
jjminus = mod((ii-2),npts)+1;
ops.grad = sparse(ii,ii,-1/dx,npts,npts) + sparse(ii,jjplus,1/dx,npts,npts);
ops.div = sparse(ii,ii,1/dx,npts,npts) + sparse(ii,jjminus,-1/dx,npts,npts);
ops.mean = sparse(ii,ii,.5,npts,npts) + sparse(ii,jjplus,.5,npts,npts);
ops.grad_staggered =  sparse(jjplus,ii,-1/dx,npts,npts) + sparse(jjplus,jjplus,1/dx,npts,npts);
ops.div_staggered = sparse(ii,ii,-1/dx,npts,npts) + sparse(ii,jjplus,1/dx,npts,npts);
ops.mean_staggered = sparse(ii,ii,.5,npts,npts) + sparse(ii,jjminus,.5,npts,npts);

end

