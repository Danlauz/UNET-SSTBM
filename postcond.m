function datasim=postcond(x,xsim,x0sim,nbsimul,ki,k0)
% function for carry out the post-conditionning (mean 0)
% syntax: datasim=postcond(x,xsim,x0sim,nbsimul,ki,k0)
%
% Input:
%   x: matrix n x 3 of data points (x,y,z(x,y))%
%   xsim: matrix n x (2+nbsimul) ?? simulated values non-conditionally at data points
%   x0sim: matrix nx0 x (2+nbsimul) simulated points not conditionally at points x0
%   model, c:  like cokri
%   nbsimul : number of realization
%   k : matrix of covariance of data points x
%   k0 : matrix of covariance of data points x with data points xsim
%
% Output:
%   datasim: matrix nx0 x (2+nbsimul) simulated values conditionally at points x0
% ATTENTION
%   that program assumes that the mean is 0. Then x, xsim, x0sim must match with fields of mean 0.
%
% Modification: D. Lauzon, 2018
% Autor: D. Marcotte, 2004
b0=ki*x(:,end);
datasim=x0sim;

if size(x,2)==4
bi=ki*xsim(:,4:end);
datasim(:,4:end)=((b0*ones(1,nbsimul)-bi)'*k0)'+x0sim(:,4:end);
end

if size(x,2)==3
bi=ki*xsim(:,3:end);
datasim(:,3:end)=((b0*ones(1,nbsimul)-bi)'*k0)'+x0sim(:,3:end);
end