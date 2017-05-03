function out=projop(theta,y,thetamax,eps)
% Projection operator for bounding adaptive control weights
%outputs a bounded thetadot,
%input: 
%theta: parameter vector (adaptive weights)
%y: thetadot
%thetamax: upper limit of theta
%eps: parameter for projection operator, nonzero typically small of the
%order of 0.001
ftheta=(theta'*theta-thetamax^2)/(eps*thetamax^2);
gradftheta=2*theta/(eps*thetamax^2);
if ftheta < 0
    out=y;
elseif ftheta>=0 &&  gradftheta'*y<=0
    out=y;
elseif ftheta>=0 &&  gradftheta'*y>0
    out=y-gradftheta/norm(gradftheta)*(gradftheta'/norm(gradftheta)*y)*ftheta;
end
