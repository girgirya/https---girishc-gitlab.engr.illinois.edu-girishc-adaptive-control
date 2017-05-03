function out=projop(theta,y,thetamax,eps)

ftheta=(theta'*theta-thetamax^2)/(eps*thetamax^2);
gradftheta=2*theta/(eps*thetamax^2);
if ftheta < 0
    out=y;
elseif ftheta>=0 &&  gradftheta'*y<=0
    out=y;
elseif ftheta>=0 &&  gradftheta'*y>0
    out=y-gradftheta/norm(gradftheta)*(gradftheta'/norm(gradftheta)*y)*ftheta;
end
