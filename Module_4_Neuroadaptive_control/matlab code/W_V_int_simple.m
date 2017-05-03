

%% W integrate RK4 based integration for W and V matrices
function [W,V,r]=W_V_int_simple(gammaW,gammaV,sigma,sigmadash,V,W,xbar,e,P,B,kappa,dt)


r=e'*P*B;


%% propogate W
xp=stat_W(gammaW,sigma,sigmadash,V,xbar,r,kappa,W);

%1

rk1=dt*xp;
x1=W+rk1/2;
%2
xp=stat_W(gammaW,sigma,sigmadash,V,xbar,r,kappa,W);
rk2=dt*xp;
x1=W+rk2/2;

%3
xp=stat_W(gammaW,sigma,sigmadash,V,xbar,r,kappa,W);
rk3=dt*xp;
x1=W+2*rk3;

%4
xp=stat_W(gammaW,sigma,sigmadash,V,xbar,r,kappa,W);
rk4=dt*xp;
W=W+(rk1+2*(rk2+rk3)+rk4)/6;



%% propogate V
clear xp;
xp=stat_V(gammaV,sigmadash,V,xbar,r,kappa,W);

%dt=0.5*dt;


rk1=dt*xp;
x1=V+rk1/2;

%2
xp=stat_V(gammaV,sigmadash,V,xbar,r,kappa,W);


rk2=dt*xp;
x1=V+rk2/2;

%3
xp=stat_V(gammaV,sigmadash,V,xbar,r,kappa,W);

rk3=dt*xp;
x1=V+2*rk3;

%4
xp=stat_V(gammaV,sigmadash,V,xbar,r,kappa,W);
rk4=dt*xp;
V=V+(rk1+2*(rk2+rk3)+rk4)/6;

%% W update model (c1)
function [Wdot]=stat_W(gammaW,sigma,sigmadash,V,xbar,r,kappa,W)


Wd=-gammaW*((sigma-sigmadash*V'*xbar)*r'+kappa*norm(r)*W); %in the ALR the V' term might not be required

Wdot=Wd;

%% V update model (c1)
function[Vdot]=stat_V(gammaV,sigmadash,V,xbar,r,kappa,W)

Vdot=-gammaV(1)*(xbar*r'*W'*sigmadash+kappa*norm(r)*V);


