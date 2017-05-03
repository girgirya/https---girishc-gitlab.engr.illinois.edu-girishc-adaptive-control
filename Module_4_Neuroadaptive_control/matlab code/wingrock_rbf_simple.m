%wingrock simulation with matched uncertainty

close all;
clear all;

%% sim params
t0=0;
tf=67;%170
dt=0.05;
t=t0:dt:tf;


%% initialization

x=[1.5;1.5];



Wstar=[0.8 0.2314 0.6918 -0.6245 0.0095 0.0214]';


%baeline linear design
Kp = 1.5;                 % proportional gain
Kd = 0.6;                 % derivative gain

A=[0 1;-Kp -Kd];
B = [0; 1];
Q = eye(2);
P = lyap(A',Q);%A' instead of A because need to solve A'P+PA+Q=0
 
[Lp,PP,E] = lqr(A,B,eye(2)*100,1);

%% RBF NN parameters
n2=10;
phi_space=[-2,2];
phid_space=[-2,2];
%Distribute centers using a uniform random distribution
for ii = 2 :n2+1
    rbf_c(1,ii)=unifrnd(phi_space(1),phi_space(2));
    rbf_c(2,ii)=unifrnd(phid_space(1),phid_space(2));
    rbf_mu(ii)=1;
end
W=zeros(n2+1,1);
sigma=zeros(n2+1,1);
bw=1;
sigma(1)=bw;


%% adaptive control initialization
gammaW=2;%learning rate
kappa=0.00; %sigma mod
zeta=0.2; %emod

%% commands
xref=0;
XREF=zeros(length(t),1);
 XREF(5/dt:7/dt)=0;
 XREF(15/dt:17/dt)=1;
 XREF(25/dt:27/dt)=-1;

%% reference model parameters
omegan_rm = 1;        % reference model natural freq (rad/s)
zeta_rm   = 1;        % reference model damping ratio
x_rm = x;
v_h=0;
v_ad=0;


%% plotting array initialization 
index=1;
t=t0:dt:tf;
T_REC        = zeros(length(t),1);
X_REC        = zeros(length(t),1);
XDOT_REC     = zeros(length(t),1);
XRM_REC      = zeros(length(t),1);
XDOTRM_REC   = zeros(length(t),1);
XERR_REC     = zeros(length(t),1);
XDOTERR_REC  = zeros(length(t),1);
DELTACMD_REC = zeros(length(t),1);%control input
DELTAERR_REC = zeros(length(t),1);
W_REC = zeros(n2+1,length(t));
NU_AD_REC = zeros(length(t),1);

%%
for t=t0:dt:tf
   
    e=x_rm-x;%compute reference model error
    xref=XREF(index); 
    
    v_crm=omegan_rm^2*(xref-x_rm(1))-2*zeta_rm*omegan_rm*x_rm(2);
    v_pd=[Kp Kd]*e;
    deltaCmd = v_crm + v_pd - v_ad;%Nu
    delta = deltaCmd;
    v_h=deltaCmd-delta;
          
    [x,x_rm,xddot,deltaErr,v_crm]=wingrock_correct(x,x_rm,v_h,delta,dt,dt,Wstar,xref,omegan_rm,zeta_rm);

   sigma= sdash_rbf(x,rbf_c,rbf_mu,bw,n2);

    
      Wd=-gammaW*e'*P*B*sigma-gammaW*kappa*W-gammaW*zeta*norm(e)*W;
    

    Wdot=Wd;
    W=W+dt*Wdot;
   
    
  
    %assign adaptive control
    v_ad=W'*sigma;
    
    
    
         
     %
    %record for plotting
    T_REC(index)        = t;
    X_REC(index)        = x(1);
    XDOT_REC(index)     = x(2);
    XRM_REC(index)     = x_rm(1);
    XDOTRM_REC(index)   = x_rm(2);
    XERR_REC(index)     = e(1);
    XDOTERR_REC(index)  = e(2);
    DELTACMD_REC(index) = delta;%control input
    DELTAERR_REC(index) = deltaErr;
    W_REC(:,index)         =W;
    NU_AD_REC(index)    = v_ad;
    
    PHI(index,:)=sigma(2:n2+1)';
    index = index+1;
end

%% plotting
figure(1);
subplot(2,1,1)
plot(T_REC,X_REC, T_REC,XRM_REC,'--');
xlabel('time (sec)');
ylabel('pi-rad');
title('roll angle');
legend('actual','ref model',0);
grid on;

subplot(2,1,2)
plot(T_REC,XDOT_REC, T_REC,XDOTRM_REC,'--');
xlabel('time (sec)');
ylabel('xDot (pi-rad/s)');
title('roll rate');
legend('actual','ref model',0);
grid on;

figure(2);
subplot(2,1,1);
plot(T_REC, XERR_REC);
xlabel('time (sec)');
ylabel('xErr (pi-rad)');
title('Position Error');
grid on;
subplot(2,1,2);
plot(T_REC, XDOTERR_REC);
xlabel('time (sec)');
ylabel('xDotErr (pi-rad/s)');
title('Angular Rate Error');
grid on;

figure(5)
plot(T_REC,W_REC(1,:),T_REC,W_REC(2,:),T_REC,W_REC(3,:),T_REC,W_REC(4,:),T_REC,W_REC(5,:),T_REC,W_REC(6,:))
legend('bw','W',0)
xlabel('time (sec)');
ylabel('W')

figure(12)
plot(X_REC,XDOT_REC)
grid on
xlabel('roll angle deg')
ylabel('roll rate deg/sec')


figure(13)
plot(T_REC,DELTAERR_REC,T_REC,NU_AD_REC)
xlabel('time')
ylabel('model error')
legend('model error', 'NN estimate',0)

