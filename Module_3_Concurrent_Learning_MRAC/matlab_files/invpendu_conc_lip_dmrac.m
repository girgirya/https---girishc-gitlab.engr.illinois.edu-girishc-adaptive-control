
%inverted pendulum like system simulation with structured uncertainty
%refer to Chowdhary G., Johnson E., Concurrent Learning for Convergence in
%Adaptive Control without Persistency of Excitation, IEEE Conference on Decision and Control, GA, USA 2010.
% Author: Girish Chowdhary, all rights reserved
% email girish.chowdhary@gatech.edu for further information
% last modified: 11/16/2010

close all;
clear all;


%% sim params
t0=0; %intial time
tf=40;%170%final time
dt=0.01; %time step
t=t0:dt:tf;%time interval


%% initialization

x=[1 ;1];%intial state

sdev=0; %noise std deviation if used

n2=3; %dimension of uncertainty
Wstar=[-1 1 0.5]';%ideal weights
W=Wstar*0;%wight initialization

Kp = 1.5;                 % proportional gain
Kd = 1.3;                 % derivative gain

A=[0 1;-Kp -Kd];           %Linear controller
B = [0; 1];                 % B matrix
Q = eye(2);                 % for Lyapunov equation
P = lyap(A',Q);%A' instead of A because need to solve A'P+PA+Q=0
 
[Lp,PP,E] = lqr(A,B,eye(2)*100,1);  % Lyap eqn

%% adaptive control initialization
gammaW=3.5;         % online learning rate
gammaWbkk=3.5;      % learning rate for recorded data updates, set to zero
                    % if don't want concurrent learning
gamWc=0;% gamWc=0 means Wc=I,  (default)
        % Set gamWc= 1 to prioratize training on current data over recorded
        % data
kappa=0;%emod gain set to zero (default)

%% commands are set in this section
xref=0;
XREF=zeros(length(t),1);
 XREF(5/dt:10/dt)=0;
XREF(20/dt:25/dt)=1;
XREF(35/dt:40/dt)=1*0;
XREF(50/dt:55/dt)=1*0;
%% reference model parameters
omegan_rm = 1;        % reference model natural freq (rad/s)
zeta_rm   = 1;        % reference model damping ratio
x_rm = x;              % initialize ref. model state at x(0)
v_h=0;                  %initialize pseudo control
v_ad=0;
%% conc parameter intialization
nbg=1;                  
Xbg=zeros(n2);%the same size as sigma

% X=zeros(max(size(freqs)),n2);


%% plotting array initialization 
index=1;        %this is a loop counter

T_REC        = zeros(length(t),1);
X_REC        = zeros(length(t),1);
XDOT_REC     = zeros(length(t),1);
XRM_REC      = zeros(length(t),1);
XDOTRM_REC   = zeros(length(t),1);
XERR_REC     = zeros(length(t),1);
XDOTERR_REC  = zeros(length(t),1);
DELTACMD_REC = zeros(length(t),1);%control input
DELTAERR_REC = zeros(length(t),1);
W_REC = zeros(n2,length(t));
VAD_REC = zeros(length(t),1);

%%
for t=t0:dt:tf
    e=x_rm-x;%compute reference model error
   
   
    xref=XREF(index); % assign command
    v_crm=omegan_rm^2*(xref-x_rm(1))-2*zeta_rm*omegan_rm*x_rm(2); %assign reference model output
    v_pd=[Kp Kd]*e; %linear control
    deltaCmd = v_crm + v_pd - v_ad;%Nu, adaptive pseudo control
    delta = deltaCmd;   %approximate inversion, nu=delta
    v_h=deltaCmd-delta; %v_h=0, this would not be zero if there were saturations, we don't consider
                        %control saturations in this sim
          
    %propagate ref. model and plant model
    [x,x_rm,xddot,deltaErr]=inv_pendu(x,x_rm,v_crm,v_h,delta,dt,dt,Wstar);
    x=x+randn(2,1)*sdev;%option to add noise, deault to no noise
    sigma= [sin(pi*x(1)) abs(x(2))*x(2) exp(-x(1)*x(2))]'; %this is the tru uncertainty
  
    
    Wd=-gammaW*e'*P*B*sigma; %online weight update (traditional)
    
    %projection matrix for prioratization, if gamWc is nonzero, then we get
    %prioratization of current training over training on recorded data
    Wc=eye(n2)-Wd*pinv(Wd'*Wd)*Wd'*gamWc;%sigma*sigma'/(sigma'*sigma)*gamWc;%;%*norm(e)%
    
    
    %store data points for concurrent learning
    oldrank=rank(Xbg*Xbg');
    if norm(sigma-Xbg(:,nbg))/norm(sigma)>0.2%&& det(Xbg*Xbg')<2 && index>1 % %norm(v_ad-deltaErr)>0.1 
        nbg = nbg +1;
        Xbg(:,nbg)=sigma;
        Delta_bg(:,nbg)=deltaErr;
         if oldrank==rank(Xbg*Xbg')&& det(Xbg*Xbg')>0.002
             nbg=nbg-1;
         end
        %if rank(Xbg)==3
         %   break
        %end
    end



%do learning on recorded data
Wb=Wd*0;
summer=0;
  for kk=2:nbg
       if norm(Wc*Xbg(:,kk)-(Xbg(:,kk)))<0.1
      Wb=-gammaWbkk*(W'*Xbg(:,kk)-Delta_bg(kk))*Xbg(:,kk)+Wb;
      summer=Xbg(:,kk)*Xbg(:,kk)'+summer;
       end
  end
      wtilde=W-Wstar;

    
     %Lyapunov derivative and lyapunov function
       LYAPDOT(index)=-wtilde'*(sigma*sigma'+Wc*summer)*wtilde;
       LYAP(index)=0.5*e'*P*e+0.5*trace(wtilde'*1/gammaW*wtilde);

    Wdot=Wd+Wc*Wb-norm(e)*W*kappa; %complete training law, note that there is no emod by default
    W=W+dt*Wdot; %integrate learning law (Wdot) using Euler integration
    v_ad=W'*sigma;  %assing adaptive element output
    
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
    Wb_REC(:,index)        = Wb;
    WSTAR_REC(:,index)     =Wstar;
    VAD_REC(index)       =v_ad;
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
plot(T_REC,W_REC(1,:),T_REC,WSTAR_REC(1,:),':',T_REC,W_REC(2,:),T_REC,W_REC(3,:))
hold on
plot(T_REC,WSTAR_REC(2,:),':',T_REC,WSTAR_REC(3,:),':')
legend('W(i)','W^*(i)',0)
xlabel('time (sec)');
ylabel('W')



figure(13)
plot(T_REC,DELTAERR_REC,T_REC,VAD_REC)
grid on
xlabel('Time')
ylabel('\Delta -\nu_{ad}')


figure(14)
plot(T_REC,LYAP)
grid on
xlabel('Time')
ylabel('Lyapunov Candidate')

