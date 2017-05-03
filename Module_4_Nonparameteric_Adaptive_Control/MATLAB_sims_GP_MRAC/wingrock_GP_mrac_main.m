%Wingrock Dynamics simulation, with unstructured uncertainty and GP-MRAC
% Code for data presented in
%Chowdhary G., Kingravi, H. Vela P., How J., Nonparameteric adaptive control
%using Gaussian Processes, submitted to IEEE Transactions on Neural Networks and
%Learning Systems
% Authors: Hassan A. Kingravi, Girish Chowdhary, all rights reserved
% this software comes with no guarantees or warranties of any kind
% email girishc@mit.edu for further information
% last modified: 08/10/2012

close all;
clear all;

tic

%% sim params
t0=0;%initial time
tf=50;%final time in seconds
dt=0.005;%simulation time step in seconds
t=t0:dt:tf;

ProgressBarFlag=1;


wn=0.01;% measurement noise covariance

%% Flags for selecting what will be simulated
% if gp_on then GP-MRAC, if gp_on=0, then projection based MRAC

gp_on=1;
postplot=1;
gp_off=not(gp_on);
sprintf('gp on is %d',gp_on)
%% initialization
%x=[1.2 ;1];%[3;,6];%
x = [3;6]; %Initial conditions
%% fixed dictionary when GP is not used
mu = 0.3; %bandwidth of the Gaussian kernel
load rbf_c_store100_minus_1_1
n2=max(size(rbf_c));
rbf_mu = ones(n2+1,1)*mu;
alpha=sdash_rbf(x,rbf_c,rbf_mu,1,n2);
%% baseline control initialization
Kp = 1.2;            % proportional gain
Kd = 1.2;              % derivative gain

A=[0 1;-Kp -Kd];
B = [0; 1];
Q = eye(2);
P = lyap(A',Q);%A' instead of A because need to solve A'P+PA+Q=0
p=0;
[Lp,PP,E] = lqr(A,B,eye(2)*100,1);

%% adaptive control initialization and control parameters

gammaW=10;          %online learning rate %use high learning rate of 5.5
proj_eps=0.001;
theta_max = 7;       %bound on state space norm for projection operator

Wstar=[0.8 0.2314 0.6918 -0.6245 0.0095 0.0214]';% the parameters of the wingrock uncertainty
Wstar_orig=Wstar;
W=zeros(n2+1,1);% NN weights when not using GP-MRAC

v_h=0; %Hedged adaptive control output
v_ad=0;

%% commands
xref=0;
XREF=zeros(length(t),1);
XREF(5/dt:7/dt)=0;
XREF(15/dt:17/dt)=2;
XREF(25/dt:27/dt)=-3;
XREF(35/dt:37/dt)=3;
XREF(45/dt:47/dt)=-2;

%% reference model parameters
omegan_rm = 3;        % reference model natural freq (deg/s)
zeta_rm   = 0.5;        % reference model damping ratio
x_rm = x;


v_crm=omegan_rm^2*(XREF(1)-x_rm(1))-2*zeta_rm*omegan_rm*x_rm(2);
cov_meas=1;
x_input=x;
mean_post=0;
var_post=wn;
deltaErr=Wstar'*[1;x(1);x(2);abs(x(1))*x(2);abs(x(2))*x(2);x(1)^3];
meas=deltaErr;

%% Gaussian Process parameter settings and initialization

bandwidth = 0.3;
noise = wn;
tol = 0.0001;
strg_indx=1;
strg_counter=1;
max_points=100; % max size of GP dictionary
pp=1;
%The following creates the gpr function object that contains the Gaussian
%Process regression functions. A call to the subroutines of gpr is as
%follows: gpr.sub_routine(.), where sub_routine is the name of the
%appropriate subroutine
gpr = onlineGP(bandwidth,noise,max_points,tol);




%% data storage array initialization
index=1;

T_REC        = zeros(length(t),1);
X_REC        = zeros(length(t),1);
XDOT_REC     = zeros(length(t),1);
XRM_REC      = zeros(length(t),1);
XDOTRM_REC   = zeros(length(t),1);
XERR_REC     = zeros(length(t),1);
XDOTERR_REC  = zeros(length(t),1);
DELTACMD_REC = zeros(length(t),1);%control input
DELTAERR_REC = zeros(length(t),1);
VAR_REC      = zeros(length(t),1);
W_REC = zeros(max(size(W)),length(t));
GP_OUT = DELTAERR_REC;
Wbdot_REC=W_REC*0;
STRG_INDX=DELTAERR_REC;
VAD_REC = zeros(length(t),1);

%% Progress Bar

LimPB=(tf-t0)/100/dt;
HI=0;
llll=0;
PBBars='/';
Progressbar='0';


%% start of sim
for t=t0:dt:tf
    
    e=x_rm-x;%compute reference model error
    xref=XREF(index); % bring in the coommands
    
    v_crm=omegan_rm^2*(xref-x_rm(1))-2*zeta_rm*omegan_rm*x_rm(2);% output of the reference model v_rm
    v_pd=[Kp Kd]*e; % the baseline linear control law
    if gp_on==1
        if t == t0
            % if at first time step, initialize GP
            gpr.process(x,meas);
            [mean_post var_post] = gpr.predict(x);
            var_post = 2-var_post; %slight hack only used at initialization
        else
            [mean_post var_post] = gpr.predict(x);      %GP prediction
        end
        gp_out_recurs=mean_post;
    else
        gp_out_recurs=0;
    end
    
    v_ad=W'*alpha*gp_off+gp_out_recurs*gp_on; % adaptive element output, projection or GP-MRAC
    %is selected based on what is
    %on. The interesting point to note is that in this simulation, both the
    %controllers are on but the output of only 1 is chosen. Note also that
    %if GP-MRAC is on, the MRAC adaptation is basically gibberish, since
    %the loop is not closed on MRAC through the tracking error e, it is
    %closed on the GP. But even when MRAC is on, the GP still performs
    %inference, since it trains on Delta and x.
    
    deltaCmd = v_crm + v_pd - v_ad;%Nu: the pseudo-control (desired acceleration)
    delta = deltaCmd;
    v_h=deltaCmd-delta; % no saturation considered here, so v_h=0 (no pseudo-control-hedging)
    
    [x,x_rm,xddot,deltaErr,v_crm]=wingrock_correct(x,x_rm,v_h,delta,dt,dt,Wstar,xref,omegan_rm,zeta_rm);
    x=x+randn(2,1)*wn; % measurement noise is added
    sigma= [1;x(1);x(2);abs(x(1))*x(2);abs(x(2))*x(2);x(1)^3]; % the (unknown) to the controller basis of uncertainty
    
    alpha=sdash_rbf(x,rbf_c,rbf_mu,1,n2); % calls the RBF basis
    deltaErrNoisy=Wstar'*sigma; % the noisey modeling error
    
    Wt = -gammaW*e'*P*B*alpha; % MRAC adaptive law
    Wd = projop(W,Wt,theta_max, proj_eps); %projection based bounds
    
    
    W=W+dt*Wd; %integration of the MRAC weights
    
    
    %% Update GP regression model (perform GP inference)
    if gp_on
        gpr.update(x,deltaErrNoisy);
    end
    
    %record for plotting
    T_REC(index)        = t;
    X_REC(index)        = x(1);
    XDOT_REC(index)     = x(2);
    XRM_REC(index)     = x_rm(1);
    XDOTRM_REC(index)   = x_rm(2);
    XERR_REC(index)     = e(1);
    XDOTERR_REC(index)  = e(2);
    DELTACMD_REC(index) = delta;%control input
    DELTAERR_REC(index) = deltaErrNoisy;
    W_REC(:,index)         =W;
    WSTAR_REC(:,index)     =Wstar;
    VAD_REC(index)       =v_ad;
    GP_OUT(index)       = gp_out_recurs;
    STRG_INDX(index)    = strg_indx;
    VAR_REC(index)      = var_post;
    index = index+1;

%% progress bar    
        cccc=index-(HI*LimPB);
    
    if cccc==LimPB
        HI=HI+1;
        llll=int2str(HI);
        clear Progressbar;
        Progressbar=[];
        for iii=1:HI
            Progressbar=[Progressbar PBBars];
        end
        clc;
        Progressbar=[Progressbar llll '%']
        %MSV_ind=old_eig%min(svd(Xbg*Xbg'))
    end

end

%% plotting
figure(1);
subplot(2,1,1)
plot(T_REC,X_REC, T_REC,XRM_REC,'--');
xlabel('time (seconds)');
ylabel('deg');
title('roll angle');
legend('actual','ref model',0);
grid on;

subplot(2,1,2)
plot(T_REC,XDOT_REC, T_REC,XDOTRM_REC,'--');
xlabel('time (seconds)');
ylabel('xDot (deg/s)');
title('roll rate');
legend('actual','ref model',0);
grid on;

figure(2);
subplot(3,1,1);
plot(T_REC, XERR_REC);
xlabel('time (seconds)');
ylabel('xErr (deg)');
title('Position Error');
grid on;
subplot(3,1,2);
plot(T_REC, XDOTERR_REC);
xlabel('time (seconds)');
ylabel('xDotErr (deg/s)');
title('Angular Rate Error');
grid on;
subplot(3,1,3);
plot(T_REC, DELTACMD_REC);
xlabel('time (seconds)');
ylabel('\delta (deg)');
title('Angular Rate Error');
grid on;

figure(3)
plot(X_REC,XDOT_REC)
grid on
xlabel('roll angle deg')
ylabel('roll rate deg/seconds')

figure(4)
plot(T_REC,DELTAERR_REC,T_REC,VAD_REC,':')
grid on
xlabel('Time')
ylabel('\nu_{ad}')
legend('\Delta','\nu_{ad}')

if gp_on
    %GP estimate of uncertainty with variance
    figure(5)
    f = [VAD_REC+1*sqrt(VAR_REC); flipdim(VAD_REC-1*sqrt(VAR_REC),1)];
    fill([T_REC; flipdim(T_REC,1)], f, [7 7 7]/8)
    hold on;
    plot(T_REC,VAD_REC,'g')
    hold on;
    plot(T_REC,DELTAERR_REC,'b')
    xlabel('Time')
    ylabel('\nu_{ad}')
end

if gp_off
    figure(6)
    plot(T_REC,W_REC)
end

if gp_off
    if postplot==1
        for ii=1:index-1;
            x(1)=X_REC(ii);
            x(2)=XDOT_REC(ii);
            sigma = sdash_rbf(x,rbf_c,rbf_mu,1,n2);
            NU_POST(ii)=W'*sigma;
        end
        
        
        figure(8)
        plot(DELTAERR_REC,'b')
        hold on
        plot(NU_POST,'r')
        title('Uncertainty tracking with learned weights')
        set(findobj(figure(8),'Type','line'),'lineWidth',1.0)
        legend('true uncertainty','estimated uncertainty')
    end
end

toc