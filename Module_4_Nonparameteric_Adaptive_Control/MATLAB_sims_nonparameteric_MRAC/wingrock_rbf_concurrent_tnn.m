
%Wingrock Dynamics simulation, with unstructured uncertainty and BKR-CL
% Code for data presented in
%Kingravi H., Chowdhary G., Vela P., Johnson E., A Reproducing Kernel 
%Hilbert Space Approach for the Online Update of Radial Bases in 
%Neuro-Adaptive Control,   IEEE Transactions on Neural Networks and 
%Learning Systems, Vol 23, no. 7, pp 1130-1141, July 2012.
% Authors: Hassan A. Kingravi, Girish Chowdhary, all rights reserved
% this software comes with no guarantees or warranties of any kind
% email girishc@mit.edu for further information
% last modified: 08/10/2012

close all;
clear all;
clc

ProgressBarFlag=1;

%% simulation params
t0 = 0;              %intial time
tf = 60;             %170%final time
dt = 0.01;           %time step
t = t0:dt:tf;        %time interval


%% initialization

x=[2 ;2];            %intial state
sdev=0.0;            %noise std deviation if used
Wstar=[0.8 0.2314 0.6918 -0.6245 0.0095 0.0214]'; %ideal weights
W=Wstar*0;           %weight initialization

Kp = 1.2;            % proportional gain
Kd = 1.2;

A=[0 1;-Kp -Kd];     %Linear controller
B = [0; 1];          %B matrix
Q = eye(2);          %for Lyapunov equation
P = lyap(A',Q);      %A' instead of A because need to solve A'P+PA+Q=0
[Lp,PP,E] = lqr(A,B,eye(2)*100,1);  % Lyap eqn

%% adaptive control initialization
gammaW=10;          %online learning rate %use high learning rate of 5.5
gammaWbkk=0.5;       %learning rate for recorded data updates, set to zero
gamWc=0;             %gamWc=0 means Wc=I (default). Set to 1 to prioritize
                     %training on current data over recorded
kappa = 0.01;        %sigmamod/emod gain (default is off)
theta_max = 2;       %bound on state space norm for projection operator 
proj_eps = 0.001;     %projection operator tolerance 

%% Flags for what is going to be simulated
%you can simulate 4 possibilities:

%1. concurrent learning with skernel (BKR-CL): for this set use_sigma_mod=0;
%use_e_mod=0;, and use_skernel=1
%2. concurrent learning without skernel (CL): for this set use_sigma_mod=0;
% use_e_mod=0;, and use_skernel=0
% 3. sigma mod, for this set use_sigma_mod=1 and set use_e_mod=0;, and use_skernel=0
% 4. e mod: for this set use_e_mod=1 and set use_sigma_mod=0;, and use_skernel=0
% 5. 3,4 can also be simulated with skernel (BKR), for this use_skernel=1
%default is BKR-CL

%flags
use_sigma_mod = 0;   %if 1 then NO CL
use_e_mod = 0;       %if 1 NO CL
use_proj_op = 0;     %if 1 NO CL
use_skernel = 1;     %flag for using skernel, can be used with any method
postplot = 1;  % when 1, performs post-simulation analysis
save_figs = 0; % when 1, saves figures as pdfs
%%
%budget sizes
dict_size = 12; % Any dictionary size can be picked, however if you wanted
                % to compare with fixed RBFs, then fixed dictionaries of
                % upto 14 RBF  centers are provided
max_nbg = 25;
old_eig = 0;    %will get overwritten after first point is added
add_new_points = 1;

%% concurrent learning parameter intialization
nbg = 0;
X_state = zeros(2);  %this stores the states for concurrent learning
Delta_bg = 0;


last_point=1;

%% commands are set in this section
xref=0;
XREF=zeros(length(t),1);
XREF(5/dt:10/dt)=1;
XREF(20/dt:25/dt)=1;
XREF(35/dt:40/dt)=1;
XREF(50/dt:55/dt)=1*0;

XREF(50/dt:60/dt)=1;
XREF(65/dt:70/dt)=1;
XREF(85/dt:90/dt)=1;
XREF(105/dt:110/dt)=1;
XREF(125/dt:130/dt)=1;
XREF(145/dt:150/dt)=1;
XREF(165/dt:170/dt)=1;
XREF(185/dt:190/dt)=1;
XREF(205/dt:210/dt)=1;
XREF(220/dt:225/dt)=1;
XREF(235/dt:240/dt)=1;
XREF(250/dt:255/dt)=1;
XREF(265/dt:270/dt)=1;
XREF(285/dt:290/dt)=1;

%% reference model parameters
omegan_rm = 1;       %reference model natural freq (deg/s)
zeta_rm   = 0.5;     %reference model damping ratio
x_rm = x;            %initialize ref. model state at x(0)
v_h=0;               %initialize pseudo control
v_ad=0;

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
DELTACMD_REC = zeros(length(t),1); %control input
DELTAERR_REC = zeros(length(t),1);
NU_AD_REC    = zeros(length(t),1);
VAD_REC      = zeros(length(t),1);
SVD_REC      = zeros(length(t),1);

%%
phi_space=[-3,3];
phid_space=[-4,4];



%% Progress Bar

LimPB=(tf-t0)/100/dt;
HI=0;
llll=0;
PBBars='/';
Progressbar='0';

%% s kernel initialization


if use_skernel==0
    %load set of saved centers
    rbf_c = [];
    if dict_size==4
        load rbf_c_store4
    elseif dict_size==8
        load rbf_c_store8
    elseif dict_size==12
        load rbf_c_store12
     elseif dict_size==14
        load rbf_c_store14
    end

    n2 =max(size(rbf_c));
    rbf_c(:,1)=[0;0]; %sets the first RBF to be at the origin
    %rbf_mu = ones(n2+1,1)*mu;
    %W = zeros(n2+1,1);
    %W_REC = zeros(n2,length(t));
end
%this section initializes the sparse kernel object

tic
mu = 2;
sparse_tol = 0.0005;

if use_skernel==1
    dictionary_size = dict_size;
    rbf_c = zeros(2,dictionary_size);    
    %The following line creates the sparse kernel function object
    %subroutines of the object can be called using the . opertor
    %e.g. the get_dictionary subroutine that returns the current set of RBF
    %dictionary is called as
    %rbf_c=o.get_dictionary;
    o = skernel(rbf_c',0.01,dictionary_size,sparse_tol);
    rbf_c = o.get_dictionary();% get current dictionary of RBF centers
    rbf_c = rbf_c';%transpose, needed for consistency
    n2 = size(rbf_c,2);%get n2 based on center allocation
end


rbf_mu = ones(n2+1,1)*mu; % the RBF bandwidth
Xbg=zeros(n2+1);%the same size as sigma

W = zeros(n2+1,1); %initialize the weights
W_REC = zeros(n2,length(t)); %the weight recording array

%initial plot of centers
figure(6)
plot(rbf_c(1,:),rbf_c(2,:),'+r')
hold on
title('Plot of RBF centers')
set(findobj(figure(6),'Type','line'),'lineWidth',1.5)

if use_sigma_mod==1 || use_e_mod==1 || use_proj_op ==1
    Wc = eye(n2+1)*0; %no concurrent learning if sigma mod is on
else
    Wc = eye(n2+1);
end
Wb=W*0;

%% Actual start of sim
tic
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
    %[x,x_rm,xddot,deltaErr]=inv_pendu(x,x_rm,v_crm,v_h,delta,dt,dt,Wstar);
    [x,x_rm,xddot,deltaErr,v_crm]=wingrock_correct(x,x_rm,v_h,delta,dt,dt,Wstar,xref,omegan_rm,zeta_rm);
    x=x+randn(2,1)*sdev;%option to add noise, deault to no noise
    
    % save current variables
    rbf_c_old=rbf_c;
    XbgOld=Xbg;
    Delta_bg_old=Delta_bg;
    X_state_old=X_state;
    n2_old=n2;
    last_point_old=last_point;
    nbg_old=nbg;
    
    if use_skernel==1
        if mod(t*10,5) == 0
            %these two elements are shifted downward
            %send current state to add_point, decided whether or not to add the
            %current state into the dictionary as rbf center. If flag then added
            %and kicked out ptNr, flag = 0 means no center added
            [flag, ptNr, old_x] = o.add_point(x');
            
            if(flag > 0)
                %create the new dictionary
                old_x = old_x';
                rbf_c = o.get_dictionary();
                rbf_c = rbf_c';
                n2 = size(rbf_c,2);
                rbf_mu = mu*ones(n2+1,1);
                new_center=1;
                
                %you now need to consider adding the center to the Xbg matrix; 
                %first, via the X_state matrix
                if flag == 1
                    %nothing to check; just add the new state to the matrix
                    %if (nbg)<max_nbg
                    nbg = nbg+1;%size(rbf_c,2);
                    end_point = size(X_state,2);
                    X_state(:,nbg) = x;
                    Delta_bg(:,nbg)=deltaErr;
                    Timebg(nbg)=index;
                    last_point=nbg;
                    %end
                    %
                elseif flag == 2
                    %find old center and replace
                    ind = find(ismember(X_state',old_x','rows') == 1);
                    if(ind)%why am I getting twice the same meas? (fixed it)
                        if(size(ind,1) > 1)
                            ind %just a debugging tool, this should not happen
                        else
                            X_state(:,ind) = x;
                            Delta_bg(:,ind(1))=deltaErr;
                            Timebg(ind(1))=index;
                            last_point=ind;
                        end
                    end
                end
                
                Xbg = zeros(n2+1,nbg); %repopulate Xbg
                for i=1:nbg
                    Xbg(:,i) = sdash_rbf(X_state(:,i),rbf_c,rbf_mu,1,n2);
                end
                if ptNr > 0 && new_center==1
                    %add a state to W, or replace the old one
                    %reset the weights for the points that were kicked out
                    % W(ptNr) = 0; %this is not recommended
                elseif ptNr==0 %or add a weight for the point that has been added
                    W = [W; 0];
                end
            end
        end %end mod if
    end %if use skernel
    sigma = sdash_rbf(x,rbf_c,rbf_mu,1,n2);
    Wd=-gammaW*e'*P*B*sigma; %online weight update (traditional)
    
    
    %store data points for concurrent learning
    flag_point=0;
    %store data points
    old_eig=(min(svd(Xbg*Xbg',0)));%sqrt(eig(Xbg(:,1:nbg)*Xbg(:,1:nbg)')));
    
    if add_new_points==1 && use_sigma_mod==0 && use_e_mod==0 && use_proj_op == 0
        if nbg==0
            nbg = nbg +1;
            flag_point=1;
            Xbg(:,nbg)=sigma;
            X_state(:,nbg) = x; %we're going to assume that these remain relevant for all time
            Timebg(nbg)=t;
            Delta_bg(:,nbg)=deltaErr;
            last_point=nbg;
        elseif nbg>0
            if norm(Xbg(:,last_point)-sigma)/norm(sigma)>=0.0001 || rank([Xbg sigma],1e-5)>nbg% ||  norm(min(svd([Xbg sigma])))>0.1%norm(v_ad-deltaErr)>0.051 %&& det(Xbg*Xbg')<1%norm(x-Xbg(:,nbg))>0.2
                if nbg<max_nbg
                    nbg = nbg +1;
                    flag_point=1;
                    last_point=nbg;
                    XbgOld=Xbg;
                    Xbg(:,nbg)=sigma;
                    X_state(:,nbg) = x; %we're going to assume that these remain relevant for all time
                    Delta_bg(:,nbg)=deltaErr;
                    Timebg(nbg)=index;
                    %new_eig=svd(Xbg');
                    %A=Xbg;
                else
                    XbgOld=Xbg;
                    [Xbg,Delta_bg,nbg,I]=data_point_remover(XbgOld,Delta_bg,sigma,old_eig,deltaErr,nbg);
                    %if min(svd(Xbg*Xbg'))<old_eig
                    %    disp('min')
                    %end
                    if norm(Xbg(:,1:nbg)-XbgOld(:,1:nbg))>0
                        X_state(:,I) = x; %we're going to assume that these remain relevant for all time
                        Timebg(I)=index;
                        flag_point=1;
                        last_point=I;
                    else
                        flag_point=0;
                    end
                    
                end
            end
        end
    end
    
    %rank(Xbg);
    %do learning on recorded data
    if use_sigma_mod == 0 && use_e_mod == 0 && use_proj_op == 0 %do only if conc is on
        Wb=Wd*0;
        summer=0;
        for kk=1:nbg
            if norm(Xbg(:,kk)-(Xbg(:,kk)))<1
                Wb=-gammaWbkk*(W'*Xbg(:,kk)-Delta_bg(kk))*Xbg(:,kk)+Wb;
                summer=Xbg(:,kk)*Xbg(:,kk)'+summer;
            end
        end
    end
    
    
    if use_sigma_mod == 1
        Wdot=Wd-W*kappa; %complete training law, note that there is no emod by default    
    elseif use_e_mod == 1
        Wdot=Wd-W*kappa*norm(e); %complete training law, note that there is no emod by default
    elseif use_proj_op == 1
        %use_proj_op
        Wdot= projop(W,Wd,theta_max, proj_eps);
    else %use concurrent learning with or without skernel        
         Wdot=Wd+Wb;
    end
    W=W+dt*Wdot; %integrate learning law (Wdot) using Euler integration
    v_ad=W'*sigma;  %assign adaptive element output
    
    
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
    W_REC(1:n2+1,index)         =W;
    VAD_REC(index)       =v_ad;
    SVD_REC(index)       =old_eig;
    SVD_SUMMER(index)    =min(svd(summer));
    index = index+1;
    
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
toc



%% plotting
figure(1);
axis([0 t -0.6 1.5]);
subplot(2,1,1)
plot(T_REC,X_REC, T_REC,XRM_REC,'--');
xlabel('time (sec)');
ylabel('pi-rad');
title('roll angle');
legend('actual','ref model',0);
grid on;
%
axis([0 t -0.6 1.5]);
subplot(2,1,2)
axis([0 t -1 1]);
plot(T_REC,XDOT_REC, T_REC,XDOTRM_REC,'--');
xlabel('time (sec)');
ylabel('xDot (pi-rad/s)');
title('roll rate');
legend('actual','ref model',0);
grid on;
%
figure(2);
axis([0 t -1 1]);
subplot(2,1,1);
plot(T_REC, XERR_REC);
xlabel('time (sec)');
ylabel('xErr (pi-rad)');
title('Position Error');
grid on;
axis([0 t -1 1]);
subplot(2,1,2);
plot(T_REC, XDOTERR_REC);
xlabel('time (sec)');
ylabel('xDotErr (pi-rad/s)');
title('Angular Rate Error');
grid on;


figure(4)
plot(T_REC,W_REC)
hold on
legend('W(i)','W^*(i)',0)
xlabel('time (sec)');
ylabel('W')
title('Convergence of weights')


figure(5)
%grid on
hold on
%for i=2:length(X_REC)
plot(DELTAERR_REC,'b')
plot(VAD_REC,'r')
hold off;
set(findobj(figure(5),'Type','line'),'lineWidth',1.5)
xlabel('time (sec)');
title('Online uncertainty tracking')
legend('true uncertainty','estimated uncertainty')


figure(6)
plot(rbf_c(1,:),rbf_c(2,:),'db')
plot(X_REC,XDOT_REC,'--g')
%hold off;
set(findobj(figure(6),'Type','line'),'lineWidth',1.5)
title('RBF centers in state space')
legend('old centers','new centers')



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

%%save some figures

if save_figs == 1
 set(figure(2),'Position',[100 100 800 600]);
 set(figure(2),'PaperOrientation','portrait','PaperSize',[8.5 6.0],'PaperPositionMode', 'auto', 'PaperType','<custom>');
 saveas(figure(2),'tracking_error','pdf')

 set(figure(6),'Position',[100 100 800 600]);
 set(figure(6),'PaperOrientation','portrait','PaperSize',[8.5 6.0],'PaperPositionMode', 'auto', 'PaperType','<custom>');
 saveas(figure(6),'state_space','pdf')

 set(figure(5),'Position',[100 100 800 600]);
 set(figure(5),'PaperOrientation','portrait','PaperSize',[8.5 6.0],'PaperPositionMode', 'auto', 'PaperType','<custom>');
 saveas(figure(5),'u1','pdf')

 set(figure(8),'Position',[100 100 800 600]);
 set(figure(8),'PaperOrientation','portrait','PaperSize',[8.5 6.0],'PaperPositionMode', 'auto', 'PaperType','<custom>');
 saveas(figure(8),'post_plot_uncertainty_BKR-CL','pdf')
end

