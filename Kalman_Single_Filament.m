
%Student Dave's tutorial on:  Object tracking in image using 2-D kalman filter
%HEXBUG TRACKING!
%Copyright Student Dave's Tutorials 2012
%if you would like to use this code, please feel free, just remember to
%reference and tell your friends! :)

%here we take the hexbug extracted coordinates and apply the kalman fiter
%to the 2-dimensions of motion  within the image
%to see if it can do better than the tracking alone
%requires matlabs image processing toolbox.

% Edited by  Christpher Solis-Ocampo
% South Dakota State University
% July 23 2014

clear all;
close all;
clc;

[FileName,PathName] = uigetfile({'*.txt','Text(*.txt)'},...
    'Select Filament Coordinates file','MultiSelect','on'); 
for k = 1:length(FileName);
data_i = importdata([PathName,FileName{k}], '\t', 3);
Data(k) = data_i;
end


time = cell(length(FileName),1);
x_position = cell(length(FileName),1);
y_position = cell(length(FileName),1);

for k = 1:length(FileName);
time_i = [];
time_i = Data(k).data(:,2);
time{k} = time_i; 
        
x_pos_i = [];
x_pos_i = Data(k).data(:,3);
x_position{k} = x_pos_i; 

y_pos_i = [];
y_pos_i = Data(k).data(:,4);
y_position{k} = y_pos_i; 
end

%%
filament = 1;
CM_idx2 = [x_position{filament} y_position{filament}];
CM_idx = CM_idx2./1000;
timeF = time{filament};

% define main variables
dt = 0.2;  %our sampling rate
v_0 = 0.15; % initial velocity (nm/sec)
S_frame = 1; %starting frame
u = .005; % define acceleration magnitude
Q= [CM_idx(S_frame,1); CM_idx(S_frame,2); v_0; v_0]; %initized state--it has four components: [positionX; positionY; velocityX; velocityY] of the hexbug
Q_estimate = Q;  %estimate of initial location estimation of where the hexbug is (what we are updating)
HexAccel_noise_mag = 0.5; %process noise: the variability in how fast the Hexbug is speeding up (stdv of acceleration: meters/sec^2)
tkn_x = 1;  %measurement noise in the horizontal direction (x axis).
tkn_y = 1;  %measurement noise in the horizontal direction (y axis).
Ez = [tkn_x 0; 0 tkn_y];
Ex = [dt^4/4 0 dt^3/2 0; ...
    0 dt^4/4 0 dt^3/2; ...
    dt^3/2 0 dt^2 0; ...
    0 dt^3/2 0 dt^2].*HexAccel_noise_mag^2; % Ex convert the process noise (stdv) into covariance matrix
P = Ex; % estimate of initial Hexbug position variance (covariance matrix)

% Define update equations in 2-D! (Coefficent matrices): A physics based model for where we expect the HEXBUG to be [state transition (state + velocity)] + [input control (acceleration)]
A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]; %state update matrice
B = [(dt^2/2); (dt^2/2); dt; dt];
C = [1 0 0 0; 0 1 0 0];  %this is our measurement function C, that we apply to the state estimate Q to get our expect next/new measurement


% initize result variables
% Initialize for speed
Q_loc = []; % ACTUAL hexbug motion path
vel = []; % ACTUAL hexbug velocity
Q_loc_meas = []; % the hexbug path extracted by the tracking algo

% initize estimation variables
Qx_loc_estimate = []; %  X-position estimate
Qy_loc_estimate = []; %  Y-position estimate
velx_estimate = []; % X-velocity estimate
vely_estimate = []; % Y-velocity estimate
P_estimate = P;
predic_state = [];
predic_var = [];

% to save the animation as a .avi file
aviobj = avifile([PathName,FileName{k}(1:end-4),'.avi'],'compression','None');
fig = figure;
for t = 1:length(CM_idx)
    
    % load the image
    %img_tmp = double(imread(f_list(t).name));
    %img = img_tmp(:,:,1);
    % load the given tracking
    Q_loc_meas(:,t) = [ CM_idx(t,1); CM_idx(t,2)];
    
    %% do the kalman filter   
    
    % Predict next state of the Hexbug with the last state and predicted motion.
    Q_estimate = A * Q_estimate + B * u;
    predic_state = [predic_state; Q_estimate(1)] ;
    %predict next covariance
    P = A * P * A' + Ex;
    predic_var = [predic_var; P] ;
    % predicted Ninja measurement covariance
    % Kalman Gain
    K = P*C'*inv(C*P*C'+Ez);
    % Update the state estimate.
    if ~isnan(Q_loc_meas(:,t))
        Q_estimate = Q_estimate + K * (Q_loc_meas(:,t) - C * Q_estimate);
    end
    % update covariance estimation.
    P =  (eye(4)-K*C)*P;
    
    %% Store data
    Qx_loc_estimate = [Qx_loc_estimate; Q_estimate(2)];
    Qy_loc_estimate = [Qy_loc_estimate; Q_estimate(1)];
    velx_estimate = [velx_estimate; Q_estimate(4)];
    vely_estimate = [vely_estimate; Q_estimate(3)];
    Velxy_estimate = (velx_estimate.^2+vely_estimate.^2).^0.5;
    %% plot the images with the  tracking
    %imagesc(img);
    %axis off
    %colormap(gray);    
    
    % take this out cause with Qx_loc_estimate you know all these positions
    % and make a looop for all the filaemnts!
    subplot(211),
    plot(CM_idx(:,2),CM_idx(:,1),'b')
    set(fig,'Color',[1 1 1]);
    hold on;
    plot(Q_loc_meas(2,t),Q_loc_meas(1,t),'*k','MarkerSize',10); % the actual tracking    
    plot(Q_estimate(2),Q_estimate(1),'or','MarkerSize',10); % the kalman filtered tracking
    axis([0.98*min(CM_idx(:,2)) 1.02*max(CM_idx(:,2)) 0.98*min(CM_idx(:,1)) 1.02*max(CM_idx(:,1))]);
    xlabel('x-coordiante (\mum)','FontSize',18,'FontName','Helvetica');
    ylabel('y-coordiante (\mum)','FontSize',18,'FontName','Helvetica');
    hold off
    %axis([1500 2500 1000 1400]) % molecule 2
    %axis([1 40000 1 40000])
    subplot(212),
    plot(timeF(1:t),Velxy_estimate(1:t));
    axis([1 max(timeF) 0 0.7]);
    xlabel('time (sec)','FontSize',18,'FontName','Helvetica');
    ylabel('speed (\mum/s)','FontSize',18,'FontName','Helvetica');
    hold off
    pause(0.01)
     F = getframe(fig);
    aviobj = addframe(aviobj,F);
    
    %time{filament}(1):dt:time{filament}(t)
end

close(fig);
aviobj = close(aviobj);

Vel_estimate = ((velx_estimate).^2+(vely_estimate).^2).^0.5;

%% Plot velocities
%figure,
%subplot(211),
%plot(abs(Vel_estimate))
%subplot(212),

%% PLot Histogram
[Counts,Bins] = hist(Vel_estimate,30);
X_Data = Bins';
Y_Data = Counts';
X_Data2 = min(X_Data):0.1*dt:max(X_Data);

CurvFit = fit(X_Data,Y_Data,'gauss1');

fig2 = figure; 
bar(X_Data,Y_Data,'hist');
hold on 
plot(X_Data2,CurvFit(X_Data2),'-','Color',[1,0,0]);

title([' Velocity ',FileName{1}],'FontSize',18,'FontName','Helvetica');
xlabel('Velocity (\mum/s)','FontSize',18,'FontName','Helvetica');
ylabel('Counts','FontSize',18);
text(0.96*min(X_Data),1.1*max(Y_Data),['< v > = ',num2str(CurvFit.b1,'%0.2f'),'\pm',num2str(CurvFit.c1,'%0.2f'), ' \mum/s'],'FontSize',18,'FontName','Helvetica');
xlim([0 1.2*max(X_Data)]);
ylim([0 1.5*max(Y_Data)]);
hold off
saveas(fig2,[PathName,FileName{k}(1:end-4),'.eps']) 

