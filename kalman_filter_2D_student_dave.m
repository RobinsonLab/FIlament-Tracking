
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
data_i = importdata([PathName,FileName{k}], '\t', 2);
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
CM_idx = [x_position{filament} y_position{filament}];

% define main variables
dt = 0.2;  %our sampling rate
v_0 = 50; % initial velocity (nm/sec)
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
Q_loc_estimate = []; %  position estimate
vel_estimate = []; % velocity estimate
P_estimate = P;
predic_state = [];
predic_var = [];
r = 5; % r is the radius of the plotting circle
j=0:.01:2*pi; %to make the plotting circle


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
    Q_loc_estimate = [Q_loc_estimate; Q_estimate(1:2)];
    vel_estimate = [vel_estimate; Q_estimate(3:4)];
   
    %% plot the images with the  tracking
    %imagesc(img);
    %axis off
    %colormap(gray);    
    plot(CM_idx(:,2),CM_idx(:,1),'b')
    hold on;
    plot(r*sin(j)+Q_loc_meas(2,t),r*cos(j)+Q_loc_meas(1,t),'.g'); % the actual tracking    
    plot(r*sin(j)+Q_estimate(2),r*cos(j)+Q_estimate(1),'.r'); % the kalman filtered tracking
    %axis([0.5*min(CM_idx(:,1)) 1.5*max(CM_idx(:,1)) 0.5*min(CM_idx(:,2)) 1.5*max(CM_idx(:,2))]);
    axis([1 40000 1 40000])
    hold off
    pause(0.01)
    
end


% Plot velocities
subplot(221),
plot(abs(vel_estimate(1:2:end)))
subplot(223),
plot(abs(vel_estimate(2:2:end)))
subplot(222),
hist(abs(vel_estimate(1:2:end)))
subplot(224),
hist(abs(vel_estimate(1:2:end)))