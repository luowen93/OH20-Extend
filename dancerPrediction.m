% Servo creel constants

% Read the csv 

clear all

% Servo creel constants
addpath('..\MATLAB Functions');

M = csvread('Ext 2.csv');

dt = 0.00146;
t = M(:,2);
feedbool = M(:,3);
cutbool = M(:,4);
clampbool = M(:,5);
lanebool = M(:,6);
xd = M(:,7)/10000*0.0254;
vd = numderiv(xd)/dt;
ad = numderiv(vd)/dt;

yd = M(:,8);
yd = -(yd-yd(1))/1000;
feedvel = M(:,9)/1000;

%% Save the main variables
close all

% Lane characteristics
d = 0.787;
l = 10.365;
ydp = numderiv(yd)/dt;
Ld = -ydp.*(l-yd)./sqrt((d-yd).^2+l.^2);

in2m = 0.0254;
plot(t,Ld*in2m);
xlim([0.7 0.9]);

Ldm = Ld*in2m;

save 'ldm2.mat' Ldm feedbool cutbool clampbool feedvel lanebool xd t
%% Numerical modelling



dt = 0.00146;
vd = numderiv(xd)/dt;

feedvel = feedvel*1000/875*0.3728;
feedvel = feedvel;

% Construct the toolpoint velocity
tpvel = zeros(size(feedvel));
risingfeed = 1+find(diff(feedbool)==1);
risingclamp = 1+find(diff(clampbool)==1);

for i=1:length(risingfeed)
    tpvel(risingfeed(i):risingclamp(i)) = feedvel(risingfeed(i):risingclamp(i));
end

%%

ke = 8.0;
% feedvel = 1600*0.0254/60;

% Measure the displacement (Initial variable update)
xk1 = 0;
% vk1 = 0.8*feed_vel;
vk1 = tpvel(1);
ak1 = 0;

% Define sample time
dt = 0.00146;
n = 17031;
X = zeros(size(tpvel));

extendVel = Ldm.*lanebool;
extendAccel= numderiv(extendVel)/dt;
setnum = 1;

for i=1:n

% Update reading
xk = xk1; % Dancer displacement from previous
% vk = vk1;
if( sum(i == risingfeed)>0 )
vk = 0.3728/2; % Velocity from previous
elseif( sum(i == risingclamp)>0)
vk = -0.3728/2;
else 
    vk = vk1;
end

% Might want some code to compare the displacement to the model
% displacement

% Calculate the acceleration 
% Assumes that we can poll for the feed acceleration and velocity
ak = -ke^2*xk-2*ke*vk+ke*(tpvel(i)+extendVel(i))+1/2*extendAccel(i);
% ak = -ke^2*xk-2*ke*vk+ke*(feed_vel);

% Predict velocity
vk1 = vk + ak*dt;

% Predict displacement
xk1 = xk + vk*dt;

%
X(i) = xk1;
end

close all


subplot(3,1,1)
plot(X+0.0254)
hold on
plot(xd);
xlim([0 3000]);
ylim([0 0.12]);
legend('Model','Actual');
title('Dancer model')

subplot(3,1,2)
plot(Ldm)
xlim([0 3000]);
title('Estimated tow speed due to extend/retract');

subplot(3,1,3);
plot(clampbool);
xlim([0 3000]);
ylim([-1,2]);
title('Clamp bool')
