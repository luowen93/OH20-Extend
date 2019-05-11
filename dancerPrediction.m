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
yd = -(yd-yd(1))/1000; % Converted to inches 
feedvel = M(:,9)/1000;

%% Save the main variables
close all

% Lane characteristics
d = 0.787;
l = 10.365;
ydp = numderiv(yd)/dt;

% Numerical hypotenuse method
Ld = -ydp.*(l-yd)./sqrt((d-yd).^2+l.^2);
in2m = 0.0254;
plot(t,Ld*in2m);
% xlim([0.7 0.9]);

Ldm = Ld*in2m;

save 'ldm2.mat' Ldm feedbool cutbool clampbool feedvel lanebool xd t
%% Numerical modelling

dt = 0.00146;
vd = numderiv(xd)/dt;

feedvel = feedvel*1000/875*0.3728;

% Construct the toolpoint velocity
tpvel = zeros(size(feedvel));
risingfeed = 1+find(diff(feedbool)==1);
risingclamp = 1+find(diff(clampbool)==1);

for i=1:length(risingfeed)
    tpvel(risingfeed(i):risingclamp(i)) = feedvel(risingfeed(i):risingclamp(i));
end

%%

ke = 8.0;
% Measure the displacement (Initial variable update)
xk1 = 0;
vk1 = 0;
ak1 = 0;

% Define sample time
dt = 0.00146;
n = 17031;
X = zeros(size(tpvel));

extendVel = Ldm.*lanebool;
extendAccel= numderiv(extendVel)/dt;
setnum = 1;
vels = zeros(n,1);
accels = zeros(n,1);

% Backwards difference
% derivCoeff = [1/4, -4/3, 3, -4, 25/12]/dt;
% derivCoeff = -fliplr(derivCoeff);
derivCoeff = [1/12 -2/3 0 2/3 -1/12]/dt;
deriv2Coeff = [-1/12 4/3 -5/2 4/3 -1/12]/dt;

for i=5:n-5

% Update reading
xk = xk1; % Dancer displacement from previous
% vk = vk1;
if( sum(i-20 == risingfeed)>0 )
    vk = 0.3728/2; % Velocity from previous tpvel(i)/2; 
elseif( sum(i-20 == risingclamp)>0)
    vk = -0.3728/2; %tpvel(i-1)/2; 
else 
    vk = vk1;
end
% Measured extend position
e = yd(i);
extendPositions = yd(i-4:i);
extendVelocity = derivCoeff*yd(i-4:i);
extendAcceleration = deriv2Coeff*yd(i-4:i);

% Derivative of the hypotenuse with respect to extension for chain rules
dH = @(e) -(d-e)./(sqrt(l^2+(d-e).^2));
ddH = @(e) (2*d*e-d^2-e.^2+l^2+(d-e).^2)./(l^2+(d-e).^2).^(3/2);

% Derivae of extension with respect to time for chain rules
towExtendVel = in2m*dH(e)*extendVelocity*lanebool(i);
towExtendAccel = in2m*(ddH(e)*extendVelocity^2+dH(e)*extendAcceleration);

towExtendVel = extendVel(i);
towExtendAccel = extendAccel(i);
% Sample it
vels(i) = extendVelocity;
accels(i) = extendAcceleration;
towvels(i) = towExtendVel;
towaccels(i)= towExtendAccel;

% Calculate the acceleration 
% Assumes that we can poll for the feed acceleration and velocity
ak = -ke^2*xk-2*ke*vk+ke*(tpvel(i)+towExtendVel)+1/2*towExtendAccel;
% Predict velocity
vk1 = vk + ak*dt;
% Predict displacement
xk1 = xk + vk*dt;
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
plot(vels)
xlim([0 3000]);
title('Estimated tow speed due to extend/retract');

subplot(3,1,3);
plot(accels);
xlim([0 3000]);
title('Accels bool')

%%
close all
subplot(2,2,1);
plot(vels)
title('RT');
subplot(2,2,3);
plot(extendVel);

subplot(2,2,2);
plot(accels);
title('accels');
subplot(2,2,4);
plot(extendAccel);


