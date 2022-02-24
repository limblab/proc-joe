function [t,uc,ui] = runArtifactModel(R, C, dt, T, stimParams)

% set up some constants
R1 = R(1);
R2 = R(2);
R3 = R(3);

C1 = C(1);
C2 = C(2);

numSteps = ceil(T/dt);
t = zeros(numSteps,1);
ui = zeros(numSteps,1);
uc = zeros(numSteps,1);

pw = stimParams(1);
inter = stimParams(2);
I = stimParams(3);
timeStart = 3.48; % fixed 
timeSwitch = timeStart + pw;
timeEnd = timeStart + 2*pw+inter;

ui(ceil(timeStart/dt):ceil(timeSwitch/dt)) = I;
ui(ceil((timeSwitch+inter)/dt):ceil(timeEnd/dt)) = -I;

for i = 2:numSteps
    t(i) = t(i-1) + dt;
    uc(i) = uc(i-1) + dt*(R1*ui(i-1)-(1/R3)*(R1+R2+R3)*uc(i-1))/(R1*C1 + R2*C1);
end

end