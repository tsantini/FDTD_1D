% FDTD_1D
% Simulates the wave propagation from a source, passing through a tissue
% and achieving a perfect electrical conductor (PEC). this is a
% demonstration of the finite-difference time-domain algorithm.
%
% ------------------------------- INPUTS: --------------------------------
%   none     
%
% ------------------------------- OPTIONS: -------------------------------
%
% Tissues are selected with the variable 'tissue'
%
% ------------------------------ OUTPUTS: --------------------------------
%   extenal video file in *.avi format
%
% 
%
% ------------------------------ Modified --------------------------------


clear all; close all; fclose('all'); clc; tic

%select one of the 3 tissues for the calculations:
% 'Bladder', 'Fat' or 'Bone Marrow'
tissue = 'Bladder';

%Time and Spatial step
delta_t = 2.5E-11; %time step
delta_z = 1.5E-2; %spatial step

%Initializing Electric and Magnetic Field matrix
E_field = zeros(62.5E-9/delta_t,1000); %Z=1:F, t=0.5:62.5n
H_field = zeros(62.5E-9/delta_t,999); %Z=1.5:F-0.5, t=0:62n

%source constant
T = 1e-9;

%Tissues properties
mu_0 = 1.2566E-6;

if  (exist('tissue')) && (strcmp(tissue,'Fat'))
    episilon_0 = 8.8542E-12; 
    episilon_z  = zeros(1,1000); 
    episilon_z(1:500) = episilon_0;
    episilon_z(501:600) = 5.9215*episilon_0;
    episilon_z(601:1000) = episilon_0;
    sigma_z = zeros(1,1000);
    sigma_z(500:600) = 0.03687;
elseif (exist('tissue')) && (strcmp(tissue,'Bone Marrow'))
    episilon_0 = 8.8542E-12; 
    episilon_z  = zeros(1,1000); 
    episilon_z(1:500) = episilon_0;
    episilon_z(501:600) = 7.2103*episilon_0;
    episilon_z(601:1000) = episilon_0;
    sigma_z = zeros(1,1000);
elseif (exist('tissue')) && (strcmp(tissue,'Bladder'))
    episilon_0 = 8.8542E-12; 
    episilon_z  = zeros(1,1000); 
    episilon_z(1:500) = episilon_0;
    episilon_z(501:600) = 20.093*episilon_0;
    episilon_z(601:1000) = episilon_0;
    sigma_z = zeros(1,1000);
else
    display('Enter with the correct tissue')
    return;
end

%source equation
Ex = @(t) exp(-(t-3*T)^2/T^2);

%Source vector in the firts collum of E_field
i=1;
for t=0:delta_t:6*T+delta_t
    E_field(i,1) = Ex(t); i = i + 1; 
end

%determing the E_field and H_field for every point and time
for i = 1:62.5E-9/delta_t-1 %time related
    for j=1:998 %space related
        E_field(i+1,j+1) = (-(H_field(i,j+1) - H_field(i,j))/delta_z - ...
            sigma_z(j+1)*E_field(i,j+1))*delta_t/episilon_z(j+1) + ...
            E_field(i,j+1);
        H_field(i+1,j) = (-(E_field(i+1,j+1)-E_field(i+1,j))/delta_z)*...
            delta_t/mu_0+H_field(i,j);
    end
        j=999; %last H_field update in the space
        H_field(i+1,j) = (-(E_field(i+1,j+1)-E_field(i+1,j))/delta_z)*...
            delta_t/mu_0+H_field(i,j);
end

%Measuring the time to run
toc

%%
%Ploting the recording
i = 1;
figure, set(gcf, 'Color','white');
bars = -2*ones(1,1000);
bars(501:600) = 2;
nFrames = 20;
vidObj = VideoWriter([tissue '.avi']);
vidObj.Quality = 100;
vidObj.FrameRate = 100;
open(vidObj);

figure(1)
plot(E_field(22.5e-9/delta_t:62.5e-9/delta_t,480));
title('')
xlabel('')
ylabel('')

figure(2)
plot(E_field(:,560));

%ploting the video
space = 1:1000;
for t=1:62.5E-9/delta_t-1
    plot(space,E_field(t,:),space,bars)
    title(['FDTD simulation for ' tissue ', time = ' sprintf('%2.2f',t*delta_t*1E9) 'ns'])
    ylim([-1 1])
    xlabel('Position')
    ylabel('Amplitude')
    writeVideo(vidObj, getframe(gcf));
end
close(gcf)
close(vidObj);
winopen([tissue '.avi'])