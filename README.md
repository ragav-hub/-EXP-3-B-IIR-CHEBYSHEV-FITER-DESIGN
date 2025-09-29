# EXP 3 : IIR-CHEBYSHEV-FITER-DESIGN

## AIM: 

 To design an IIR Chebyshev filter  using SCILAB. 

## APPARATUS REQUIRED: 
PC installed with SCILAB. 

## PROGRAM (LPF): 

clc ;
close ;
wp=input('Enter the pass band frequency (Radians )= ' );
ws=input('Enter the stop band frequency (Radians )= ' );
alphap=input( ' Enter the pass band attenuation (dB)=' );
alphas=input( ' Enter the stop band attenuation(dB)=' );
T=input('Enter the Value of sampling Time=');

omegap=(2/T)*tan(wp/2);
disp(omegap,'omegap=');
omegas=(2/T)*tan(ws/2);
disp(omegas,'omegas=');
N=acosh(sqrt(((10^(0.1*alphas))-1)/((10^(0.1*alphap))-1)))/(acosh(omegas/omegap));
disp(N,'N=');
N=ceil(N);
disp(N,'Round off value of N=');
omegac=omegap/(((10^(0.1*alphap)) -1)^(1/(2* N)));
disp(omegac,'omegac=');
Epsilon = sqrt ((10^(0.1*alphap))-1);
disp(Epsilon,'Epsilon=');
[pols ,gn] = zpch1(N, Epsilon,omegap );
disp(gn,'Gain');
disp(pols,'Poles');
hs=poly(gn,'s','coeff')/real(poly(pols,'s'));
disp(hs,'Analog Low pass Chebyshev Filter Transfer function');
z=poly(0,'z');
Hz=horner(hs,(2/ T)*((z -1)/(z+1)))
disp(Hz,'Digital LPF Transfer function H(Z)=');
HW=frmag(Hz,512); 
w=0:%pi/511:%pi ;
plot(w/%pi,abs(HW));
xlabel(' Normalized Digital Frequency w');
ylabel('Magnitude ');
title(' Frequency Response of Chebyshev IIR LPF');


## PROGRAM (HPF): 

clc;
close;

wp = input('Enter the pass band frequency (Radians )= ');  
ws = input('Enter the stop band frequency (Radians )= ');
alphap = input('Enter the pass band attenuation (dB)= ');
alphas = input('Enter the stop band attenuation (dB)= ');
T = input('Enter the Value of sampling Time= ');

// ============================
// Pre-warping (Bilinear Transformation)
// ============================
omegap = (2/T)*tan(wp/2);   // Passband edge (analog)
disp(omegap,'omegap=');
omegas = (2/T)*tan(ws/2);   // Stopband edge (analog)
disp(omegas,'omegas=');

// ============================
// Order of the HPF
// ============================
// For HPF: use acosh(omegap/omegas)
N = acosh(sqrt(((10^(0.1*alphas))-1)/((10^(0.1*alphap))-1))) / acosh(omegap/omegas);
disp(N,'N=');
N = ceil(N);
disp(N,'Round off value of N=');

// ============================
// Cutoff frequency
// ============================
omegac = omegap / (((10^(0.1*alphap))-1)^(1/(2*N)));
disp(omegac,'omegac=');

// ============================
// Chebyshev Prototype (Normalized LPF)
// ============================
Epsilon = sqrt((10^(0.1*alphap))-1);
disp(Epsilon,'Epsilon=');

[pols, gn] = zpch1(N, Epsilon, 1);   // normalized LPF prototype at 1 rad/s
disp(gn,'Gain');
disp(pols,'Poles');

s = poly(0,'s');   // Laplace variable
hs = poly(gn,'s','coeff') / real(poly(pols,'s'));
disp(hs,'Analog Normalized Chebyshev LPF');

// ============================
// LPF → HPF Transformation (s -> omegac/s)
// ============================
sh = horner(hs, omegac/s);
disp(sh,'Analog Chebyshev High-Pass Filter');

// ============================
// Bilinear Transformation: s -> (2/T)*((z-1)/(z+1))
// ============================
z = poly(0,'z');   // z-domain variable
Hz = horner(sh, (2/T) * ((z-1)/(z+1)));
disp(Hz,'Digital HPF Transfer function H(Z)=');

// ============================
// Frequency Response
// ============================
HW = frmag(Hz, 512);
w = 0:%pi/511:%pi;

plot(w/%pi, abs(HW));
xlabel('Normalized Digital Frequency ×π rad/sample');
ylabel('Magnitude');
title('Frequency Response of Chebyshev IIR High-Pass Filter');


## OUTPUT (LPF) : 

<img width="1919" height="1199" alt="Screenshot 2025-09-29 152821" src="https://github.com/user-attachments/assets/72596626-14f8-472e-bcd7-c3c64e9cdf17" />

<img width="1919" height="1199" alt="Screenshot 2025-09-29 152834" src="https://github.com/user-attachments/assets/555bff56-5ec8-4bd3-bd8d-33ae3e792cdd" />

<img width="1919" height="1199" alt="Screenshot 2025-09-29 152855" src="https://github.com/user-attachments/assets/0529d43e-4895-4126-ade7-6b9abb31e89e" />

## OUTPUT (HPF) : 

<img width="1914" height="1199" alt="Screenshot 2025-09-29 153716" src="https://github.com/user-attachments/assets/0374a475-a885-4c6f-9d71-da703346677a" />

<img width="1919" height="1199" alt="Screenshot 2025-09-29 153735" src="https://github.com/user-attachments/assets/6d250437-8aa6-4489-b2bd-6bb3638831e6" />

<img width="1919" height="1199" alt="Screenshot 2025-09-29 153805" src="https://github.com/user-attachments/assets/c4e45460-36ca-46e4-99be-e86fee6a2076" />


## RESULT: 
Thus, the IIR Chebysheve filter (LPF & HPF) using SCILAB was successfully executed and the output is verified.
