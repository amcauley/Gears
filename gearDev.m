w = 1; %toothed gear rotation speed, rad/s
r = 6; %toothed gear radius
theta0 = pi; %starting angular offset

rho = 1; %cut gear rotation speed, rad/s
a = 10; %cut gear radius
phi0 = 0; %starting angular offset

d = 15; %separation of gears (center to center)

Nt = 512; %number of time samples

rotSpeed = min(w,rho);
if (rotSpeed == 0)
    rotSpeed = max(w,rho);
    if (rotSpeed == 0)
        display('Error: at least one of the gears must have non-zero rotation speed.')
    endif
endif

t = linspace(-pi/rotSpeed,pi/rotSpeed,Nt);

gamma = pi-theta0-w*t;

ang = atan(r*sin(gamma)./(d-r*cos(gamma))) + phi0 - rho*t;
mag = sqrt((d - r*cos(gamma)).^2 + r^2*sin(gamma).^2);
mag(mag > a) = a;

subplot(2,1,1);
plot(ang, mag);
axis([t(1) t(end) 0 a]);
xlabel('angle (rad)');
ylabel('radius');
subplot(2,2,3);
polar(ang, mag);
subplot(2,2,4);
polar(ang, mag);