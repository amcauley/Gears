w = 2; %toothed gear rotation speed, rad/s. Cut gear is assumed to be 1 rad/s.
r = 6; %toothed gear radius
theta0 = pi; %starting angle, toothed gear

a = 10; %cut gear radius
phi0 = 0; %starting angle, cut gear

d = 15; %separation of gears (center to center)

Nt = 4096; %number of time samples
Np = 3; %number of angular points

phQ = .01; %angle quantization size in radians

t = linspace(-pi*1.5,pi*1.5,Nt)';

%Create time/phase matrix, g. Up/Down within a column moves in time, left/right
%within a row changes the phase of the point we're considering on the toothed gear.
g = linspace(t, t+2*pi, Np+1)(:,1:end-1);
g = pi - theta0 - w*g;

ang = atan(r*sin(g)./(d-r*cos(g))) + phi0 + t*ones(1,Np);
mag = zeros(size(ang));
mag = sqrt((d - r*cos(g)).^2 + r^2*sin(g).^2);
mag(mag > a) = a;

ang = fix(ang/phQ)*phQ;
minAng = min(min(ang));
maxAng = max(max(ang));
ang2 = minAng:phQ:maxAng;
mag2 = zeros(size(ang2));
for n = 1:length(ang2)
    thisAng = ang2(n);
    angColCnt = 0;
    for i = 1:size(ang)(2)
        if(sum(ang(:,i) == thisAng))
            angColCnt = angColCnt + 1;
        end
    end
    if(angColCnt == size(ang)(2))
        mags = mag(ang == thisAng);
        mag2(n) = min(mags);
    end
end

ang2 = ang2(mag2 > 0);
mag2 = mag2(mag2 > 0);

%plotting functions treat each column as a separate function
figure(1);
subplot(2,1,1);
plot(ang, mag);
%hold on
%scatter(ang2, mag2);
%hold off;
axis([minAng maxAng 0 a]);
subplot(2,2,3);
polar(ang, mag);
subplot(2,2,4);
polar(ang2, mag2);
figure(2); %figure(2) used for debugging
subplot(1,1,1);