w = 1; %toothed gear rotation speed, rad/s. Cut gear is assumed to be 1 rad/s.
rOut = 3; %toothed gear outer radius
rIn = rOut*.76; %toothed gear inner radius
theta0 = pi; %starting angle, toothed gear

a = 3; %cut gear radius
phi0 = 0; %starting angle, cut gear

d = a+rIn; %separation of gears (center to center)

Nt = 128; %number of time samples
Np = 18; %number of angular points per replication unit of toothed gear

t = linspace(-pi*1.5,pi*1.5,Nt)';

%INPUT GEAR DEFINITION (FOR A SINGLE TOOTH)
nTeeth = 10;
maxAngIn = 2*pi/nTeeth;
angIn = linspace(0, maxAngIn, Np);
radIn = ones(1,Np)*rIn;
idx = angIn <= maxAngIn*1/6;
radIn(idx) = rIn + (rOut-rIn)*angIn(idx)/(maxAngIn/6);
idx = (angIn > maxAngIn*1/6) & (angIn <= maxAngIn*3/6);
radIn(idx) = rOut;
idx = (angIn > maxAngIn*3/6) & (angIn <= maxAngIn*4/6);
radIn(idx) = rIn + (rOut-rIn)*(maxAngIn*4/6-angIn(idx))/(maxAngIn/6);

%plot(radIn)
%return

%replicate input profile around the entire 2*pi radian gear
angIn = linspace(0, 2*pi, Np*nTeeth);
radIn = repmat(radIn, 1, nTeeth);


%Create time/phase matrix, g. Up/Down within a column moves in time, left/right
%within a row changes the phase of the point we're considering on the toothed gear.
g = linspace(t, t+2*pi, Np*nTeeth);
g = pi - theta0 - w*g;

r = ones(Nt,1)*radIn;
ang = atan(r.*sin(g)./(d-r.*cos(g))) + phi0 + t*ones(1,Np*nTeeth);
mag = zeros(size(ang));
mag = sqrt((d - r.*cos(g)).^2 + r.^2.*sin(g).^2);
mag(mag > a) = a;

%plotting functions treat each column as a separate function
figure(1);
subplot(1,2,1);
polar(angIn, radIn, 'k');
subplot(1,2,2);
polar(ang, mag, 'k');
figure(2); %figure(2) used for debugging
subplot(1,1,1);