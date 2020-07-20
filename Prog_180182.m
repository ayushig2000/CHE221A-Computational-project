clc; clear all; close all
a =  7.356;
b = 0.06425;
R = 0.0821;
Tc = (8*a)/(27*R*b);
Pc = a/(27*b*b);
sp=[63.33 63.91 64.57 65.22 65.87]; %using antoine for vanderwaals
% Vcoef = [Pc -(Pc*b)-(R*Tc) a -a*b];
% r = roots(Vcoef);
% r = real(r);
% plot(r(1),Pc,'*');
% roooot = r(1);
% hold on;
% V = linspace(0,10);
% figure(1);
p=70; % Critical pressure from literature
fplot(p);
for T= 409:1:413
    i=T-408;
hold on;
syms v;
fplot(R*T/(v-b) - a/(v^2),[0 0.4]);
xlim([0 0.4])
ylim([58 70])
V = linspace(0.1,0.4); 
hold on
A = 4.28176;
B = 959.43;
C = -29.758;
p_sat = 57.669;
fplot(sp(i));
hold on
Prfunc = @(V) (R*T./(V - b) - a./(V.^2)-sp(i));
P = Prfunc(V);
cs = P.*circshift(P,-1,2);        % Product Negative At Zero-Crossings
xc = V(cs <= 0);                    % Values Of ‘x’ Near Zero Crossings
for k1 = 1:length(xc)
    zeros(k1) = fzero(Prfunc, xc(k1));      % Use ‘xc’ As Initial Zero Estimate
end
sv(i*2)=zeros(1); %saving the first and third root in the saturation region
sv(i*2+1)=zeros(3);
% plot(V,P);
zeros

hold on;
end

y=[62 63.32 63.32 63.91 63.91 64.57 64.57 65.22 65.22 65.87];
sv(1)=0.15
sv(11)=[];
x=sv;
plot(sv,y,'o');
xx = 0.15:.001:0.25;
yy = spline(x,y,xx);
p=plot(x,y,'o',xx,yy);
p(2).LineWidth=2;
hold on
xlabel ('Volume in L');
ylabel ('Pressure in Bar');
title ('P-V DIAGRAM OF FORMALDEHYDE FOR VANDERWAAL EQUATION')
%for i= 0:10
    %T = 165 + 9*i;
%A = 4.28176;
%B = 959.43;
%C = -29.758;
%P(i) = power(10,(A-(B/(C+T))));
%y = [P -(b*P-R*T) a -a*b];
%r = roots(y);
%v(1) = min(r);
%v(2) = max(r);
%plot([V(1), v(2)], [P,P])
%end

% Antoine equation parameters
% A = 4.28176;
% B = 959.43;
% C = -29.758;
% p_plot  = zeros(1,1000000);
% v1_plot = zeros(1,1000000);
% v2_plot = zeros(1,1000000);
% count = 1;
% for t = 164: 1.0: Tc
%     temp = (A - (B/(C + t)));
%     p_sat = power(10,(A - (B/(C + t))));
%     % 1, c2, c1, c0 are the coefficients of the cubic equation in volume i.e.
%     % V^3 + c2 * V^2 + c1 * V + c0 = 0
%      y = [p_sat -(b*p_sat+R*t) a -a*b];
%      r = roots(y)
%     v1_plot(count) = min(r);
%     v2_plot(count) = max(r);
%     r
%     p_plot(count) = p_sat;
%     count = count + 1;
% end
% fprintf('%d', count);
% x1 = zeros(1, count);
% x2 = zeros(1, count);
% y1 = zeros(1, count);
% y2 = zeros(1, count);
% for i = 1:1: count
%     x1(i) = v1_plot(i);
%     x2(i) = v2_plot(i);
%     y1(i) = p_plot(i);
%     y2(i) = p_plot(count -i + 1);
%     fprintf('x1 = %f x2 = %f\n', x1(i), x2(count - i +1));
% end
% plot(x1, y1);
% hold on;
% plot(x2, y2);