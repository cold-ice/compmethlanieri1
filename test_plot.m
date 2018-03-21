clear all
close all
clc

test=0;
arch=0;

if arch==0
    archs='16';
else
    archs='32';
end
As= sprintf('~/Documents/Eurecom/CompMeth/assignment1/work/tests/16_%d_%d.txt', test, arch);
Bs= sprintf('~/Documents/Eurecom/CompMeth/assignment1/work/tests/64_%d_%d.txt', test, arch);
Cs= sprintf('~/Documents/Eurecom/CompMeth/assignment1/work/tests/256_%d_%d.txt', test, arch);
Ds= sprintf('~/Documents/Eurecom/CompMeth/assignment1/work/tests/1024_%d_%d.txt', test, arch);
Es= sprintf('~/Documents/Eurecom/CompMeth/assignment1/work/tests/4096_%d_%d.txt', test, arch);

A=readtable(As);
B=readtable(Bs);
C=readtable(Cs);
D=readtable(Ds);
E=readtable(Es);

x1=A{:,1};
y1=A{:,2};

x2=B{:,1};
y2=B{:,2};

x3=C{:,1};
y3=C{:,2};

x4=D{:,1};
y4=D{:,2};

x5=E{:,1};
y5=E{:,2};

limit=ones(length(x5)).*30;

fig=figure;
plot(x1, y1)

if(test==0)
    title('Distortion plots - Cosine signal')
    tests='cos';
elseif(test==1)
    title('Distortion plots - QPSK signal')
    tests='QPSK';
elseif(test==2)
    title('Distortion plots - 16QAM signal')
    tests='16QAM';
elseif(test==3)
    title('Distortion plots - WGN signal')
    tests='WGN';
end
    
xlabel('Input magnitude [dB]')
ylabel('Max SNR [dB]')
hold on
grid on

plot(x2, y2)
plot(x3, y3)
plot(x4, y4)
plot(x5, y5)
plot(x5, limit, 'black')
legend1=legend('16 points', '64 points', '256 points', '1024 points', '4096 points')
set(legend1,...
    'Position',[0.144047622879346 0.130158735480573 0.21964285331113 0.1928571375353]);
hold off

john=sprintf('%s_%s', tests, archs);
print(john, '-dpng')