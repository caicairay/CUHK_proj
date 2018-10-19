import flmtdtct.*

clear all
close all
% addpath('/data/zitan/960/run960-6');
addpath('/data/zitan/240/240-4');
fnum=30;
flnm=['hdfaa.' sprintf('%03d',fnum)];
rho=hdf5read(flnm,'gas_density');
rhos=squeeze(rho(:,:,120));

result=flmtdtct(rhos);

sizeA=size(rhos);
x=1:sizeA(1);
y=1:sizeA(2);
[X,Y]=meshgrid(x,y);
%% figures
figure
contourf(rhos,100)
hold on
scatter(X(result),Y(result),'.r')
axis equal
axis tight
hold off
%isosurface(rhos,500);