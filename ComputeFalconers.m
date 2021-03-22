close all;
clear all;
clc;



%directory strucutre


cd(dir_base)
%load data divided in twinA and twinB (and diveded by zygosity)

dz = xlsread('');
mz = xlsread('');
data(:,1)=mz;
data(:,2)=dz;
data(:,3)=2*(data(:,1)-data(:,2)); 
data(:,4)=2*data(:,2)-data(:,1); 
data(:,5)=1-data(:,1);

