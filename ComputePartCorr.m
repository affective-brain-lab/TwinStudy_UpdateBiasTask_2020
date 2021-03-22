close all;
clear all;
clc;

%directory strucutre
dir_base = '';

cd(dir_base)
%load data divided in twinA and twinB (and diveded by zygosity)

data1 = xlsread('');
data2 = xlsread('');

%select the variable you want to correlate
twin1=data1(:,5);
twin2=data2(:,5);

%set the seed so that every time you run the code you have the same random
%numbers
rng('default');
rng(42);
N=length(twin1);
A = randi([1,2],10000*N,1);


age1=data1(:,3);
age2=data2(:,3);

EE1=data1(:,10); 
EE2=data2(:,10); 

fe1=data1(:,4);
fe2=data2(:,4);

mem1=data1(:,13);
mem2=data2(:,13);

exp1=data1(:,6);
exp2=data2(:,6);

fam1=data1(:,7);
fam2=data2(:,7);

viv1=data1(:,8);
viv2=data2(:,8);

neg1=data1(:,9);
neg2=data2(:,9);


AAllr=[];
AAllp=[];

%change position of twinA and twinB 1000 times
%each time change also the controlling variables
 for a=1:10000
 twin1dz = [];
 twin2dz = [];
%  
 EE1 = [];
 EE2 = [];

 Afe1 = [];
 Afe2 = [];
 
 Amem1 = [];
 Amem2 = [];
 
 Aexp1=[];
 Aexp2=[];
 
 Afam1=[];
 Afam2=[];
 
 Aviv1=[];
 Aviv2=[];
 
 Aneg1=[];
 Aneg2=[];
 
 

 %twinA and twinB and covariates are assigned to x or y 
    for i=1:size(twin1)
           x=A(i+N*(a-1));
            
    	if x==1
             twin1dz =[twin1dz,twin1(i)];
             twin2dz =[twin2dz,twin2(i)];
          
            Aexp1 =[Aexp1,exp1(i)];
            Aexp2 =[Aexp2,exp2(i)];
            
            Aneg1 =[Aneg1,neg1(i)];
            Aneg2 =[Aneg2,neg2(i)];
            
            Aviv1 =[Aviv1,viv1(i)];
            Aviv2 =[Aviv2,viv2(i)];
            
            Afam1 =[Afam1,fam1(i)];
            Afam2 =[Afam2,fam2(i)];
            
            EE1 =[EE1,EE1(i)];
            EE2 =[EE2,EE2(i)];
         
            Afe1 =[Afe1,fe1(i)];
            Afe2 =[Afe2,fe2(i)];
         
             Amem1 =[Amem1,mem1(i)];
             Amem2 =[Amem2,mem2(i)];
             
             
         elseif x==2
             twin1dz =[twin1dz,twin2(i)];
             twin2dz =[twin2dz,twin1(i)];
             
             EE1 =[EE1,EE2(i)];
             EE2 =[EE2,EE1(i)];

             Afe1 =[Afe1,fe2(i)];
             Afe2 =[Afe2,fe1(i)];

             Amem1 =[Amem1,mem2(i)];
             Amem2 =[Amem2,mem1(i)];
             
             Aexp1 =[Aexp1,exp2(i)];
            Aexp2 =[Aexp2,exp1(i)];
            
            Aneg1 =[Aneg1,neg2(i)];
            Aneg2 =[Aneg2,neg1(i)];
            
            Aviv1 =[Aviv1,viv2(i)];
            Aviv2 =[Aviv2,viv1(i)];
            
            Afam1 =[Afam1,fam2(i)];
            Afam2 =[Afam2,fam1(i)];
        end
    end    
 twin1dz = twin1dz';
 twin2dz = twin2dz';  
 
 EE1 = EE1';
 EE2 = EE2';
 
 Afe1 = Afe1';
 Afe2 =Afe2';
 
 Amem1 = Amem1';
 Amem2 = Amem2';
 
 Aexp1=Aexp1';
 Aexp2=Aexp2';
 
 Afam1=Afam1';
 Afam2=Afam2';
 
 Aviv1=Aviv1';
 Aviv2=Aviv2';
  
 Aneg1=Aneg1';
 Aneg2=Aneg2';
 
 %compute differece between twinA and twinB in the covariates  
 %first estimate of twin1-twin2
 Afe=Afe1-Afe2;
  
 AEE=EE1-EE2;
 Amem=Amem1-Amem2;
 Aexp=Aexp1-Aexp2;
 Afam=Afam1-Afam2;
 Aviv=Aviv1-Aviv2;
 Aneg=Aneg1-Aneg2;
 
 %run partial correlation
 a=[twin1dz twin2dz];
 z=[Afe age1 Amem Aexp  Afam Aviv Aneg AEE];
 
 [r1,p1] = partialcorr(a,z,'rows','pairwise');

 r=r1(1,2);
 p=p1(1,2);
 
 %print r coefficient
 AAllr=[AAllr;r];
 AAllp=[AAllp;p];
 AAmeanR=mean(AAllr);
 AAmeanP=mean(AAllp);
 


 end
  writematrix(AAllr,'Update_mz.xlsx') 
  