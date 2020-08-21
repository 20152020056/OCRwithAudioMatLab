clear all
close all
clc

img = imread('BOGOTA.png');

figure, imshow(img)
grises=rgb2gray(img);

umbral=0.4;
BW = im2bw(grises,umbral);


m=BW;
[fil,col]=size(m);
t=1:256;
PENDIENTE=-5;
C=10;
L=(PENDIENTE*(t))+C;

ma=max(L);
L=(L./ma)*255;
L=round(L);
HA=L;



R=zeros(fil,col);
for f=1:fil
    for c=1:col
        g=m(f,c)+1;
        R(f,c)= HA(g);
    end
end



[L Ne]=bwlabel(R);
figure,imshow(label2rgb(L));

prop=regionprops(L);

length(prop)

hold on

for n=1:length(prop)
rectangle('position',prop(n).BoundingBox,'EdgeColor','g','LineWidth',2)
end

centroids=cat(1,prop.Centroid);
plot(centroids(:,1),centroids(:,2),'r*')

hold off

cont=length(prop);
imagen=R;
for nn=1:length(prop)
B1=prop(nn).BoundingBox;

A1=imagen(B1(2):B1(2)+B1(4),B1(1):B1(1)+B1(3));
figure, imshow(A1)

c=dct2(A1);
figure, imshow(uint8(c))
i=idct2(c(1:27,1:27));
vmax=max((max(i)));
vmin=min((min(i)));
d=vmax-vmin;
in=(((i-vmin)./d).*255);
figure,imshow(uint8(in))

[fil,col]=size(in);

v(fil,col)=0;
for i=1:1:col
    for j=1:1:fil
        v(j)=in(j,i)+v(j);
    end
end

v=v*(1/length(in));
t=(1:length(in));
length(t);
length(v);

AA1(fil,col)=0;
for i=1:1:fil
    for j=1:1:col
        AA1(j,i)=in(i,j);
    end
end

v2(fil,col)=0;
for i=1:1:col
    for j=1:1:fil
        v2(j)=AA1(j,i)+v2(j);
    end
end
v2=v2*(1/length(AA1));

W(1:length(v))=0;
LETRA=v+v2;
figure,plot(LETRA)

load A,load B,load C,load D,load E,load F 
load G, load H, load J, load K, load L
load M, load N, load O, load P, load Q 
load R, load S, load T, load U, load V
load W, load X, load Y, load Z

la = corrcoef(LETRA,A);PRO1=mean(la);co1=(PRO1(1)+PRO1(2))/2;
lb = corrcoef(LETRA,B);PRO2=mean(lb);co2=(PRO2(1)+PRO2(2))/2;
lc = corrcoef(LETRA,C);PRO3=mean(lc);co3=(PRO3(1)+PRO3(2))/2;
ld = corrcoef(LETRA,D);PRO4=mean(ld);co4=(PRO4(1)+PRO4(2))/2;
le = corrcoef(LETRA,E);PRO5=mean(le);co5=(PRO5(1)+PRO5(2))/2;
lf = corrcoef(LETRA,F);PRO6=mean(lf);co6=(PRO6(1)+PRO6(2))/2;
lg = corrcoef(LETRA,G);PRO7=mean(lg);co7=(PRO7(1)+PRO7(2))/2;
lh = corrcoef(LETRA,H);PRO8=mean(lh);co8=(PRO8(1)+PRO8(2))/2;
lj = corrcoef(LETRA,J);PRO9=mean(lj);co9=(PRO9(1)+PRO9(2))/2;
lk = corrcoef(LETRA,K);PRO10=mean(lk);co10=(PRO10(1)+PRO10(2))/2;
ll = corrcoef(LETRA,L);PRO11=mean(ll);co11=(PRO11(1)+PRO11(2))/2;
lm = corrcoef(LETRA,M);PRO12=mean(lm);co12=(PRO12(1)+PRO12(2))/2;
ln = corrcoef(LETRA,N);PRO13=mean(ln);co13=(PRO13(1)+PRO13(2))/2;
lo = corrcoef(LETRA,O);PRO14=mean(lo);co14=(PRO14(1)+PRO14(2))/2;
lp = corrcoef(LETRA,P);PRO15=mean(lp);co15=(PRO15(1)+PRO15(2))/2;
lr = corrcoef(LETRA,R);PRO16=mean(lr);co16=(PRO16(1)+PRO16(2))/2;
ls = corrcoef(LETRA,S);PRO17=mean(ls);co17=(PRO17(1)+PRO17(2))/2;
lt = corrcoef(LETRA,T);PRO18=mean(lt);co18=(PRO18(1)+PRO18(2))/2;
lu = corrcoef(LETRA,U);PRO19=mean(lu);co19=(PRO19(1)+PRO19(2))/2;
lv = corrcoef(LETRA,V);PRO20=mean(lv);co20=(PRO20(1)+PRO20(2))/2;
lw = corrcoef(LETRA,W);PRO21=mean(lw);co21=(PRO21(1)+PRO21(2))/2;
lx = corrcoef(LETRA,X);PRO22=mean(lx);co22=(PRO22(1)+PRO22(2))/2;
ly = corrcoef(LETRA,Y);PRO23=mean(ly);co23=(PRO23(1)+PRO23(2))/2;
lz = corrcoef(LETRA,Z);PRO24=mean(lz);co24=(PRO24(1)+PRO24(2))/2;
lq = corrcoef(LETRA,Q);PRO25=mean(lq);co25=(PRO25(1)+PRO25(2))/2;

if co1>=0.998
    Se(nn)=strcat('A');
    disp('A')
    
 
   [y,fs]=audioread('AudioA.wav');
     sound(y,fs)
plot(y)
end
if co2>=0.998
    Se(nn)=strcat('B');
    disp('B')
     [y,fs]=audioread('AudioB.wav');
     sound(y,fs)
     plot(y)
end
if co3>=0.998
    Se(nn)=strcat('C');
    disp('C')
     [y,fs]=audioread('AudioC.wav');
     sound(y,fs)
     plot(y)
end
if co4>=0.998
    Se(nn)=strcat('D');
    disp('D')
     [y,fs]=audioread('AudioD.wav');
     sound(y,fs)
     plot(y)
end
if co5>=0.998
    Se(nn)=strcat('E');
    disp('E')
     [y,fs]=audioread('AudioE.wav');
     sound(y,fs)
     plot(y)
end
if co6>=0.998
    Se(nn)=strcat('F');
    disp('F')
     [y,fs]=audioread('AudioF.wav');
     sound(y,fs)
     plot(y)
end
if co7>=0.998
    Se(nn)=strcat('G');
    disp('G')
     [y,fs]=audioread('AudioG.wav');
     sound(y,fs)
     plot(y)
end
if co8>=0.998
    Se(nn)=strcat('H');
    disp('H')
     [y,fs]=audioread('AudioH.wav');
     sound(y,fs)
     plot(y)
end
if co9>=0.998
    Se(nn)=strcat('J');
    disp('J')
     [y,fs]=audioread('AudioJ.wav');
     sound(y,fs)
     plot(y)
end
if co10>=0.998
    Se(nn)=strcat('K');
    disp('K')
     [y,fs]=audioread('AudioK.wav');
     sound(y,fs)
     plot(y)
end
if co11>=0.998
    Se(nn)=strcat('L');
    disp('L')
     [y,fs]=audioread('AudioL.wav');
     sound(y,fs)
     plot(y)
end
if co12>=0.998
    Se(nn)=strcat('M');
    disp('M')
     [y,fs]=audioread('AudioM.wav');
     sound(y,fs)
     plot(y)
end
if co13>=0.998
    Se(nn)=strcat('N');
    disp('N')
     [y,fs]=audioread('AudioN.wav');
     sound(y,fs)
     plot(y)
end
if co14>=0.998
    Se(nn)=strcat('O');
    disp('O')
     [y,fs]=audioread('AudioO.wav');
     sound(y,fs)
     plot(y)
     
end
if co15>=0.998
    Se(nn)=strcat('P');
    disp('P')
     [y,fs]=audioread('AudioP.wav');
     sound(y,fs)
     plot(y)
end
if co16>=0.998
    Se(nn)=strcat('R');
    disp('R')
     [y,fs]=audioread('AudioR.wav');
     sound(y,fs)
     plot(y)
end
if co17>=0.998
    Se(nn)=strcat('S');
    disp('S')
     [y,fs]=audioread('AudioS.wav');
     sound(y,fs)
     plot(y)
end
if co18>=0.998
    Se(nn)=strcat('T');
    disp('T')
     [y,fs]=audioread('AudioT.wav');
     sound(y,fs)
     plot(y)
end
if co19>=0.998
    Se(nn)=strcat('U');
    disp('U')
     [y,fs]=audioread('AudioU.wav');
     sound(y,fs)
     plot(y)
end
if co20>=0.998
    Se(nn)=strcat('V');
    disp('V')
     [y,fs]=audioread('AudioV.wav');
     sound(y,fs)
     plot(y)
end
if co21>=0.998
    Se(nn)=strcat('W');
    disp('W')
     [y,fs]=audioread('AudioW.wav');
     sound(y,fs)
     plot(y)
end
if co22>=0.998
    Se(nn)=strcat('X');
    disp('X')
     [y,fs]=audioread('AudioX.wav');
     sound(y,fs)
     plot(y)
end
if co23>=0.998
    Se(nn)=strcat('Y');
    disp('Y')
     [y,fs]=audioread('AudioY.wav');
     sound(y,fs)
     plot(y)
end
if co24>=0.998
    Se(nn)=strcat('Z');
    disp('Z')
     [y,fs]=audioread('AudioZ.wav');
     sound(y,fs)
     plot(y)
end
if co25>=0.998
    Se(nn)=strcat('Q');
    disp('Q')
     [y,fs]=audioread('AudioQ.wav');
     sound(y,fs)
     plot(y)
end
end 

palabra=char(Se)

 








