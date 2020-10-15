RGB = imread('Zdj2.jpg');
[A,B] = MaskaNiebieski(RGB);
% Pewnie bedzie trzeba tutaj jakies przeksztalcenia, pytanie czy sama
% rotacja wystarczy, poki co po prostu odczytam pod jakim katem jest ta
% linia i o taki kat obroce

RGB = imrotate(RGB,2);
A = imrotate(A,2);



A=imopen(A,strel('cube',17));
A=imclose(A,strel('cube',17));
%imshow(A);figure;imshow(B);figure;imshow(C);
%imshow(A);
[H,T,R] = hough(A);
%houghMatViz(H,T,R);
peaks = houghpeaks(H,20,'NHoodSize',[401 31]);
size(peaks,1);
lines = houghlines(A,T,R,peaks);

imshow(RGB)
hold on
numLines = length(lines);
for k = 1:numLines  %5:6  %roboczo, nie jestem pewien jak wybraæ linie 
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','red','Marker', 'x', 'MarkerEdgeColor', 'g');
end

%%
GornaLiniaLewo = lines(6).point1;
GornaLiniaPrawo = lines(6).point2;
DolnaLiniaLewo = lines(5).point1;
DolnaLiniaPrawo = lines(5).point2;

DlKratkiX_Gora = (GornaLiniaPrawo(1,1) - GornaLiniaLewo(1,1))/11;
DlKratkiX_Dol = (DolnaLiniaPrawo(1,1) - DolnaLiniaLewo(1,1))/11;
DlKratkiY_Lewo = (DolnaLiniaLewo(1,2) - GornaLiniaLewo(1,2))/3;
DlKratkiY_Prawo = (DolnaLiniaPrawo(1,2) - GornaLiniaPrawo(1,2))/3;
X = zeros(4,12);
Y = zeros(4,12);
for i = 1:12
    XKratki = GornaLiniaLewo(1,1) + (i-1)*DlKratkiX_Gora;
    for j = 1:4
        X(j,i) = XKratki;
    end
end
for k = 1:4
    YKratki = GornaLiniaLewo(1,2)+ (k-1)*DlKratkiY_Lewo;
    for l = 1:12
         Y(k,l) = YKratki;
    end
end


imshow(RGB)
hold on
plot(X,Y,'LineWidth',2,'Color','red','Marker', 'x', 'MarkerEdgeColor', 'g');
hold on
plot(X',Y','LineWidth',2,'Color','red','Marker', 'x', 'MarkerEdgeColor', 'g');
%%
RGB2 = imread('Zdj22.jpg');
RGB2=imrotate(RGB2,2);
RGB2Gray = rgb2gray(RGB2);
I = imbinarize(RGB2Gray,0.99);
I=imopen(I,strel('cube',15));
imshow(RGB2)
hold on
plot(X,Y,'LineWidth',2,'Color','red','Marker', 'x', 'MarkerEdgeColor', 'g');
hold on
plot(X',Y','LineWidth',2,'Color','red','Marker', 'x', 'MarkerEdgeColor', 'g');

[im_lb, num] = bwlabel(I,4);
feats = regionprops(im_lb,'all');
hold on
for k = 1:length(feats)
 x=feats(k).Centroid(1,1);
 y=feats(k).Centroid(1,2);
 plot(x,y,'*')
end
%%
for k = 1:length(feats)
 %plot(feats(k).Centroid(1),feats(k).Centroid(2),'or'); hold on
A=feats(k).Image;
 
m00 = 0;
m01 = 0;
m11 = 0;
m10 = 0;
m02 = 0;
m20 = 0;
m12 = 0;
m21 = 0;
m03 = 0;
m30 = 0;

width = size(A,2);
height = size(A,1);
for x = 1:width
 for y = 1:height
 m00 = m00 + A(y,x)*x^0*y^0;
 m01 = m01 + A(y,x)*x^0*y^1;
 m11 = m11 + A(y,x)*x^1*y^1;
 m10 = m10 + A(y,x)*x^1*y^0;
 m02 = m02 + A(y,x)*x^0*y^2;
 m20 = m20 + A(y,x)*x^2*y^0;
 m12 = m12 + A(y,x)*x^1*y^2;
 m21 = m21 + A(y,x)*x^2*y^1;
 m03 = m03 + A(y,x)*x^0*y^3;
 m30 = m30 + A(y,x)*x^3*y^0;

Xs=m10/m00;
Ys=m01/m00;
M20=m20-(m10^2/m00);
M02=m02-(m01^2/m00);
M11=m11-(m01*m10/m00);
M30=m30-3*m20*Xs+2*m10*Xs^2;
M03=m03-3*m02*Ys+2*m01*Ys^2;
M21=m21-2*m11*Xs-m20*Ys+2*m01*Xs^2;
M12=m12-2*m11*Ys-m02*Xs+2*m10*Ys^2;

N00 = 1;
N01 = 0;
N10 = 0;
N20 = M20 / m00^2;
N02 = M02 / m00^2;
N12 = M12 / m00^((3/2)+1);
N21 = M21 / m00^((3/2)+1);
N03 = M03 / m00^((3/2)+1);
N30 = M30 / m00^((3/2)+1);
N11 = M11 / m00^2;


I0=N20+N02;
I1=(N20-N02)^2+4*N11^2;
I2=(N30-3*N12)^2+(3*N21-N03)^2;
I3=(N30+N12)^2+(N21+N03)^2;
I4=(N30-3*N12)*(N30+N12)*[(N30+N12)^2-3*(N21*N03)^2];
I5=(N20-N02)*[(N30+N12)^2-(N21+N03)^2]+4*N11*(N30+N12)*(N21+N03);

%feats(k).M00 = M00;
%feats(k).M01 = M01;
feats(k).M11 = M11;
%feats(k).M10 = M10;
feats(k).M02 = M02;
feats(k).M20 = M20;
feats(k).M12 = M12;
feats(k).M21 = M21;
feats(k).M03 = M03;
feats(k).M30 = M30;

feats(k).N00 = N00;
feats(k).N01 = N01;
feats(k).N11 = N11;
feats(k).N10 = N10;
feats(k).N02 = N02;
feats(k).N20 = N20;
feats(k).N12 = N12;
feats(k).N21 = N21;
feats(k).N03 = N03;
feats(k).N30 = N30;

feats(k).I0 = I0;
feats(k).I1 = I1;
feats(k).I2 = I2;
feats(k).I3 = I3;
feats(k).I4 = I4;
feats(k).I5 = I5;

end
end
end
%%
WspolSrKratki=[0 0 0 0];
n=1;
for i=1:length(feats)
    for k=1:3
        for j=1:11
            if (feats(i).Centroid(1,1) > X(1,j) && feats(i).Centroid(1,1) < X(1,j+1) && feats(i).Centroid(1,2) > Y(k,1) && feats(i).Centroid(1,2) < Y(k+1,1))
                WspolSrKratki(n,1) = X(1,j) + (X(1,j+1)-X(1,j))/2;
                WspolSrKratki(n,2) = Y(k,1)+ (Y(k+1,1)-Y(k,1))/2;
                WspolSrKratki(n,3) = feats(i).Area;
                n=n+1;
            end
        end
    end
    
    ang=0:0.01:2*pi;
    if(WspolSrKratki(i,3)>2500 && WspolSrKratki(i,3)<4000) %Jakas cecha decydujaca co to jest
xp=((X(1,j+1)-X(1,j))/2 - 20)*cos(ang);
yp=((X(1,j+1)-X(1,j))/2 - 20)*sin(ang);
plot(WspolSrKratki(i,1)+xp,WspolSrKratki(i,2)+yp,'r');
    end
    
    if(WspolSrKratki(i,3)>6000 && WspolSrKratki(i,3)<7000) %Jakas cecha decydujaca co to jest
xp=((X(1,j+1)-X(1,j))/2 - 20)*cos(ang);
yp=((X(1,j+1)-X(1,j))/2 - 20)*sin(ang);
plot(WspolSrKratki(i,1)+xp,WspolSrKratki(i,2)+yp,'b');
    end
    
    if(WspolSrKratki(i,3)>8000 && WspolSrKratki(i,3)<12000) %Jakas cecha decydujaca co to jest
xp=((X(1,j+1)-X(1,j))/2 - 20)*cos(ang);
yp=((X(1,j+1)-X(1,j))/2 - 20)*sin(ang);
plot(WspolSrKratki(i,1)+xp,WspolSrKratki(i,2)+yp,'g');
    end
end
