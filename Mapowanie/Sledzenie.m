RGB = imread('11.jpg');
[A1,B] = MaskaNiebieskaSledzenie(RGB);
%imshow(A1);figure;
A=imclose(A1,strel('cube',25));
CC = imerode(A,strel('cube',23));
CC = imerode(CC,strel('cube',13));
CC = imerode(CC,strel('cube',5));
A = CC;
%imshow(A);figure;
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
Lewa = [0 0];
Prawa = [0 0];
j=1;
k=1;
for i= 1:length(lines)
    if (lines(i).rho < size(A,2)/2)
        Lewa(j,1) = lines(i).theta;
        Lewa(j,2) = lines(i).rho;
        j=j+1;
    end
    if (lines(i).rho > size(A,2)/2)
        Prawa(k,1) = lines(i).theta;
        Prawa(k,2) = lines(i).rho;
        k=k+1;
    end
end
