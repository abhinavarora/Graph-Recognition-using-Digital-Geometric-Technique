function [Alpha1, ReduBW_Alpha, rystart, rxstart, ryend, rxend, ryint, rxint, Beta] = plotAreaIdentification(Alpha0, Alpha, BW_Alpha)

%%Find axes here; won't work if two lines closer to 90 deg than axes are
anglestart = -90;
[count, s]=radon(1-BW_Alpha,anglestart:0.5:90);
regmax = imregionalmax(count);
i=0; yt = 0;
while (length(yt)<3)	%Can be increased if expecting graphs with multiple lines, although code must be edited to do so
	[yt, angt] = find(regmax ~= 0 & count > max(max(count))*(0.75-i));
	i = i+0.01;
end
angt = deg2rad((angt - 1)/2 + anglestart);
for i = 1:length(angt)
	for j = 1:length(angt)
		angdifm(i,j)=abs(abs(angt(i)-angt(j))-deg2rad(90));
	end
end

%Determine which line is horizontal/vertical
[i, j] = find(angdifm==min(min(angdifm)));	%i(1) in case there is more than one minimum and i is a vector
if (abs(angt(i(1))) < abs(angt(j(1))))
	temp = i(1); i(1) = j(1); j(1) = temp; clear temp
end
a_ho = angt(i(1)); a_ve = angt(j(1));
s_ho = s(yt(i(1))); s_ve = s(yt(j(1)));
if (a_ho > 0) 
	a_ho = a_ho - pi; 
	s_ho = -s_ho;
end

clear count regmax s yt angt angdifm Alpha2

%If the lines are close to parallel, i.e. less than 30 degrees close, something's wrong and we should pause
if (abs(a_ve-a_ho) < pi/6)
	disp('Error, cannot find axes')
	pause
end

%%Projective transform on image
tp = (-sin(a_ho)*(s_ve*sin(a_ve)-s_ho*sin(a_ho))+cos(a_ho)*(s_ho*cos(a_ho)-s_ve*cos(a_ve)))./(sin(a_ve)*cos(a_ho)-sin(a_ho)*cos(a_ve));
t = (tp*sin(a_ve)+s_ve*cos(a_ve)-s_ho*cos(a_ho))./(sin(a_ho));
centers=(size(Alpha)+1)/2; yc=centers(1); xc=centers(2);

%(tp-100) should be in right direction for tilts on both side of axis, as long as they're between 45 deg of horiz/vert
xint = xc + tp*sin(a_ve)+s_ve*cos(a_ve);yint = yc + tp*cos(a_ve)-s_ve*sin(a_ve);
x_ve = xc + (tp-100)*sin(a_ve)+s_ve*cos(a_ve);y_ve = yc + (tp-100)*cos(a_ve)-s_ve*sin(a_ve);
x_ho = xc + (t-100)*sin(a_ho)+s_ho*cos(a_ho);y_ho = yc + (t-100)*cos(a_ho)-s_ho*sin(a_ho);
xint2 = x_ve + x_ho - xint; yint2 = y_ve + y_ho - yint;

RotBW_Alpha = 1-imtransform(1-BW_Alpha,maketform('projective',[x_ve y_ve; xint yint; xint2 yint2; x_ho y_ho], [xint yint-100; xint yint; xint+100 yint-100; xint+100 yint]),'nearest');
Alpha1 = 1-imtransform(1-Alpha0,maketform('projective',[x_ve y_ve; xint yint; xint2 yint2; x_ho y_ho], [xint yint-100; xint yint; xint+100 yint-100; xint+100 yint]),'nearest');

%This will actually fail if the lines are really close to parallel, but program should have stopped before this
Temp = zeros(size(Alpha)); Temp(round(yint),round(xint)) = 1;
Temp = imtransform(Temp,maketform('projective',[x_ve y_ve; xint yint; xint2 yint2; x_ho y_ho], [xint yint-100; xint yint; xint+100 yint-100; xint+100 yint]),'nearest');
[yint,xint]=find(Temp ~=0);	%Finds location of new intersection
yint=yint(1); xint=xint(1);	%in case there's more than one

clear y_ve x_ve y_ho x_ho xint2 yint2 Temp s_ho a_ho a_ve s_ve xc yc t tp Temp Alpha Alpha0

%Find bounds
ystart = 1;
yend = size(RotBW_Alpha,1);
xstart = 1;
xend = size(RotBW_Alpha,2);
for i = 0:(size(RotBW_Alpha,1)-yint-11)
	if (nnz(1-RotBW_Alpha((yint+i-10):(yint+i+10),(xint-10):(xint+10))) < 5)
		yend = yint + i;
		break;
	end
end
for i = 0:(yint-11)
	if (nnz(1-RotBW_Alpha((yint-i-10):(yint-i+10),(xint-10):(xint+10))) < 5)
		ystart = yint - i
		break;
	end
end
for i = 0:(size(RotBW_Alpha,2)-xint-11)
	if (nnz(1-RotBW_Alpha((yint-10):(yint+10),(xint+i-10):(xint+i+10))) < 5)
		xend = xint + i;
		break;
	end
end
for i = 0:(xint-11)
	if (nnz(1-RotBW_Alpha((yint-10):(yint+10),(xint-i-10):(xint-i+10))) < 5)
		xstart = xint - i
		break;
	end
end

%100 is the buffer on the sides, in case only one quadrant is included, to get the axis labels
ReduBW_Alpha = RotBW_Alpha(max((ystart-100),1):min((yend+100),size(RotBW_Alpha,1)),max((xstart-100),1):min((xend+100),size(RotBW_Alpha,2)));
Alpha1 = Alpha1(max((ystart-100),1):min((yend+100),size(RotBW_Alpha,1)),max((xstart-100),1):min((xend+100),size(RotBW_Alpha,2)),:);
rystart = (ystart > 100)*(100-1)+1;
rxstart = (xstart > 100)*(100-1)+1;
ryend = yend - max((ystart-100),1) + 1;
rxend = xend - max((xstart-100),1) + 1;
rxint = xint - max((xstart-100),1) + 1;
ryint = yint - max((ystart-100),1) + 1;

%%Find characters
Beta=1-ReduBW_Alpha;
CC = bwconncomp(Beta);
for i = 1:length(CC.PixelIdxList)
	if ((length(CC.PixelIdxList{i}) > 1000) || (length(CC.PixelIdxList{i}) < 100))
		Beta(CC.PixelIdxList{i})=0;
	end
end
Beta(1:(ryint+10),(rxint-10):size(Beta,2))=0;	%Only include characters close to lines -> scale should be changed to depend on image size...
Beta(:,1:(rxint-100))=0;
Beta((yint+200):size(Beta,1),:)=0;