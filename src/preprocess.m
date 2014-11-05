function [Alpha0, NumGraph, f, sfit, Alpha, BW_Alpha, endpoint] = preprocess(myFileName, bfuncdef, bLower, bUpper, bStart, endpoint)

%   Inputs: 
%	myFileName: name of file containing plot image.
%	bfuncdef{k}: function definitions, in string form.
%	bLower{k}: lower bounds for parameters in function. Can be empty.
%	bUpper{k}: upper bounds for parameters in function. Can be empty.
%	bStart{k}: starting guesses for parameters in function. Can be empty.

clc
tic

%%Pre-processing
%Read in Image
Alpha0=im2double(imread(myFileName));
mycbcr=rgb2ycbcr(Alpha0);

NumGraph = length(bfuncdef);
IsLinedPaper = (std2(mycbcr(:,:,3)) > 0.008);

%Fill in missing inputs if user didn't specify
clear sfit f 
for k = 1:length(bfuncdef)
	eval(['f{' int2str(k) '} = fittype(''' bfuncdef{k} ''',''independent'',''x'');'])
	if (length(bStart{k}) == 0)
		bStart{k}=ones(1,numargs(f{k})-1);
	end
	if (length(bUpper{k}) == 0)
		bUpper{k}=ones(1,numargs(f{k})-1)*Inf;
	end
	if (length(bLower{k}) == 0)
		bLower{k}=ones(1,numargs(f{k})-1)*-Inf;
	end
	sfit(k) = fitoptions('Method','NonlinearLeastSquares','Lower',bLower{k},'Upper',bUpper{k},'Startpoint',bStart{k});
end

%A bit of error checking
if (length(endpoint > 0))
	if endpoint(3) > endpoint(4) 
		Temp = endpoint(3);
		endpoint(3) = endpoint(4);
		endpoint(4) = Temp;
	end
	if endpoint(1) > endpoint(2)
		Temp = endpoint(3);
		endpoint(3) = endpoint(4);
		endpoint(4) = Temp;
	end
end

%Conversion to B/W
if (IsLinedPaper)	%lined paper
	IsRed = ~(mycbcr(:,:,3) > mean2(mycbcr(:,:,3))*1.05);
	Alpha=((1-mycbcr(:,:,3)).^3+mycbcr(:,:,1));						%Get rid of blue
	Alpha=IsRed.*Alpha + ~IsRed.*graythresh(mycbcr(:,:,1))*1.15;	%Get rid of red
	t=5;		
else
	Alpha=rgb2gray(Alpha0);			%blank paper
	t = 3.5;	%Low number for low expected noise. 2 is too low, 5 is too high; 3.5 works well.
end
sbox=32;
for i = 0:0.5:(size(Alpha,1)/sbox-1)
for j = 0:0.5:(size(Alpha,2)/sbox-1)
	if std2(Alpha(1+i*sbox:sbox+i*sbox,1+j*sbox:sbox+j*sbox)) < t/255
		AThresh(2*i+1,2*j+1)=0;	
		unif_Alpha(2*i+1,2*j+1) = 0;
	else
		AThresh(2*i+1,2*j+1)=graythresh(Alpha(1+i*sbox:sbox+i*sbox,1+j*sbox:sbox+j*sbox));
		unif_Alpha(2*i+1,2*j+1) = 1;
	end
end
end
AThresh2 = imresize(AThresh,[size(Alpha,1) size(Alpha,2)], 'bilinear');
unif_Alpha = imresize(unif_Alpha,[size(Alpha,1) size(Alpha,2)],'nearest');
BW_Alpha = (Alpha >= AThresh2);

clear AThresh AThresh2 unif_Alpha IsRed

%Conservative small-region removal
BW_Alpha = 1-bwareaopen(1-BW_Alpha,100,8);