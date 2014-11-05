%%This is not part of image processing function. 
%However, this demonstrates one possible way to obtain the inputs to the function.

prompt = {'Enter file name','How many equations will you be entering?'};
dlg_title = 'Input equation #';
num_lines = 1;
def = {'graph_3.jpg','2'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
myFileName = answer{1};
NumGraph = str2num(answer{2});
endpoint = [];

for k = 1:NumGraph
	prompt = {'Enter equation','Enter lower bounds (optional)','Enter upper bounds (optional)','Enter starting point (optional)'};
	dlg_title = 'Input equation #';
	num_lines = 1;
	def = {'A*x+B','-Inf -Inf','Inf Inf','1 1'};
	answer = inputdlg(prompt,dlg_title,num_lines,def);

	bfuncdef{k} = answer{1};
	bLower{k} = str2num(answer{2});
	bUpper{k} = str2num(answer{3});
	bStart{k} = str2num(answer{4});
end
clear NumGraph answer def dlg_title k num_lines prompt

% Step 1 - Preprocessing
[Alpha0, NumGraph, f, sfit, Alpha, BW_Alpha, endpoint] = preprocess(myFileName, bfuncdef, bLower, bUpper, bStart, endpoint)

% Step 2 - Plot Area Identification
[Alpha1, ReduBW_Alpha, rystart, rxstart, ryend, rxend, ryint, rxint, Beta] = plotAreaIdentification(Alpha0, Alpha, BW_Alpha)

%%Clean up image 
%Change image into thin 4-connected line segments with branch points identified
Alpha1=Alpha1(rystart:ryend,rxstart:rxend,:);

Clean_Alpha=bwareaopen((ReduBW_Alpha+Beta),20,4);
Clean_Alpha=Clean_Alpha(rystart:ryend,rxstart:rxend);
rxint = rxint - rxstart + 1;
ryint = ryint - rystart + 1;
Clean_Alpha = bwmorph(1-Clean_Alpha,'thin',Inf);

BTemp=bwmorph(Clean_Alpha,'branchpoints');
for i = 1:11
	ETemp=bwmorph(Clean_Alpha,'endpoints');
	Clean_Alpha = Clean_Alpha - (Clean_Alpha & ETemp);
end
Clean_Alpha = bwareaopen(Clean_Alpha,20);

%Remove corners
for i = 2:(size(Clean_Alpha,1)-1)
for j = 2:(size(Clean_Alpha,2)-1)
if (Clean_Alpha(i,j) == 1)
	ZTemp=Clean_Alpha((i-1:i+1),(j-1:j+1));
	if (mean2(ZTemp == [0 1 0; 1 1 0; 0 0 0])==1) Clean_Alpha(i,j) = 0; end
	if (mean2(ZTemp == [0 1 0; 0 1 1; 0 0 0])==1) Clean_Alpha(i,j) = 0; end
	if (mean2(ZTemp == [0 0 0; 1 1 0; 0 1 0])==1) Clean_Alpha(i,j) = 0; end
	if (mean2(ZTemp == [0 0 0; 0 1 1; 0 1 0])==1) Clean_Alpha(i,j) = 0; end
	
	if (mean2(ZTemp == [0 0 0; 0 1 1; 1 1 0])==1) Clean_Alpha(i,j) = 0; end
	if (mean2(ZTemp == [0 0 0; 1 1 0; 0 1 1])==1) Clean_Alpha(i,j) = 0; end
	if (mean2(ZTemp == [1 0 0; 1 1 0; 0 1 0])==1) Clean_Alpha(i,j) = 0; end
	if (mean2(ZTemp == [0 1 0; 1 1 0; 1 0 0])==1) Clean_Alpha(i,j) = 0; end
end
end
end
clear ETemp BTemp ZTemp


%Separate axes and curves
%%WILL NOT REMOVE BRANCH POINT IF AXIS ENDS PRECISELY AT BRANCH POINT (very unlikely)
BTemp = zeros(size(Clean_Alpha));
BTemp(ryint,rxint) = 1; %To separate x- and y-axes from being connected
BTemp = imdilate(BTemp,strel(ones(ceil(size(Alpha1,1)/40),ceil(size(Alpha1,2)/40))));

BTemp=BTemp | bwmorph(Clean_Alpha,'branchpoints');
BTemp = imdilate(BTemp,strel([1 1 1; 1 1 1; 1 1 1]));

%Remove axes
Label_Alpha=zeros(size(Clean_Alpha));
cc = bwconncomp(Clean_Alpha - (Clean_Alpha & BTemp));
for k = 1:cc.NumObjects
	[i,j]=ind2sub(size(Clean_Alpha),cc.PixelIdxList{k});
	%abs(mean(i-(ryint-rystart+1)))
	if (abs(mean(j-rxint)) > 20) && (abs(mean(i-ryint)) > 20)
		Label_Alpha(cc.PixelIdxList{k}) = k;
		%abs(mean(j-(rxint-rxstart+1)))
	end
end
Label_Alpha = bwlabel(Label_Alpha & ones(size(Label_Alpha)));

clear ReduBW_Alpha RotBW_Alpha BW_Alpha Alpha0 Beta BTemp mycbcr cc
% Step 4 - Data Labelling
Label_Alpha = datalabelling(NumGraph, Alpha1, Label_Alpha);
% Step 5 - Plotting
TotRMSE = plotgraph(NumGraph, f, sfit, Alpha1, ryint, rxint, Clean_Alpha, Label_Alpha, endpoint);

TotRMSE