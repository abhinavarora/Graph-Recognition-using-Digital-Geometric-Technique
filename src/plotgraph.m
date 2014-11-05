function TotRMSE = plotgraph(NumGraph, f, sfit, Alpha1, ryint, rxint, Clean_Alpha, Label_Alpha, endpoint)

%Build possibility matrix
Unique_Alpha = unique(Label_Alpha);
if (length(Unique_Alpha) == (NumGraph+1))
	PMatrix = ones(NumGraph,NumGraph);
else
	PMatrix = eye(length(Unique_Alpha)-1,length(Unique_Alpha)-1);
	for k = 2:size(Label_Alpha,2)
		T_unique = unique(Label_Alpha(:,k));
		if (length(setxor(T_unique,unique(Label_Alpha(:,k-1)))) > 0) && (length(T_unique)>2)
			for i = 2:(length(T_unique)-1)
				for j = (i+1):length(T_unique)
					PMatrix(T_unique(j),T_unique(i)) = 1;
				end
			end
		end
	end
end

%Set scales
if (length(endpoint) == 0)
	yScale = 100;
	xScale = 100;
else
	yScale = size(Clean_Alpha,1)/(endpoint(2)-endpoint(1));
	xScale = size(Clean_Alpha,2)/(endpoint(4)-endpoint(3));
end

%%Get data
Attempts = [];
minrmse = Inf;
for n = 1:min(5*factorial(NumGraph)+30,1000)
	Attempted = 0;
	
	Groups = [];
	Lengths = ones(NumGraph,1);
	for i = 1:(length(Unique_Alpha)-1)
		Temp = randperm(NumGraph);
		for j = 1:(i-1)
			if (PMatrix(i,j) == 1)
				[Ti,Tj] = find(Groups == j);
				Temp(find(Temp==Ti)) = [];
			end
		end
		
		%Ideally changes should be made to the algorithm here so that a Temp of length 0 never happens, 
		%i. e. don't place the labels in fits so that later labels have nowhere to go, but such a change would be really tedious.
		%This should work almost as well as long as NumGraph < 10 or so.
		if length(Temp ~= 0)
			Groups(Temp(1),Lengths(Temp(1))) = i;
			Lengths(Temp(1)) = Lengths(Temp(1)) + 1;
		else
			Attempted = 1	%Wasn't actually attempted before, but this is a failed attempt so we don't want it to do anything
			break
			n = n - 1;		%Don't count failed attempts
		end
	end
	
	%Get a hash value for this attempt
	Temp = [];
	for i = 1:NumGraph
		Temp = strcat(Temp,int2str(Groups(i,:)));
	end
	%See if this combination was already attempted
	for i = 1:length(Attempts)
		if ((length(Attempts{i}) == length(Temp)) && all(Attempts{i} == Temp))
			Attempted = 1;
			break;
		end
	end
	
	%Try next possibility if this combination already attempted
	if (Attempted == 1)
		%Skip fitting
	else
		Attempts{end+1}	= Temp;
		%Do fitting
		for i = 1:NumGraph
			py = [];
			px = [];
			for j = 1:(Lengths(i)-1)
				[ty, tx]=find(Label_Alpha==Groups(i,j));
				%Get rid of repeated pixels in order to not weight too heavily
				ty = (ryint - ty) / yScale;
				tx = (tx - rxint) / xScale;
				py = [py; ty];
				px = [px; tx];
			end
			%Matlab can generate some weird errors with fit(), but they can usually be ignored (just means iteration failed)
			try
				[c2{i},gof2{i}] = fit(px,py,f{i},sfit(i));
			catch
				disp('Warning: Fitting encountered at least one error.')
			end
		end
		Temp = 0;
		for i = 1:length(gof2)
			Temp = Temp + gof2{i}.rmse;
		end
		
		if (Temp < minrmse)
			minrmse = Temp;
			for i = 1:length(gof2)
				gof2min{i} = gof2{i};
				c2min{i} = c2{i};
				Groupsmin = Groups;
				Lengthsmin = Lengths;
			end
		end
	end
end
clear Attempted Attempts n Groups Lengths PMatrix c2 gof2

%Get colours of lines in plot
Color = zeros(NumGraph,3);
Temp = rgb2ycbcr(Alpha1);
for k = 1:NumGraph
	FTempi = [];
	FTempj = [];
	for i = 1:(Lengthsmin(k)-1)
		[FTempii,FTempjj] = find(Label_Alpha==Groupsmin(k,i));
	end
	FTempi = [FTempi; FTempii];
	FTempj = [FTempj; FTempjj];
	for j = 1:length(FTempi)
		Color(k,1) = Color(k,1) + Temp(FTempi(j),FTempj(j),1);
		Color(k,2) = Color(k,2) + Temp(FTempi(j),FTempj(j),2);
		Color(k,3) = Color(k,3) + Temp(FTempi(j),FTempj(j),3);
	end
	Color(k,:) = Color(k,:)/length(FTempi);
end
Color = ycbcr2rgb([0.45*ones(NumGraph,1) Color(:,2:3)])

toc	%Program finished computing

%%Plot figure
plot(1,1)	%Reset plot
hFig = figure(1);
set(hFig, 'Position', [1 1 size(Label_Alpha)])
hold off
axis([(1-rxint)/xScale (size(Label_Alpha,2)-rxint)/xScale -(size(Label_Alpha,1)-ryint)/yScale -(1-ryint)/yScale])
grid off
hold on
for i = 1:NumGraph
	%Can be uncommented to plot the data as well as lines
	%py = [];
	%px = [];
	%for j = 1:(Lengthsmin(i)-1)
	%	[ty, tx]=find(Label_Alpha==Groupsmin(i,j));
	%	ty = (ryint - ty) / yScale;
	%	tx = (tx - rxint) / xScale;
	%	py = [py; ty];
	%	px = [px; tx];
	%end
	p=plot(c2min{i})
	set(p,'Color',Color(i,:))
	set(p,'LineWidth',3)
	%p=plot(px,py,'.')
	%set(p,'Color',Color(i,:))
	%set(p,'MarkerSize',1)
end
if (length(endpoint) ~= 0)	%If bounds are provided, plot axes if bounds aren't 0
	if ~((endpoint(3)==0) | (endpoint(4) == 0))
		Temp=(((1-rxint)/xScale)-5):1:(((size(Label_Alpha,2)-rxint)/xScale)+5);	%Plot x-axis
		p=plot(Temp,zeros(1,length(Temp)),'k')
		set(p,'LineWidth',2)
	end
	if ~((endpoint(1)==0) | (endpoint(2) == 0))
		Temp=((-(size(Label_Alpha,1)-ryint)/yScale)-5):1:((-(1-ryint)/yScale)+5); %Plot y-axis
		p=plot(zeros(1,length(Temp)),Temp,'k')
		set(p,'LineWidth',2)
	end
else (length(endpoint) == 0)	%If graph is normalised, always manually plot axes
	axis off
	Temp=(((1-rxint)/xScale)-5):1:(((size(Label_Alpha,2)-rxint)/xScale)+5);	%Plot x-axis
	p=plot(Temp,zeros(1,length(Temp)),'k')
	set(p,'LineWidth',2)
	Temp=((-(size(Label_Alpha,1)-ryint)/yScale)-5):1:((-(1-ryint)/yScale)+5); %Plot y-axis
	p=plot(zeros(1,length(Temp)),Temp,'k')
	set(p,'LineWidth',2)
end
legend('off')
hold off
tempaxes = gca;
F = getframe;
pause

%Super-impose figure on original image
Temp = imresize(F.cdata,size(Label_Alpha),'nearest');
Temp2 = Temp==255;
Temp = (255-(uint8(Temp.*uint8(1-Temp2))));	%Invert the colours so they show up better
Preview = uint8(double(Temp).*(Temp2==0) + (Temp2==1).*Alpha1*255);
imshow(Preview(10:(end-10),10:(end-10),:))

%Give parameters
for k = 1:length(c2min)
	c2min{k}
	gof2min{k}
end

%Give total rmse
TotRMSE = 0;
for k = 1:length(c2min)
	c2min{k}
	gof2min{k}
	TotRMSE = TotRMSE + gof2min{k}.rmse;
end