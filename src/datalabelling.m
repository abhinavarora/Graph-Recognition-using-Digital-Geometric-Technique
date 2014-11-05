function Label_Alpha = datalabelling(NumGraph, Alpha1, Label_Alpha)

%%Matching 
%Matching when only one curve breaks and rest are present (works for n graphs)
%Conservative, should never go wrong
if ((length(unique(Label_Alpha))-1) > NumGraph)
	PActivePoints = [];
	for i = 1:size(Label_Alpha,2)
		if (length(unique(Label_Alpha(:,i))) == NumGraph+1)
			ActivePoints = unique(Label_Alpha(:,i));
			if ((length(PActivePoints) ~= 0) && ~all(ActivePoints == PActivePoints) && (length(setxor(ActivePoints,PActivePoints)) == 2))
				Temp=setxor(ActivePoints,PActivePoints);
				Label_Alpha(find(Label_Alpha==Temp(1)))=Temp(2);
				PActivePoints = [setxor(ActivePoints,PActivePoints); Temp(2)];
			else
				PActivePoints = ActivePoints;
			end
		end
    end
end

%Match colors in ycbcr domain
if ((length(unique(Label_Alpha))-1) > NumGraph)
	Unique_Alpha=unique(Label_Alpha);
	Color = zeros((length(Unique_Alpha)-1),3);
	Temp = rgb2ycbcr(Alpha1);
	for k = 2:length(Unique_Alpha)
		[FTempi,FTempj] = find(Label_Alpha==Unique_Alpha(k));
		for j = 1:length(FTempi)
			Color(k-1,1) = Color(k-1,1) + Temp(FTempi(j),FTempj(j),1);
			Color(k-1,2) = Color(k-1,2) + Temp(FTempi(j),FTempj(j),2);
			Color(k-1,3) = Color(k-1,3) + Temp(FTempi(j),FTempj(j),3);
		end
		Color(k-1,:) = Color(k-1,:)/length(FTempi);
	end

	%Only match if colors aren't all the same
	if ((std(Color(:,2)) < 0.04) && (std(Color(:,3)) < 0.04))
		%Do nothing because all lines are the same colour
		IsColor = 0;
	else
		IsColor = 1;
		%Match together lines based on color. Kmeans is actually a terrible way to do this since it often screws up for small data sets...
		[idx1,C] = kmeans(Color(:,2:3),NumGraph,'replicates',50)
		for k = 1:length(Color)
			sum((Color(k,2:3)-C(idx1(k),:)).^2);
			if (sum((Color(k,2:3)-C(idx1(k),:)).^2) < 0.0004)	%Each distance max of 0.0141 apart, e.g. 0.6 and 0.6141
				Label_Alpha(find(Label_Alpha == Unique_Alpha(k+1))) = idx1(k)+max(Unique_Alpha);
			end
		end
	end
end

%Match endpoint slopes and positions
%Can go wrong for special cases
if ((length(unique(Label_Alpha))-1) > NumGraph)
	ETemp = bwmorph(Label_Alpha,'endpoints');
	[Ti,Tj] = find(ETemp ~= 0); %Find the positions of every endpoint
	Matched = zeros(length(Ti),1);
	Angle = zeros(length(Ti),1);

	for k = 1:length(Ti)	
		%See if endpoint is already matched to another endpoint in the same vicinity
		Value = Label_Alpha(Ti(k),Tj(k));
		LTemp = Label_Alpha((Ti(k)-5):(Ti(k)+5),(Tj(k)-5):(Tj(k)+5));
		LTemp(find(LTemp~=Value)) = 0;
		LTemp = bwlabel(LTemp);
		Matched(k)= (length(unique(LTemp)) > 2);
		
		%Find the angles of every endpoint
		for j = 5:-1:1
			LTemp = bwlabel(Label_Alpha((Ti(k)-j):(Ti(k)+j),(Tj(k)-j):(Tj(k)+j)));
			Value=LTemp(j+1,j+1);
			LTemp(2:(end-1),2:(end-1))=0;
			[FTempi FTempj]=find(LTemp==Value);
			if (length(FTempi) > 0)
				Angle(k) = atan(-(mean(FTempi)-j-1)/(mean(FTempj)-j-1));
				if ((mean(FTempj)-j-1) < 0) 
					Angle(k) = Angle(k) - sign(Angle(k))*pi - (Angle(k) == 0)*pi; 
				end
				break
			end
		end
	end
	for k = 1:length(Ti)
		minangle = pi;
		for j = (k+1):length(Ti)	%For each other endpoint
			%If this is the smallest angle, and the other endpoint wasn't already matched
			if ((abs(abs(Angle(j) - Angle(k))-pi)<minangle) && (abs(Ti(k) - Ti(j)).^2 + abs(Tj(k) - Tj(j)).^2 < 15^2) && (Matched(j) ~= 1))
				jmin = j;
				minangle=abs(abs(Angle(j) - Angle(k))-pi);
			end
		end	
		if ((minangle < pi/4));
			Matched(jmin) = 1;
			Label_Alpha(find(Label_Alpha == Label_Alpha(Ti(jmin),Tj(jmin)))) = Label_Alpha(Ti(k),Tj(k));
		end
	end
end

CC = bwconncomp(Label_Alpha);
for i = 1:length(CC.PixelIdxList)
	if (length(CC.PixelIdxList{i}) < 40)
		Label_Alpha(CC.PixelIdxList{i})=0;
	end
end

if (length(unique(Label_Alpha))-1 < NumGraph)
	disp('Error, not enough graphs to fit your functions')
	pause
end

%Reduce elements to range of 1 to {# of different segments}
Unique_Alpha = unique(Label_Alpha);
for i = 2:length(Unique_Alpha)
	Label_Alpha(find(Label_Alpha==Unique_Alpha(i))) = (i-1);
end

clear Ti Tj Angle minangle jmin Unique_Alpha ActivePoints Matched Color