function [Clean_Alpha, Label_Alpha] = cleanNsplit(Alpha1, ReduBW_Alpha, rystart, rxstart, ryend, rxend ryint, rxint, Beta)

%%Clean up image 
%Change image into thin 4-connected line segments with branch points identified
Alpha1=Alpha1(rystart:ryend,rxstart:rxend,:);

Clean_Alpha=bwareaopen((ReduBW_Alpha+Beta),20,4);
Clean_Alpha=Clean_Alpha(rystart:ryend,rxstart:rxend);
rxint = rxint - rxstart + 1;
ryint = ryint - rystart + 1;
Clean_Alpha = bwmorph(1-Clean_Alpha,'thin',Inf);
Label_Alpha = [];