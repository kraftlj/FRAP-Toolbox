%%% FRAP Toolbox
%%
%     FRAP Toolbox is designed to be a modular software program designed
%     for the purposes of analyzing Fluorescence Recovery After
%     Photobleaching (FRAP) data. Copyright (C) 2011  Lewis J. Kraft
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or (at
%     your option) any later version.
%
%     This program is distributed in the hope that it will be useful, but
%     WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [r, pbp]=InitialConditions_Diffusion(basicinput,data)
% Inputs: basicinput=Basic user defined input from the main GUI; data=image
% data from the current stack
classname=class(data{1,1}{1,1});
for index1=1:basicinput{1,8} % Load prebleach images into variable img
    index2=basicinput{1,6}-1-basicinput{1,8}:basicinput{1,6}-1;
    img(:,:,index1)=double(data{1,1}{index2(index1),1});
end
img2=mean(img,3); % Load mean prebleach image into variable img2
pbi=double(data{1,1}{basicinput{1,6},1}); % Load postbleach image into variable pbi
normpbi=pbi./img2; % Normalize the postbleach image
normpbi(normpbi>2)=2; % Suppress extreme values to enhance visualization
width=size(normpbi,2); % Determine the width of postbleach image
height=size(normpbi,1); % Determine the height of postbleach image
originX=basicinput{1,7}(1); % Load the X coordinates for the center of the bleaching ROI
originY=basicinput{1,7}(2); % Load the Y coordinates for the center of the bleaching ROI
[x,y]=meshgrid(1:width,1:height);
r=sqrt((x-originX).^2+(y-originY).^2); % radial distances from the center of the bleaching ROI
r(pbi==intmax(classname))=NaN; % Over-saturated pixels are not useful
r(pbi==0)=NaN; % Under-saturated pixels are not useful
r(img2==intmax(classname))=NaN; % Over-saturated pixels are not useful
r(img2==0)=NaN; % Under-saturated pixels are not useful
A=[r(~isnan(r)),normpbi(~isnan(r))]; % Average intensities at the same radial distances
[r,~,id2] = unique(A(:,1),'rows');
pbp = accumarray(id2,A(:,2))./accumarray(id2,1);

end