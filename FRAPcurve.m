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

function [frap time cell adjacent]=FRAPcurve(data,basicinput,omeMeta,bleachroimask,name,cellroimask,adjacentroimask)
% Inputs: data=image planes from image stack; basicinput=basic user
% inputs from the main GUI; metadata=meta data from image stack;
% bleachroimask=mask for bleaching; name=name of stack file;
% region; cellroimask=mask for cell region.  Outputs: frap=mean intensity
% inside the bleach region; time=time information from the image stack;
% cell=mean intensity inside the cell region.

for index1=1:length(data{1,1})
    
    img=data{1,1}{index1,1}; % Load the current image plane in the stack
    
    frap(index1)=mean(img(bleachroimask))-basicinput{1,5}; % Load the mean intensity inside bleach ROI to variable frap
    
    if basicinput{1,9}==1;
        adjacent(index1)=mean(img(adjacentroimask))-basicinput{1,5}; % Load the mean intensity inside adjacent ROI to variable adjacent
    else
        adjacent=0;
    end
    %     time(index1)=metadata.get([name,' TimeStamp',num2str(index1-1)]); % Load the time information to variable time
    time(index1)=double(omeMeta.getPlaneDeltaT(0,index1-1));
    if basicinput{1,4}==2 % IF normalizing by whole cell intensity
        cell(index1)=mean(img(cellroimask)); % Load the mean intensity inside the cell ROI to variable cell
    else
        cell=0;
    end
    
    
    
end

end