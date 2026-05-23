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

function [frap, time, cell]=FRAPcurve_Reaction(data,basicinput,omeMeta,bleachroimask,cellroimask)
% Inputs:
% data - image planes from FRAP dataset.
% basicinput - basic user inputs from the main GUI.
% metadata - meta data from FRAP dataset.
% bleachroimask - mask for the bleaching region.
% name - name of FRAP dataset file.
% cellroimask - mask for cell region.


% Outputs:
% frap - mean intensity inside the bleach region.
% time - time information from the FRAP dataset.
% cell - mean intensity inside the cell region.


for index1=1:length(data{1,1})
    
    img=data{1,1}{index1,1}; % Load the current image plane in the stack
    
    frap(index1)=mean(img(bleachroimask))-basicinput{1,5}; % Load the mean intensity inside bleach ROI to variable frap; subtract the background
        
    time(index1)=double(omeMeta.getPlaneDeltaT(0,index1-1)); % Load the time information from the meta data using Bio-Formats
    
    if basicinput{1,4}==2 % IF normalizing by whole cell intensity
        cell(index1)=mean(img(cellroimask)); % Load the mean intensity inside the cell ROI to variable cell
    else
        cell=0;
    end
    
    
    
end

end