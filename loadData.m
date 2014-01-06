%%% FRAP Toolbox
%%
%     FRAP Toolbox is designed to be a modular software program
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

function [data]=loadData(fullfilepaths,basicinput)
% Inputs: basicinput=basic user input from main GUI; fullfilepaths=cell
% array of the full paths to user selected files in the main GUI.
% Outputs:
% Function initialization--------------------------------------------------
h=waitbar(0/length(fullfilepaths),['Loading Data (Initializing)']);
for index1=1:length(fullfilepaths) % Loop through each file one by one
    waitbar((index1-1)/length(fullfilepaths),h,['Loading Data (% Complete): ', num2str((index1-1)/length(fullfilepaths)*100)])
    
    [pathstr, name, ext] = fileparts(fullfilepaths{index1}); % define the file parts
    
    imgdata=bfopen(fullfilepaths{index1}); % Load the image file using open Bio-Formats
    
%     metadata=imgdata{1,2}; % Load the meta data information into variable metadata
    omeMeta=imgdata{1,4};
    img=imgdata{1,1}{1,1}; % Load the first image plane in the stack for use in user defining ROI
    % End of Function initialization-------------------------------------------
    
    % ROI initialization-------------------------------------------------------
    [bleachroimask cellroimask adjacentroimask]=ROIinitialization(img,basicinput);
    % End of ROI intialization-------------------------------------------------
    
    % Record FRAP curve information--------------------------------------------
    [data(index1).frap data(index1).time data(index1).cell data(index1).adjacent]=FRAPcurve(imgdata,basicinput,omeMeta,bleachroimask,name,cellroimask,adjacentroimask);
    % End of Record FRAP curve information-------------------------------------
    
    % Normalize the FRAP curve-------------------------------------------------
    [data(index1).normfrap data(index1).correctMF data(index1).adjacent]=NormalizeFRAP(data(index1).frap,data(index1).cell,basicinput,data(index1).adjacent);
    % End of Normalize the FRAP curve------------------------------------------
    
    % Get post-bleach profile information
    % We only need this function if we are using the 2D diffusion model from
    % Kang et al. 2008
    data(index1).voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).getValue();
    if basicinput{1,2}==1
        [data(index1).r data(index1).pbp]=InitialConditions(basicinput,imgdata);
        data(index1).r=data(index1).r.*data(index1).voxelSizeX;
    else
        data(index1).r=0;
        data(index1).pbp=0;
    end
    % End of Get post-bleach profile information
end
close(h)
end
