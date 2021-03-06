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

function [data]=loadData_Diffusion(fullfilepaths,basicinput)
% Inputs:
% basicinput - basic user input from main GUI
% fullfilepaths - cell array of the full paths to user selected files in the main GUI.

% Outputs:
% data - the normalized FRAP curves and intial conditions

% Function initialization--------------------------------------------------
h=waitbar(0/length(fullfilepaths),['Loading Data (Initializing)']);
for index1=1:length(fullfilepaths) % Loop through each file one by one
    waitbar((index1-1)/length(fullfilepaths),h,['Loading Data: ', num2str((index1-1)/length(fullfilepaths)*100),'%'])
    
    [pathstr, name, ext] = fileparts(fullfilepaths{index1}); % define the file parts
    
    try
        imgdata=bfopen(fullfilepaths{index1}); % Load the image file using open Bio-Formats
    catch errObj
        errordlg('Bioformats does not support the filetype you selected.  Please select a microscopy image.');
        close(h);
    end
    
    omeMeta=imgdata{1,4}; % Load the meta data information from the raw microscope file
    img=imgdata{1,1}{1,1}; % Load the first image plane in the stack for use in user defining ROI
    % End of Function initialization-------------------------------------------
    
    % ROI initialization-------------------------------------------------------
    [bleachroimask, cellroimask, adjacentroimask]=ROIinitialization_Diffusion(img,basicinput);
    % End of ROI intialization-------------------------------------------------
    
    % Record FRAP curve information--------------------------------------------
    [data(index1).frap, data(index1).time, data(index1).cell, data(index1).adjacent]=FRAPcurve_Diffusion(imgdata,basicinput,omeMeta,bleachroimask,cellroimask,adjacentroimask);
    % End of Record FRAP curve information-------------------------------------
    
    % Normalize the FRAP curve-------------------------------------------------
    [data(index1).normfrap, data(index1).correctMF, data(index1).adjacent]=NormalizeFRAP_Diffusion(data(index1).frap,data(index1).cell,basicinput,data(index1).adjacent);
    % End of Normalize the FRAP curve------------------------------------------
    
    % Get post-bleach profile information----------------------------------
    data(index1).voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).getValue(); % Read the pixel size from the image metadata.
    [data(index1).r data(index1).pbp]=InitialConditions_Diffusion(basicinput,imgdata); % Calculate the normalize post-bleach profile.
    data(index1).r=data(index1).r.*data(index1).voxelSizeX; % convert units from pixels.
    % End post-bleach profile information----------------------------------
    
end
close(h)
end
