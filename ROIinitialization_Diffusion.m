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

function [bleachroimask cellroimask adjacentroimask]=ROIinitialization_Diffusion(img,basicinput)
% Inputs: 
% img - image used for user defining ROIs.
% basicinput - basic user inputs from the main GUI.  

% Outputs: 
% bleachroimask - mask for bleaching region.
% cellroimask - mask for cell region.
% adjacentroimask - mask for the ROI adjacent to the bleaching region (only
% if correcting the mobile fraction).
   
    x0=basicinput{1,7}(1); % X coordinates for the center of a circular bleaching ROI
    y0=basicinput{1,7}(2); % Y coordinates for the center of a circular bleaching ROI
    R0=basicinput{1,7}(3); % Radius of circular bleaching ROI
    t = 0:pi/20:2*pi;
    xi = R0*cos(t)+x0; % X coordinates for the circle
    yi = R0*sin(t)+y0; % Y coordinates for the circle
    bleachroimask = poly2mask(xi,yi, size(img,1),size(img,2)); % Create a mask using the circular coordinates
    
    if basicinput{1,9}==1; % IF using an adjacent ROI to calculate a corrected mobile fraction
    xi = R0*cos(t)+x0+R0*2.5; % X coordinates for the circle
    yi = R0*sin(t)+y0; % Y coordinates for the circle
    adjacentroimask = poly2mask(xi,yi, size(img,1),size(img,2));
    else
        adjacentroimask = NaN;
    end
    
    if basicinput{1,4}==2 % IF normalizing by whole cell intensity
        h=msgbox({'Define the cell region using the mouse.  Right click and press create mask when finished.'});
        movegui(h,'northwest');
        cellroimask = roipoly(imadjust(img,stretchlim(img,[0,1]))); % Create a mask using the user defined ROI
        h2=gcf;
        close(h);
        close(h2);
    else
        cellroimask = NaN;
    end

end