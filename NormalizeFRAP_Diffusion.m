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

function [normfrap, correctMF, normadjacent]=NormalizeFRAP_Diffusion(frap,cell,basicinput,adjacent)
% Inputs:
% frap - mean intensity inside bleach region.
% cell - mean intensity inside cell region.
% basicinput - basic user inputs from the main GUI.

% Outputs:
% normfrap - normalized FRAP curve.
% correctMF - the mobile fraction corrected for loss of fluorescence in the
% compartment due to the bleach.
% normadjacent - normalized intensity inside the adjacent ROI

if basicinput{1,4}==2 % IF normalizing by whole cell intensity
    normfrap=frap./cell; % This takes into account loss of fluorescence in the compartment and photodecay.
    normfrap=normfrap./mean(normfrap(basicinput{1,6}-1-basicinput{1,8}+1:basicinput{1,6}-1)); % This normalizes to the pre-bleach intensity
    
    correctMF=NaN; % no reason to correct the MF.
    normadjacent=NaN;
    
    
else % IF normalizing by the pre-bleach steady-state intensity
    normfrap=frap./mean(frap(basicinput{1,6}-1-basicinput{1,8}+1:basicinput{1,6}-1)); % normalize to the pre-bleach intensity.
    
    
    if basicinput{1,9}==1; % If using an adjacent ROI to calculate a corrected Mobile fraction
        normadjacent=adjacent./mean(adjacent(basicinput{1,6}-1-basicinput{1,8}+1:basicinput{1,6}-1));
        correctMF=1-(mean(normadjacent(end-basicinput{1,8}*3:end))-mean(normfrap(end-basicinput{1,8}*3:end)));
    else
        correctMF=NaN;
        normadjacent=NaN;
    end
    
end

end