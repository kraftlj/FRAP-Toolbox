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

function [decayrate, data]=PhotoDecay_Reaction(data,basicinput,usrinputs,val)
% Input:
% data - output from loadData_Reaction
% basicinput - Basic user inputs from the Main_GUI
% usrinputs - Basic user inputs from the Figure_GUI_Reaction
% val - The datasets that were selected by the user for plotting.

% Output:
% decayrate - This is a decay rate obtained by fitting some portion of the end of the frapdata to a single exponential
% data - This is the data which was inputed amended with the FRAP curves
% corrected for photodecay.

options=optimset('lsqcurvefit');
options.Display='off';

%% Error checking
% Check to make sure time is not different across data sets. This is the
% primary way that the program ensures the FRAP datasets are identical when
% doing batch processing of multiple datasets at once.
for index1=1:length(val)
    t(:,index1)=data(val(index1)).time-data(val(index1)).time(basicinput{1,6})';
    t(:,index1)=t(:,index1)-t(1,index1);
    f(:,index1)=data(val(index1)).normfrap;
end
if size(t,2)>1
    isequalRel = @(x,y,tol) ( abs(x-y) <= ( tol*max(abs(x),abs(y)) + eps) );
    for index1=1:length(val)-1
        x=isequalRel(t(:,1),t(:,index1+1),1e-1);
        a=find(~x);
        if isempty(a)
            b(index1)=NaN;
        else
            b(index1)=index1;
        end
    end
    if sum(isnan(b))~=length(val)-1
        h = errordlg(['Time vectors 1 and ',num2str(b(~isnan(b))+1),' are not the same, check to make sure the frap data was acquired using the same settings']);
        %     dec
        return
    end
    %--------------------------------------------------------------------------
    %% Find the decay constant of the averaged FRAP datasets
    % Calculate the mean of all the equivalent FRAP datasets (it easier to find
    % the correct decay constant by fitting the average if multiple datasets)
    f=mean(f,2);
end

fun=@(p,t) p(1)*exp(-p(2)*t); % This is the single exponential fitting function
p0=[.9,usrinputs{4,1}]; % Define the initial fitting parameters
% solve the non-linear least squares problem
p=lsqcurvefit(fun,p0,t(usrinputs{6,1}:usrinputs{6,2},1),f(usrinputs{6,1}:usrinputs{6,2}),[0,usrinputs{4,2}],[2,usrinputs{4,3}],options);
decayrate=p(2);
for index1=1:length(val)
    data(val(index1)).correctfrap=data(val(index1)).normfrap./exp(-decayrate*t(:,1)'); % Correct the data using the optimized decayrate
end
end