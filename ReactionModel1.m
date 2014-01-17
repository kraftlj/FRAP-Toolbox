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

function [data, avg]=ReactionModel1(data, basicinput, usrinputs, val)
% Inputs:
% data - this is the data output from PhotoDecay_Reaction
% basicinput - this is basic user input from the Main_GUI
% usrinputs - this is basic user input from Figure_GUI_Reaction
% val - the datasets that were selected by the user for plotting

options=optimset('lsqnonlin');
options.Display='off';
%% Fit the averaged data;
if usrinputs{7,1}==2
    for index1=1:length(val)
        avgf(index1,:)=data(val(index1)).correctfrap;
    end
    f=mean(avgf,1);
    avg.f=f;
    f=f(basicinput{1,6}+usrinputs{5,1}-1:usrinputs{5,2});
    t=data(1).time;
    avg.time=t;
    t=t(basicinput{1,6}+usrinputs{5,1}-1:usrinputs{5,2})-data(1).time(basicinput{1,6}+usrinputs{5,1}-1);
    avg.t=t+t(usrinputs{5,1});
    wf=f./(t+sum(f));
    %----------------------------------------------------------------------
    
    var_indx=strcmp('Fixed',usrinputs(1:3,4));
    fun = @(p) Reaction1(p,[usrinputs{1:3,1}],var_indx,t,f);
    p=lsqnonlin(fun,[usrinputs{1:3,1}],[usrinputs{1:3,2}],[usrinputs{1:3,3}],options);
    avg.a=p(1);
    avg.b=p(2);
    avg.c=p(3);
    
    frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)) ./(t+sum(f));
    frapfit=frapfun([avg.a,avg.b,avg.c],t).*(t+sum(f)); % unweight the optimized fit
    avg.frapfit=frapfit;
    avg.frapres=wf-frapfun([avg.a,avg.b,avg.c],t); % weighted residuals
    avg.SS=sum(avg.frapres.^2); % Sum of square errors
    %----------------------------------------------------------------------
    
    %% Fit the individual FRAP datasets
else
    for index1=1:length(val)
        avgf(index1,:)=data(val(index1)).correctfrap;
        t=data(val(index1)).time(basicinput{1,6}+usrinputs{5,1}-1:usrinputs{5,2})-data(val(index1)).time(basicinput{1,6}+usrinputs{5,1}-1);
        f=data(val(index1)).correctfrap(basicinput{1,6}+usrinputs{5,1}-1:usrinputs{5,2});
        data(val(index1)).t=t+t(usrinputs{5,1});
        wf=f./(t+sum(f)); % weight the beginning of the FRAP curve more
        
        var_indx=strcmp('Fixed',usrinputs(1:3,4));
        fun = @(p) Reaction1(p,[usrinputs{1:3,1}],var_indx,t,f);
        p=lsqnonlin(fun,[usrinputs{1:3,1}],[usrinputs{1:3,2}],[usrinputs{1:3,3}],options);
        data(val(index1)).a=p(1);
        data(val(index1)).b=p(2);
        data(val(index1)).c=p(3);
        
        frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)) ./(t+sum(f));
        frapfit=frapfun([data(val(index1)).a,data(val(index1)).b,data(val(index1)).c],t).*(t+sum(f)); %unweight the optimized FRAP fit
        data(val(index1)).frapfit=frapfit;
        data(val(index1)).frapres=wf-frapfun([data(val(index1)).a,data(val(index1)).b,data(val(index1)).c],t); % weighted residuals
        data(val(index1)).SS=sum(data(val(index1)).frapres.^2); % Sum of squared errors
    end
    
    f=mean(avgf,1);
    avg.f=f;
    f=f(basicinput{1,6}+usrinputs{5,1}-1:usrinputs{5,2});
    t=data(1).time;
    avg.time=t;
    t=t(basicinput{1,6}+usrinputs{5,1}-1:usrinputs{5,2})-data(1).time(basicinput{1,6}+usrinputs{5,1}-1);
    wf=f./(t+sum(f));
    avg.t=t+t(usrinputs{5,1});
    a=mean([data.a]);
    avg.a=a;
    b=mean([data.b]);
    avg.b=b;
    c=mean([data.c]);
    avg.c=c;
    frapfun=@(p,t) (p(1)-p(2)*exp(-p(3)*t)) ./(t+sum(f));
    frapfit=frapfun([a,b,c],t).*(t+sum(f)); % unweighted optimized fit using the averaged parameters
    avg.frapfit=frapfit;
    avg.frapres=wf-frapfun([a,b,c],t); % Weighted residuals
    avg.SS=sum(avg.frapres.^2); % Sum of squared errors
end

%% This is the theoretical Reaction 1 FRAP function

    function out = Reaction1(p_all,p_var,var_indx,t,f)
        
        p = p_all;
        p(var_indx) = p_var(var_indx);
        F = p(1)-p(2)*exp(-p(3)*t);
        out = (F-f)./(t+sum(f));
        
    end


end
