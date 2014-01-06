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

function [data avg]=NCtransportModel2(data, basicinput, usrinputs, val)
options=optimset('lsqcurvefit');
options.Display='off';
if usrinputs{5,1}==2 %Fit the averaged data;
    for index1=1:length(val)
        avgf(index1,:)=data(val(index1)).correctfrap;
    end
    f=mean(avgf,1);
    avg.f=f;
    f=f(basicinput{1,6}+usrinputs{4,1}-1:usrinputs{4,2});
    t=data(1).time;
    avg.time=t;
    t=t(basicinput{1,6}+usrinputs{4,1}-1:usrinputs{4,2})-data(1).time(basicinput{1,6}+usrinputs{4,1}-1);
    avg.t=t+t(usrinputs{4,1});
    %----------------------------------------------------------------------
    wf=f./(t+sum(f));
    switch [usrinputs{1,4},usrinputs{2,4}]
        case 'AdjustableAdjustable'
            frapfun=@(p,t) ((p(3)-(p(3)-f(1))*exp(-p(1)*t)).*exp(-p(2)*t))  ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1},usrinputs{3,1}],t,wf,[usrinputs{1,2},usrinputs{2,2},usrinputs{3,2}],[usrinputs{1,3},usrinputs{2,3},usrinputs{3,3}],options);
            avg.k=p(1);
            avg.kdecay=p(2);
            avg.asymp=p(3);
        case 'FixedAdjustable'
            avg.k=usrinputs{1,1};
            frapfun=@(p,t) ((1-(1-f(1))*exp(-avg.k*t)).*exp(-p(1)*t))  ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{2,1}],t,wf,[usrinputs{2,2}],[usrinputs{2,3}],options);
            avg.kdecay=p(1);
        case 'AdjustableFixed'
            avg.kdecay=usrinputs{2,1};
            frapfun=@(p,t) ((1-(1-f(1))*exp(-p(1)*t)).*exp(-avg.kdecay*t))  ./(t+sum(f));
            p=lsqcurvefit(frapfun,[usrinputs{1,1}],t,wf,[usrinputs{1,2}],[usrinputs{1,3}],options);
            avg.k=p(1);
        case 'FixedFixed'
            avg.kdecay=usrinputs{2,1};
            avg.k=usrinputs{1,1};
    end
    frapfun=@(p,t) ((p(3)-(p(3)-f(1))*exp(-p(1)*t)).*exp(-p(2)*t))  ./(t+sum(f));
    frapfit=frapfun([avg.k,avg.kdecay,avg.asymp],t).*(t+sum(f));
    avg.frapfit=frapfit;
    avg.frapres=wf-frapfun([avg.k,avg.kdecay,avg.asymp],t);
    avg.SS=sum(avg.frapres.^2);
    %----------------------------------------------------------------------
else %Fit the individual curves
    for index1=1:length(val)
        avgf(index1,:)=data(val(index1)).correctfrap;
        t=data(val(index1)).time(basicinput{1,6}+usrinputs{3,1}-1:usrinputs{3,2})-data(val(index1)).time(basicinput{1,6}+usrinputs{3,1}-1);
        f=data(val(index1)).correctfrap(basicinput{1,6}+usrinputs{3,1}-1:usrinputs{3,2});
        data(val(index1)).t=t+t(usrinputs{3,1});
 
        wf=f./(t+sum(f));
        switch [usrinputs{1,4},usrinputs{2,4}]
            case 'AdjustableAdjustable'
                frapfun=@(p,t) ((1-(1-f(1))*exp(-p(1)*t)).*exp(p(2)*t))  ./(t+sum(f));
                p=lsqcurvefit(frapfun,[usrinputs{1,1},usrinputs{2,1}],t,wf,[usrinputs{1,2},usrinputs{2,2}],[usrinputs{1,3},usrinputs{2,3}],options);
                data(val(index1)).k=p(1);
                data(val(index1)).kdecay=p(2);
            case 'FixedAdjustable'
                data(val(index1)).k=usrinputs{1,1};
                frapfun=@(p,t) ((1-(1-f(1))*exp(-usrinputs{1,1}*t)).*exp(-p(1)*t))  ./(t+sum(f));
                p=lsqcurvefit(frapfun,[usrinputs{2,1}],t,wf,[usrinputs{2,2}],[usrinputs{2,3}],options);
                data(val(index1)).kdecay=p(1);
            case 'AdjustableFixed'
                data(val(index1)).kdecay=usrinputs{2,1};
                frapfun=@(p,t) ((1-(1-f(1))*exp(-p(1)*t)).*exp(-usrinputs{2,1}*t))  ./(t+sum(f));
                p=lsqcurvefit(frapfun,[usrinputs{1,1}],t,wf,[usrinputs{1,2}],[usrinputs{1,3}],options);
                data(val(index1)).k=p(1);
            case 'FixedFixed'
                data(val(index1)).k=usrinputs{1,1};
                data(val(index1)).kdecay=usrinputs{2,1};
        end
        frapfun=@(p,t) ((1-(1-f(1))*exp(-p(1)*t)).*exp(-p(2)*t))  ./(t+sum(f));
        frapfit=frapfun([data(val(index1)).k,data(val(index1)).kdecay],t).*(t+sum(f));
        data(val(index1)).frapfit=frapfit;
        data(val(index1)).frapres=wf-frapfun([data(val(index1)).k,data(val(index1)).kdecay],t);
        data(val(index1)).SS=sum(data(val(index1)).frapres.^2);
    end
    
    f=mean(avgf,1);
    avg.f=f;
    f=f(basicinput{1,6}+usrinputs{3,1}-1:usrinputs{3,2});
    t=data(1).time;
    avg.time=t;
    t=t(basicinput{1,6}+usrinputs{3,1}-1:usrinputs{3,2})-data(1).time(basicinput{1,6}+usrinputs{3,1}-1);
    avg.t=t+t(usrinputs{3,1});
    k=mean([data.k]);
    avg.k=k;
    kdecay=mean([data.kdecay]);
    avg.kdecay=kdecay;
    wf=f./(t+sum(f));
    frapfun=@(p,t) ((1-(1-f(1))*exp(-p(1)*t)).*exp(-p(2)*t))  ./(t+sum(f));
    frapfit=frapfun([k,kdecay],t).*(t+sum(f));
    avg.frapfit=frapfit;
    avg.frapres=wf-frapfun([k,kdecay],t);
    avg.SS=sum(avg.frapres.^2);
end

end
