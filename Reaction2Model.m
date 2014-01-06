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

function [data]=Reaction2Model(data, basicinput)
if basicinput{1,13}==1
    for index1=1:length(data)
        avgf(index1,:)=data(index1).correctfrap;
    end
    f=mean(avgf,1);
    data(end+1).correctfrap=f;
    f=f(basicinput{1,6}:basicinput{1,11}+basicinput{1,6});
    t=data(1).time;
    data(end).time=t;
    t=t(basicinput{1,6}:basicinput{1,11}+basicinput{1,6})-data(1).time(basicinput{1,6});
    data(end).t=t;
    wf=f./(t+sum(f));
    frapfun=@(p,t) ((1-p(1)*exp(p(2)*t)-(1-p(1)-f(1))*exp(p(3)*t)).*p(4)+(1-p(4))*f(1))./(t+sum(f));
    p=lsqcurvefit(frapfun,[.5,-1,-1,1],t,wf,[0,-Inf,-Inf,0],[f(1),0,0,Inf]);
    data(end).ceq1=p(1);
    data(end).k1=p(2);
    data(end).ceq2=(1-p(1)-f(1));
    data(end).k2=p(3);
    data(end).MF=p(4);
    frapfit=frapfun(p,t).*(t+sum(f));
    data(end).frapfit=frapfit;
    data(end).frapres=f-frapfit;
    data(end).r2=1-(sum(data(end).frapres.^2)/((length(f)-1)*var(f)))*((length(f)-1)/(length(f)-length(p)-1));
    
else
    for index1=1:length(data)
        avgf(index1,:)=data(index1).correctfrap;
        t=data(index1).time(basicinput{1,6}:basicinput{1,11}+basicinput{1,6})-data(index1).time(basicinput{1,6});
        f=data(index1).correctfrap(basicinput{1,6}:basicinput{1,11}+basicinput{1,6});
        data(index1).t=t;
        wf=f./(t+sum(f));
        frapfun=@(p,t) ((1-p(1)*exp(p(2)*t)-(1-p(1)-f(1))*exp(p(3)*t)).*p(4)+(1-p(4))*f(1))./(t+sum(f));
        p=lsqcurvefit(frapfun,[.5,-1,-1,1],t,wf,[0,-Inf,-Inf,0],[f(1),0,0,Inf]);
        data(index1).ceq1=p(1);
        data(index1).k1=p(2);
        data(index1).ceq2=(1-p(1)-f(1));
        data(index1).k2=p(3);
        data(index1).MF=p(4);
        frapfit=frapfun(p,t).*(t+sum(f));
        data(index1).frapfit=frapfit;
        data(index1).frapres=f-frapfit;
        data(index1).r2=1-(sum(data(index1).frapres.^2)/((length(f)-1)*var(f)))*((length(f)-1)/(length(f)-length(p)-1));
    end
    f=mean(avgf,1);
    data(end+1).correctfrap=f;
    f=f(basicinput{1,6}:basicinput{1,11}+basicinput{1,6});
    t=data(1).time;
    data(end).time=t;
    t=t(basicinput{1,6}:basicinput{1,11}+basicinput{1,6})-data(1).time(basicinput{1,6});
    data(end).t=t;
    ceq1=mean([data.ceq1]);
    data(end).ceq1=ceq1;
    k1=mean([data.k1]);
    data(end).k1=k1;
    ceq2=mean([data.ceq2]);
    data(end).ceq2=ceq2;
    k2=mean([data.k2]);
    data(end).k2=k2;
    MF=mean([data.MF]);
    data(end).MF=MF;
    frapfun=@(p,t) ((1-p(1)*exp(p(2)*t)-ceq2*exp(p(3)*t)).*p(4)+(1-p(4))*f(1))./(t+sum(f));
    frapfit=frapfun([ceq1,k1,k2,MF],t).*(t+sum(f));
    data(end).frapfit=frapfit;
    data(end).frapres=f-frapfit;
    data(end).r2=1-(sum(data(end).frapres.^2)/((length(f)-1)*var(f)))*((length(f)-1)/(length(f)-length(p)-1));
end

end
