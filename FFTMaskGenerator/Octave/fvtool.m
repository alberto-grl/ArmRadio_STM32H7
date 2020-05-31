function []=fvtool(varargin)
% function []=fvtool(varargin)
% Input:
% - fvtool(h) - plotting frequency response of FIR filter having impulse response 'h'
% - fvtool(num,den) - plotting frequency response of IIR filter having numerator 'num' and denominator 'den'
% - fvtool(num1,den1,num2,den2) - plotting several IIR frequency responses
% - setting logarithmic frequency scale (default is 'Linear' ('Lin')
%     fvtool(num,den,'FreqScale','Logarithmic'), fvtool(num,den,'Freq','Log')
% - setting logarithmic magnitude axis (default is 'Linear' ('Lin'))
%     fvtool(num,den,'MagScale','Logarithmic'), fvtool(num,den,'Mag','Log'), fvtool(num,den,'Magnitude','dB')
% - setting Sampling Rate for X axis (default is 2 - normalized X axis)
%     fvtool(num,den,'Fs',samplerate)
% - setting x axis boundaries (y bounaries set automatically)
%     fvtool(num,den,'AxisX',[xmin xmax])
% - setting y axis boundaries (x bounaries set automatically)
%     fvtool(num,den,'AxisY',[ymin ymax])
% - setting x and y axis boundaries simultaneously
%     fvtool(num,den,'AxisXY',[xmin xmax ymin ymax]),   fvtool(num,den,'Axes',[xmin xmax ymin ymax])

%% Init
FreqscaleLog_flag=0;   % linear default
MagscaleLog_flag=0;   % linear default
SampleRate=-1;         % ie 1 is default Nyquist freq
Fs_flag=0;
Xbnd=[]; Ybnd=[];
%clrord='brgcmygk';    clrlen=length(clrord);
%clridx=0;

lowbnd=180;
highbnd=-180;

Disp_Mag=0;
Disp_Phase=1;
Disp_PhaseDelay=2;
Disp_GroupDelay=3;
Disp_ImpulseResponse=4;

disp_flag= Disp_Mag;

% parse parameters
k=0;i=1;
while (i<=nargin),
  par=varargin{i};
  if ischar(par),
    if ~isempty(strfind('FreqScale',par)),
      i=i+1;
      if (i<=nargin),
        par=varargin{i};
        if ~isempty( strfind ('Logarithmic',par)),
          FreqscaleLog_flag=1;
        end
      end
    elseif ~isempty( strfind ('MagScale',par)) || ~isempty( strfind ('Magnitude',par)),
      i=i+1;
      if (i<=nargin),
        par=varargin{i};
        if ~isempty( strfind ('Logarithmic',par)) || ~isempty( strfind ('dB',par)),
          MagscaleLog_flag=1;
        end
      end
    elseif ~isempty( strfind ('Display',par)),
      i=i+1;
      if (i<=nargin),
        par=varargin{i};
        if ~isempty( strfind ('Phase',par)),        disp_flag = Disp_Phase ;
        elseif ~isempty( strfind ('PhaseDelay',par)),        disp_flag = Disp_PhaseDelay ;
        elseif ~isempty( strfind ('GroupDelay',par)),        disp_flag = Disp_GroupDelay ;
        elseif ~isempty( strfind (' ImpulseResponse ',par)),        disp_flag = Disp_ImpulseResponse ;
        else disp_flag = Disp_Mag ;
        end
      end
    elseif ~isempty( strfind ('Fs',par)),
      i=i+1;
      if (i<=nargin),
        par=varargin{i};
        if ~ischar(par),
          if ~isempty(par),
            SampleRate=par;
            Fs_flag=1;
          end
        end
      end
    elseif ~isempty( strfind ('AxisX',par)),
      i=i+1;
      if (i<=nargin),
        par=varargin{i};
        if ~ischar(par),
          if length(par)==2,
            Xbnd=par;
          end
        end
      end
    elseif ~isempty( strfind ('AxisY',par)),
      i=i+1;
      if (i<=nargin),
        par=varargin{i};
        if ~ischar(par),
          if length(par)==2,
            Ybnd=par;
          end
        end
      end
    elseif (~isempty( strfind ('AxisXY',par))) || (~isempty( strfind ('Axes',par))),
      i=i+1;
      if (i<=nargin),
        par=varargin{i};
        if ~ischar(par),
          if length(par)==4,
            Xbnd=par(1:2); Ybnd=par(3:4);
          end
        end
      end
    end
  elseif rem(k,2)==0,
    num{floor(k/2)+1}=par; k=k+1;
  else
    den{floor(k/2)+1}=par; k=k+1;
  end
  i=i+1;
end

%% fill last denominator (if not present)
if  rem(k,2),
  den{floor(k/2)+1}=1;
end

%% preparefrequency grid
if disp_flag==Disp_ImpulseResponse,	sr_dflt=1;	% default samplerate
else								sr_dflt=2;
end
if SampleRate<=0,
	SampleRate=sr_dflt;
end
Fn=SampleRate/2;

Nf=500; % 8192;
if FreqscaleLog_flag,
  f=logspace(log10(pi/1000),log10(pi),Nf);
%  cm='semilogx(';
  xmin=f(1)/pi*Fn;   xmax=f(end)/pi*Fn;
else
  f=linspace(0,pi,Nf);
%  cm='plot(';
  xmin=0;   xmax=Fn;
end

%% compute complex frequency response
if disp_flag~=Disp_ImpulseResponse,
  s=zeros(length(f),length(num));
  for i=1:length(num),
    s(:,i)=freqz( num{i} , den{i} ,f);
  end
end

%% compute y-axis data according to disp_flag
switch disp_flag,

case Disp_Mag,
  s=abs(s);
  if MagscaleLog_flag, s=20*log10(s); end
  for i=1:length(num),
%    cm=[cm 'f/pi*Fn,s(:,' num2str(i) '),''' clrord(clridx+1) ''','];
%    clridx=rem(clridx+1,clrlen);
    t=s(:,i);
    t=t(~isnan(t));
    if MagscaleLog_flag,
      t=sort(t);    l=t(round(2*length(t)/100));
    else
      l=min(t);
    end
    lowbnd=min(lowbnd,l);
    h=max(t);
    highbnd=max(highbnd,h);
  end

  plot1( f/pi*Fn, s,FreqscaleLog_flag);

  if MagscaleLog_flag,
    highbnd=highbnd+3;
  else
    highbnd=highbnd*1.1;
  end

  xlabelstr=['Frequency (' getFreqAxisType(FreqscaleLog_flag) ')'];
  ylabelstr=['Magnitude (' getMagAxisType(MagscaleLog_flag) ')'];
  titlestr=['Frequency response'];

case Disp_Phase,
  for i=1:length(num),
    s(:,i)=phase(s(:,i));
%    cm=[cm 'f/pi*Fn,s(:,' num2str(i) '),''' clrord(clridx+1) ''','];
%    clridx=rem(clridx+1,clrlen);

    t=s(:,i);
    t=t(~isnan(t));
    l=min(t);
    lowbnd=min(lowbnd,l);
    h=max(t);
    highbnd=max(highbnd,h);
  end

  plot1( f/pi*Fn, s,FreqscaleLog_flag);

  xlabelstr=['Frequency (' getFreqAxisType(FreqscaleLog_flag) ')'];
  ylabelstr=['Phase (rad)'];
  titlestr=['Phase response'];

case Disp_PhaseDelay,
  for i=1:length(num),
    t=-phase1(s(:,i))./f';
    if isnan(t(1)), t(1)=t(2); end
    if Fs_flag, t=t/SampleRate; end
    s(:,i)=t;
    t=t(~isnan(t));
    l=min(t);
    lowbnd=min(lowbnd,l);
    h=max(t);
    highbnd=max(highbnd,h);
  end

  plot1( f/pi*Fn, s,FreqscaleLog_flag);

  xlabelstr=['Frequency (' getFreqAxisType(FreqscaleLog_flag) ')'];
  ylabelstr=['Phase delay (' getTimeAxisType(SampleRate,sr_dflt) ')'];
  titlestr=['Phase delay response'];

case Disp_GroupDelay,
  for i=1:length(num),
    t=-phase1(s(:,i));
    t=diff(t)./diff(f');
    t=[t(1); t];
    if Fs_flag, t=t/SampleRate; end
    s(:,i)=t;
    t=t(~isnan(t));
    l=min(t);
    lowbnd=min(lowbnd,l);
    h=max(t);
    highbnd=max(highbnd,h);
  end

  plot1( f/pi*Fn, s,FreqscaleLog_flag);

  xlabelstr=['Frequency (' getFreqAxisType(FreqscaleLog_flag) ')'];
  ylabelstr=['Group delay (' getTimeAxisType(SampleRate,sr_dflt) ')'];
  titlestr=['Group delay response'];

case Disp_ImpulseResponse,
  N= 8192;
  maxidx=0;
  y=zeros(N,length (num));
  for i=1:length(num),
    t=filter (num{i},den {i}, [1 zeros(1,N-1)]);
    y(:,i)=t';     m=max (t);
    idx=find(t >0.005*m);
    maxidx=max (maxidx,idx (end));
  end
  y=y (1:maxidx,:);

  %plot(kron((0:size(y,1)-1)'/SampleRate,ones(1,size(y,2))),y);
  x=(0:size(y,1)-1)/SampleRate;
  plot(x,y);

  xlabelstr=['Time (' getTimeAxisType(SampleRate,sr_dflt) ')'];
  ylabelstr=['Magnitude'];
  titlestr=['Impulse response'];

  xmin=x(1); xmax=x(end);
  lowbnd=min(min(y)); highbnd=max(max(y));
end

%% replace automatic boundaries by manual ones (if any)
if ~isempty(Xbnd),
  xmin= Xbnd(1); xmax= Xbnd(2);
end
if ~isempty(Ybnd),
  lowbnd=Ybnd(1); highbnd=Ybnd(2);
end

%% plot data
%cm=[cm(1:end-1) ');'];
%eval(cm);
axis([xmin xmax lowbnd highbnd])
grid on

xlabel(xlabelstr);
ylabel(ylabelstr);
title(titlestr);

end


function b=phase(a)
  b=unwrap(angle(a));
end

function b=phase1(a)
  b=unwrappi(angle(a));
end

function b=unwrappi(a)
  c=round(diff(a)/pi)*pi;
  b=a-cumsum([0; c]);
end

function []=plot1 (x,y,logflag)
  if logflag,
    semilogx (x,y);
  else
    plot (x,y);
  end
end

function s=getMagAxisType(t)
	if t==1, s='dB';
	else s='linear';
	end
end

function s=getFreqAxisType(t)
	if t==1, s='logarithmic';
	else s='linear';
	end
end

function s=getTimeAxisType(t,tdflt)
	if t==tdflt, s='samples';
	else s='secs';
	end
end