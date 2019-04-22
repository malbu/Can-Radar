function varargout = daqradar(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @daqradar_OpeningFcn, ...
                   'gui_OutputFcn',  @daqradar_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code


% --- Executes just before daqradar is made visible.
function daqradar_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to daqradar (see VARARGIN)

% Choose default command line output for daqradar
handles.output = hObject;


%  I N P U T S  

%Selecting a source.
handles.adaptor     = 'winsound';
handles.id          = 0;
handles.chanup      = 1;
handles.chandown    = 2;
handles.samplesPerTrigger = 8192;
handles.FFTlength=8192;
handles.samplesPerChirp=56;
handles.sampleRate = 48000;
handles.numTraces = 30;   % number of traces to show in the intensity plot.
handles.HanningToggle = 0;       %Hanning window off by default
handles.RangingToggle =0; %Not ranging by default
handles.ChirpToggle=0;   %Faster trigger for edge extraction
handles.ParseNumber=1; %One parse display by default
handles.AverageParses=0; %Averaging off by default





%  F I G U R E 

%set(handles.figure1,'Color'.get
%Time Domain Axis
axes(handles.axes1);
handles.hLine1 = plot(zeros(1,handles.samplesPerTrigger)'); 
set(handles.hLine1,'Color', [.1 .1 0.5]);
set(handles.axes1,'Color',[235/255 255/255 235/255])
set(handles.axes1,'XGrid','on','YGrid','on')
t=title('Data','Color',[.05 .05 .25],'FontWeight','Bold','FontSize',9);
xlabel('Time (s)','FontSize',8);
ylabel('Voltage (V)','FontSize',8);

axes(handles.axes5);
handles.hLine5 = plot(zeros(1,handles.samplesPerTrigger)'); 
set(handles.hLine5,'Color', [.1 .1 0.5]);
set(handles.axes5,'Color',[235/255 255/255 235/255])
set(handles.axes5,'XGrid','on','YGrid','on')
t=title('Sync','Color',[.05 .05 .25],'FontWeight','Bold','FontSize',9);
xlabel('Time (s)','FontSize',8);
ylabel('Voltage (V)','FontSize',8);

%Frequency Domain Axis
axes(handles.axes2);
if(handles.RangingToggle==1)
    set(handles.axes2,'Visible','off');
    set(handles.hLine2,'Visible','off');
end
handles.hLine2 = plot(zeros(1,handles.samplesPerTrigger/2)'); 
set(handles.hLine2,'Color', [.1 0.5 .1]);
set(handles.axes2,'Color',[235/255 255/255 235/255])
set(handles.axes2,'XGrid','on','YGrid','on')
t=title('Frequency Domain','Color',[.05 0.25 .05],'FontWeight','Bold','FontSize',9);
xlabel('Frequency (Hz)','FontSize',8);
ylabel('Magnitude (dB)','FontSize',8);

axes(handles.axes3);
set(handles.axes3,'View',[90 -90]);
set(handles.axes3,'Color',[255/255 255/255 255/255]);
grid(handles.axes3,'on');
h = get(handles.axes3,'title');
set(h,'string','Intensity Plot','FontWeight','Bold','Color',[.25 .05 .05],'FontSize',9);
h = get(handles.axes3,'ylabel');
set(h,'string','Frequency (Hz)','FontSize',8);
h = get(handles.axes3,'zlabel');
set(h,'string','Magnitude (dB)','FontSize',8);

axes(handles.axes6);
handles.hLine6 = plot(zeros(1,handles.samplesPerTrigger)'); 
%handles.hLine6b = plot(zeros(1,handles.samplesPerTrigger)'); 
if(handles.RangingToggle==0)
    set(handles.axes6,'Visible','off');
    set(handles.hLine6,'Visible','off');
    %set(handles.hLine6b,'Visible','off');
end
set(handles.hLine6,'Color', [.1 0.5 .1]);
%set(handles.hLine6b,'Color', [0 0 0]);
set(handles.axes6,'Color',[235/255 255/255 235/255])
set(handles.axes6,'XGrid','on','YGrid','on')
t=title('First Parse','Color',[.05 .05 .25],'FontWeight','Bold','FontSize',9);
xlabel('Time (s)','FontSize',8);
ylabel('Voltage (V)','FontSize',8);

axes(handles.axes7);
handles.hLine7 = plot(zeros(1,handles.samplesPerTrigger)');  
if(handles.RangingToggle==0)
    set(handles.axes7,'Visible','off');
    set(handles.hLine7,'Visible','off');
end
set(handles.hLine7,'Color', [1 0 0]);
set(handles.axes7,'Color',[235/255 255/255 235/255])
set(handles.axes7,'XGrid','on','YGrid','on')
t=title('Subtracted Frames','Color',[.05 .05 .25],'FontWeight','Bold','FontSize',9);
xlabel('Time (s)','FontSize',8);
ylabel('Voltage (V)','FontSize',8);

axes(handles.axes8);
handles.hLine8 = plot(zeros(1,handles.samplesPerTrigger/2)');
if(handles.RangingToggle==0)
    set(handles.axes8,'Visible','off');
    set(handles.hLine8,'Visible','off');
end
set(handles.hLine8,'Color', [.1 0.5 .1]);
set(handles.axes8,'Color',[235/255 255/255 235/255])
set(handles.axes8,'XGrid','on','YGrid','on')
t=title('FFT','Color',[.05 0.25 .05],'FontWeight','Bold','FontSize',9);
xlabel('Frequency (Hz)','FontSize',8);
ylabel('Magnitude (dB)','FontSize',8);

set(hObject,'RendererMode','Manual')  %  If you don't do this, the surface plot
set(hObject,'Renderer','OpenGL')      %    will draw VERY slowly.


set(handles.poSampleRate,'String',[{'48000'},{'44100'},{'22000'},{'8000'}]);
set(handles.poBufferSize,'String',[{'8192'},{'4096'},{'1024'},{'512'}]);
set(handles.poFFTSize,'String',[{'8192'},{'4096'},{'1024'},{'512'}]);
set(handles.rbOnePulse,'Enable','off');
set(handles.rbTwoPulse,'Enable','off');
set(handles.cbAverageParses,'Enable','off');


ai=localSetupAI(handles);
handles.ai = ai;

% Update handles structure
guidata(hObject, handles);

localStartAI(ai);




function localStartAI(ai)

% S T A R T  A I 

start(ai);
trigger(ai);

%end localStartAI



function localStopAI(ai)

%  S T O P  A I 

stop(ai);
delete(ai);

% end localStopAI



function ai=localSetupAI(handles)

%  S E T U P   T H E   A N A L O G   I N P U T 


% Define object and add channels
ai = analoginput(handles.adaptor, handles.id);
addchannel(ai, handles.chanup);
addchannel(ai, handles.chandown);

% Configure the callback to update the display.
set(ai, 'TimerFcn', @localfftShowData);

% Configure the analog input object.
set(ai, 'SampleRate', handles.sampleRate);

% Configure the analog input object to trigger manually twice.
%  We do this because we are using peekdata to acquire the data in
%   a timer callback function.
%  The first trigger will fill the buffer with handles.samplesPerTrigger
%   number of samples.  We'll know we have enough samples to start 
%   processing data when the analog input object's SamplesAvailable property
%   is equal to handles.samplesPerTrigger.
%  The analog input object will then wait for 
%   another manual trigger, and while it is waiting the object will still be 
%   in its running state, which means the timer event will run. To keep the
%   object in the running state, we need only never manually trigger this
%   second trigger.  
%  Had we set the TriggerRepeat to 0, the analog input object would stop 
%   after the first trigger and the timer functions would stop running.
%
set(ai, 'SamplesPerTrigger', handles.samplesPerTrigger);
disp(handles.samplesPerTrigger);
set(ai, 'TriggerRepeat', 1);
set(ai, 'TriggerType', 'manual');

% Initialize callback parameters.  The TimerAction is initialized 
% after figure has been created.
set(ai, 'TimerPeriod', 0.01);  
set(ai, 'BufferingConfig',[handles.samplesPerTrigger*2,20]);

% Initialize time and frequency plots with lines of y=0
d=zeros(1,handles.samplesPerTrigger);
time = 1:handles.samplesPerTrigger;
f=1:handles.samplesPerTrigger/2;
mag=zeros(1,handles.samplesPerTrigger/2);

% Store state information in the analog input objects UserData area.
data.storedFFTsIndex = 1;
data.plotSurf        = 0;
data.ai              = ai;
data.getdata         = [d time];
data.daqfft          = [f mag];
data.handle          = [];
data.figureHandles   = handles;
%data.view            = [103 10];
%data.rotateStep      = 4;
data.counter         = 0;

% Set the object's UserData to data.
set(data.ai, 'UserData', data);
%end localSetupAI(handles)



% S E T U P   T H E   A N A L O G   I N P U T 


% --- Outputs from this function are returned to the command line.
function varargout = daqradar_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pbExit or when you press
%     the figure close 'X' button (I set this function to
%     the figures CloseRequestFcn in GUIDE)
function pbExit_Callback(hObject, eventdata, handles)
% hObject    handle to pbExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
localStopAI(handles.ai);
closereq;


% --- Executes during object creation, after setting all properties.
function poSampleRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poSampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



%  C H A N G E   T H E   S A M P L E   R A T E  

% --- Executes on selection change in poSampleRate.
function poSampleRate_Callback(hObject, eventdata, handles)

% hObject    handle to poSampleRate
% eventdata  
% handles    structure with handles and user data



% First, stop and delete the current analog input object
localStopAI(handles.ai);

% Extract the new samplerate.
v=get(handles.poSampleRate,'Value');
s=get(handles.poSampleRate,'String');
handles.sampleRate = str2num(s{v});

% Create a new analog input with the new sample rate.
handles.ai = localSetupAI(handles);

% Update handles structure
guidata(hObject, handles);

% Restart the analog input
localStartAI(handles.ai);

% end poSampleRate_Callback



% C H A N G E   T H E   SAMPLES/TRIGGER   

% --- Executes on selection change in poBufferSize.
function poBufferSize_Callback(hObject, eventdata, handles)

% First, stop and delete the current analog input object
localStopAI(handles.ai);

% Extract the new size.
v=get(handles.poBufferSize,'Value');
s=get(handles.poBufferSize,'String');
handles.samplesPerTrigger = str2num(s{v});

% Create a new analog input with the new size.
handles.ai = localSetupAI(handles);

% Update handles structure
guidata(hObject, handles);

% Restart the analog input
localStartAI(handles.ai);
% end poSampleRate_Callback


% C H A N G E   T H E   FFT size 

% --- Executes on selection change in poFFTSize.
function poFFTSize_Callback(hObject, eventdata, handles)
% hObject    handle to poFFTSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% First, stop and delete the current analog input object
localStopAI(handles.ai);

% Extract the new size.
v=get(handles.poBufferSize,'Value');
s=get(handles.poBufferSize,'String');
handles.FFTlength = str2num(s{v});

% Create a new analog input with the new size.
handles.ai = localSetupAI(handles);

% Update handles structure
guidata(hObject, handles);

% Restart the analog input
localStartAI(handles.ai);
% end poFFTSize_Callback

% --- Executes during object creation, after setting all properties.
function poFFTSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poFFTSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% ***********************************************************************  
% Calculate the fft of the data.  (Copied from demoai_fft.m)
function [f, mag] = localDaqfft(data,Fs,blockSize)

% Calculate the fft of the data.

%xFFT = fft(data(:,2)); %channel 2
%xFFT = fft(data(:,1)); %channel 1
xFFT = fft(data);
xfft = abs(xFFT);

% Avoid taking the log of 0.
index = find(xfft == 0);
xfft(index) = 1e-17;

mag = 20*log10(xfft);
mag = mag(1:floor(blockSize/2));

f = (0:length(mag)-1)*Fs/blockSize; 
f = f(:);


% ***********************************************************************  
% Update the plot. This routine is a Timer callback, it is called
%  automatically at a preset time interval. See line 144 for where
%  this routine is assigned as a callback
function localfftShowData(obj,event)

if (get(obj,'SamplesAvailable') >= obj.SamplesPerTrigger)
	
	% Get the handles.
	data = obj.UserData;
   
    global x; 
    global f;
    global mag; 
    global parse;
   global parsediff;
   global ind;
	
	handles = data.figureHandles;
	
	% Execute a peekdata.
	x = peekdata(obj, obj.SamplesPerTrigger);
    
   
    
    % Applying Hanning window if checked
	
     if(handles.HanningToggle==1 && handles.ChirpToggle==0)
        win=hanning(get(obj, 'SamplesPerTrigger'));
         x(:,2)=win(:).*x(:,2); %chandown
        %disp('Hello')
        
     end 
     
    Fs = obj.SampleRate;
	blockSize = obj.SamplesPerTrigger;
    % disp(sprintf('FS: %d',Fs));
    
    %disp(sprintf('Block Size: %d',blockSize));
    %disp(sprintf('Chirp Size: %d',handles.samplesPerChirp)); 
    
    % Edge Detection for Ranging
    %disp(sprintf('Ranging Toggle: %d',handles.RangingToggle));
    if(handles.ChirpToggle==1)
        global signal;
        sync=x(:,1);
        signal=x(:,2);
        
        %h= fdesign.lowpass('N,Fst,Ast', 10, 500,80, Fs);
        %Hd=design(h,'cheby2');
        %signal=filter(Hd,signal);
        
        
        
        tp=.5; %clock is good enough for ~.8. Anything highers risks false edges
        tn=-.5;
        atp=sync(:,1)>tp;
        atn=sync(:,1)<tn;
        atp=atp';
        atn=atn';
        global risingedge;
        global fallingedge;
        risingedge=diff(atp)>0;
        risingedge=risingedge';
        fallingedge=diff(atn)>0;
        fallingedge=fallingedge';
       
        global indrising;
        global indfalling;
        if(risingedge==[0]) %prevent initialization break because of empty array
            return;
        end
        indrising=find(risingedge); 
        indfalling=find(fallingedge);
        % index position of start
       % end
        global fsize;
        fsize=zeros(1,handles.samplesPerChirp+1);
        if(indrising(1)<indfalling(1))          %first edge is rising
            for k=1:1:length(indfalling)-1
                fsize(:,k)=indfalling(k)-indrising(k);
            end
        elseif(indrising(1)>indfalling(1))      %first edge is falling
            for k=2:1:length(indfalling)
                fsize(:,k-1)=indfalling(k)-indrising(k-1);
            end
        end
        fsize=fsize(find(fsize));
        
        global chirpsize;
        chirpsize=mode(fsize);
        
        %disp(sprintf('Index: %d',length(ind))); 

        parse=zeros(length(fsize),mode(fsize));
        
        global j;
        for j=1:1:length(fsize)-1;       %parsing data from positive slope 
            if(fsize(:,j)==fsize(:,j+1) && fsize(:,j)==chirpsize) 
                parse(j,:)=signal(indrising(j,:):indrising(j,:)+fsize(:,j)-1);
            end
        end
         
        global bprune;
        bprune=parse;
        
        parse(all(parse==0,2),:)=[]; %deleting zero rows
        
        if(size(parse,1)<2) %need at least 2 same size consecutive frames
            return;
        end
        parse=parse';
        
        
        [nrows ncols]=size(parse);
        splitpoint=floor(ncols/2);
        global parse1;
        global parse2;
        
        parse1=parse(:,1:splitpoint);
        parse2=parse(:,splitpoint+1:2*splitpoint);
        parse1=mean(parse1,2);
        parse2=mean(parse2,2);
        
        
        %subtracting consecutive frames
        global start;
        global stop;
        global pad;
        global startz;
        global stopz;
        global padded2;
        global padded1;
        global time_padded1;
        global time_padded2;
        global padded_diff;
        if(handles.AverageParses==0) %no average 
            time_padded1=parse(:,1)';
            time_padded2=parse(:,2)';
        elseif(handles.AverageParses==1)
            time_padded1=parse1';
            time_padded2=parse2';
        end
        
        %h= fdesign.lowpass('N,F3dB', 10, 5000, Fs);
        %Hd=design(h,'butter');
        %time_padded1=filter(Hd,time_padded1);
        %time_padded2=filter(Hd,time_padded2);
        if(handles.HanningToggle==1)
            win2=hanning(chirpsize)';
            time_padded1=win2.*time_padded1;
            time_padded2=win2.*time_padded2;
        end 
        
       
        
        %Zero-Padding
        pad=handles.FFTlength-chirpsize;
        start=floor(pad/2);
        stop=handles.FFTlength-(pad-start);
        %stop=pad-start;
        parsediff=parse(:,2)-parse(:,1);
        startz=zeros(1,start);
        stopz=zeros(1,handles.FFTlength-stop);
        padded1=[startz time_padded1(1:end) stopz];
        padded2=[startz time_padded2(1:end) stopz];
        padded_diff=padded2-padded1;
   
        %disp(sprintf('Parse Ranging: %d',length(parse)));
       %disp(sprintf('F Ranging: %d',length(f)));
       % disp(sprintf('Mag Ranging: %d',length(mag)));
       [f,mag] = localDaqfft(x(:,2),Fs,blockSize);
       global f1;
       global mag1;
       global f2;
       global mag2;
       global f3;
       global mag3;
       global f4;
       global mag4;
       global f5;
       global mag5;
       global f6;
       global mag6;
       [f1,mag1]=localDaqfft(parse(:,1),Fs,chirpsize); % FFT of first parsed frame
       [f2,mag2]=localDaqfft(parse(:,2),Fs,chirpsize); % FFT of second parsed frame
       [f3,mag3]=localDaqfft(parsediff,Fs,chirpsize); %FFT of difference 
       [f4,mag4]=localDaqfft(padded1,Fs,handles.FFTlength); %padded FFT 1
       [f5,mag5]=localDaqfft(padded2,Fs,handles.FFTlength); %padded FFT 2
       [f6,mag6]=localDaqfft(padded_diff,Fs,handles.FFTlength); %FFT of padded diff
       if(handles.ParseNumber==1)
           f6=f4;
           mag6=mag4;
       end
    else %not ranging
        [f,mag] = localDaqfft(x(:,2),Fs,blockSize);
        %disp(size(f))
        disp(sprintf('F Regular: %d',length(f)));
        disp(sprintf('Mag Regular: %d',length(mag)));
    end
    %Dynamically modify Analog Chan 1 (chanup=sync)
	maxX=max(x(:,1));
	minX=min(x(:,1));
    yax3=get(handles.axes5,'YLim');
    yax3(1)=minX-.0001;
    yax3(2)=maxX+.0001;
    set(handles.axes5,'YLim',yax3)
	set(handles.axes5,'XLim',[0 (obj.SamplesPerTrigger-1)/obj.SampleRate])
    
	% Dynamically modify Analog Chan 2 (chandown=data)
	maxX2=max(x(:,2));
	minX2=min(x(:,2));
	yax1=get(handles.axes1,'YLim');
	yax1(1)=minX2 - .0001; % need to subtract a value to make sure yax(1) never equals yax(2)
	yax1(2)=maxX2 + .0001; 
	set(handles.axes1,'YLim',yax1)
	set(handles.axes1,'XLim',[0 (obj.SamplesPerTrigger-1)/obj.SampleRate])
	
	% Dynamically modify Frequency axis as we go.
	maxF=max(f);
	minF=min(f);
	xax=get(handles.axes2,'XLim');
    xax(1)=minF;
    xax(2)=maxF;
	set(handles.axes2,'XLim',xax)
	
	% Dynamically modify Magnitude axis as we go.
	maxM=max(mag);
	minM=min(mag);
	yax2=get(handles.axes2,'YLim');
	yax2(1)=minM - .0001;
	yax2(2)=maxM + .0001;
	set(handles.axes2,'YLim',yax2)
    
    
    
    if(handles.ChirpToggle==1)
    %Dynamically modify Parse Frame 
	maxX3=max(parse(:,1));
	minX3=min(parse(:,1));
    yax4=get(handles.axes6,'YLim');
    yax4(1)=minX3-.0001;
    yax4(2)=maxX3+.0001;
    set(handles.axes6,'YLim',yax4)
	set(handles.axes6,'XLim',[0 (chirpsize-1)/obj.SampleRate])
    %disp(obj.SamplesPerTrigger);
    
    %Dynamically modify Difference Frame 
	maxX4=max(parsediff);
	minX4=min(parsediff);
    yax5=get(handles.axes7,'YLim');
    yax5(1)=minX4-.0001;
    yax5(2)=maxX4+.0001;
    set(handles.axes7,'YLim',yax5)
	set(handles.axes7,'XLim',[0 (chirpsize-1)/obj.SampleRate])
    
    % Dynamically modify Frequency axis DifferenceFFT
	maxF3=max(f6);
	minF3=min(f6);
	xax2=get(handles.axes8,'XLim');
    xax2(1)=minF3;
    xax2(2)=maxF3;
	set(handles.axes8,'XLim',xax2)
    
    % Dynamically modify Mag axis Difference FFT
    maxM3=max(mag6);
	minM3=min(mag6);
	yax3=get(handles.axes8,'YLim');
	yax3(1)=minM3 - .0001;
	yax3(2)=maxM3 + .0001;
	set(handles.axes8,'YLim',yax3)
    end
    
    
	
	% Update the line plots.
	set(handles.hLine1, 'XData', [0:(obj.SamplesPerTrigger-1)]/obj.SampleRate, 'YData', x(:,2)); 
    set(handles.hLine5, 'XData', [0:(obj.SamplesPerTrigger-1)]/obj.SampleRate, 'YData', x(:,1));
    set(handles.hLine2, 'XData', f(:,1), 'YData', mag(:,1));
    if(handles.ChirpToggle==1)
    set(handles.hLine6, 'XData', [0:(chirpsize-1)]/obj.SampleRate, 'YData', parse(:,1)); 
    set(handles.hLine7, 'XData', [0:(chirpsize-1)]/obj.SampleRate, 'YData', parsediff); 
    set(handles.hLine8, 'XData', f6, 'YData', mag6);
    %set(handles.hLine6b, 'XData', [0:(handles.samplesPerChirp-1)]/obj.SampleRate, 'YData', parse(:,2));
    end
    
    

    % Find the frequency at which the max signal strength is at.
    %Store FFT into array for intensity plot
    if(handles.ChirpToggle==1)
        [yrange,rangeindex] = max(mag6);
        set(handles.tFreq,'String',sprintf('%4.3d Hz',f6(rangeindex)));
        data.storedFFTs(data.storedFFTsIndex,:) = mag6';
    elseif (handles.RangingToggle==0)    
        [ymax,maxindex] = max(mag);
        set(handles.tFreq,'String',sprintf('%4.3d Hz',f(maxindex)));
        data.storedFFTs(data.storedFFTsIndex,:) = mag';
    end
   
	
    % This circular shift is used so that when we display the 3D plot, the 
    %  newest FFT will appear in 'front' and the oldest in 'back'. 
    % To understand this, note how the plotting routines are using this fftOrder 
    %  array to reorder the FFTs stored in data.storedFFTs and also note
    %  how data.storedFFTsIndex is used to store FFTs in data.storedFFTs.
    %
	fftOrder = 1:handles.numTraces;
	fftOrder = circshift(fftOrder,[ 1 -data.storedFFTsIndex ]);  
	
	data.storedFFTsIndex = data.storedFFTsIndex + 1;
	if (data.storedFFTsIndex > handles.numTraces)
        data.storedFFTsIndex = 1;
        data.plotSurf        = 1; % Indicates a full history is stored.
	end
	
	% Update the surface plot if we have a full history.
	if(data.plotSurf)
        cla(handles.axes3);
        data.view  = [180 -90];
        if(handles.RangingToggle==0) 
            data=localClassic(handles,data,f,fftOrder);
        elseif(handles.ChirpToggle==1)
            data=localClassic(handles,data,f6,fftOrder);
        end
    end
	
	set(data.ai, 'UserData', data);
	
	drawnow;
end



%  I N T E N S I T Y   G R A P H

function data=localClassic(handles,data,f,fftOrder)
    [X,Y] = meshgrid(1:handles.numTraces,f(1:end));
    surf(X,Y,data.storedFFTs(fftOrder,:)','parent',handles.axes3,'EdgeColor','none');                
  	set(handles.axes3,'XLim',[1 handles.numTraces],'YLim',[0 f(end)])
    shading(handles.axes3,'interp');
    set(handles.axes3,'View',data.view)



    
% --- Executes during object creation, after setting all properties.
function poPlotType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poPlotType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function poBufferSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poBufferSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in cbHanningToggle.
function cbHanningToggle_Callback(hObject, eventdata, handles)
% hObject    handle to cbHanningToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% hObject    handle to pbHanningToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pbHanningToggle

% First, stop and delete the current analog input object
localStopAI(handles.ai);

%
toggle=get(handles.cbHanningToggle,'Value');
%disp(toggle)

handles.HanningToggle = toggle;

% Create a new analog input with the new state.
handles.ai = localSetupAI(handles);

% Update handles structure
guidata(hObject, handles);

% Restart the analog input
localStartAI(handles.ai);

% end cbHanningToggle_Callback




% --- Executes on button press in cb_ranging.
function cb_ranging_Callback(hObject, eventdata, handles)
% hObject    handle to cb_ranging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_ranging
localStopAI(handles.ai);
ranging=get(handles.cb_ranging,'Value');
handles.RangingToggle=ranging;

if(ranging==1)
    set(handles.text_chirpf,'Enable','on'); %enable text Box
    set(handles.rbOnePulse,'Enable','on');
    set(handles.rbTwoPulse,'Enable','on');
    set(handles.cbAverageParses,'Enable','on');
    set(handles.axes6,'Visible','on');
    set(handles.axes7,'Visible','on');
    set(handles.hLine6,'Visible','on');
    set(handles.hLine7,'Visible','on');
    set(handles.axes8,'Visible','on');
    set(handles.hLine8,'Visible','on');
    %set(handles.hLine6b,'Visible','on');
    set(handles.axes2,'Visible','off');
    set(handles.hLine2,'Visible','off');
else
    set(handles.text_chirpf,'Enable','off');
    set(handles.rbOnePulse,'Enable','off');
    set(handles.rbTwoPulse,'Enable','off');
    set(handles.cbAverageParses,'Enable','off');
    set(handles.axes6,'Visible','off');
    set(handles.axes7,'Visible','off');
    set(handles.hLine6,'Visible','off');
    set(handles.hLine7,'Visible','off');
    %set(handles.hLine6b,'Visible','off');
    set(handles.axes2,'Visible','on');
    set(handles.hLine2,'Visible','on');
    set(handles.axes8,'Visible','off');
    set(handles.hLine8,'Visible','off');
    handles.ChirpToggle=0;
    
end   
handles.ai = localSetupAI(handles);

% Update handles structure
guidata(hObject, handles);

% Restart the analog input
localStartAI(handles.ai);
% end cb_ranging_Callback


function text_chirpf_Callback(hObject, eventdata, handles)
% hObject    handle to text_chirpf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_chirpf as text
%        str2double(get(hObject,'String')) returns contents of text_chirpf as a double

sample_rate=handles.sampleRate;
%disp(sample_rate);

localStopAI(handles.ai);

chirp_freq=str2double(get(handles.text_chirpf,'String'));
chirpSamples=round(1/(2*chirp_freq)*sample_rate);  
%num_samples=chirpSamples*6;   %new window size want at least 2 positive edge for difference
%disp(num_samples)

%handles.samplesPerTrigger =num_samples ;
handles.samplesPerChirp=chirpSamples;
handles.ChirpToggle=1; %faster trigger than RangingToggle
disp(sprintf('ChirpToggle inside text_chripf: %d',handles.ChirpToggle));

% Create a new analog input with the new size.
handles.ai = localSetupAI(handles);

% Update handles structure
guidata(hObject, handles);

% Restart the analog input
localStartAI(handles.ai);
% end text_chirpf_Callback



% --- Executes during object creation, after setting all properties.
function text_chirpf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_chirpf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in cbAverageParses.
function cbAverageParses_Callback(hObject, eventdata, handles)
% hObject    handle to cbAverageParses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
localStopAI(handles.ai);

average=get(handles.cbAverageParses,'Value');
handles.AverageParses=average;
% Create a new analog input with the new state.
handles.ai = localSetupAI(handles);

% Update handles structure
guidata(hObject, handles);

% Restart the analog input
localStartAI(handles.ai);


% --- Executes when selected object is changed in uipanel6.
function uipanel6_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel6 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

localStopAI(handles.ai);

%
parsenumber=get(eventdata.NewValue,'tag');

switch parsenumber
    case 'rbOnePulse'
    handles.ParseNumber = 1; 
    case 'rbTwoPulse'
    handles.ParseNumber = 0;
end

% Create a new analog input with the new state.
handles.ai = localSetupAI(handles);

% Update handles structure
guidata(hObject, handles);

% Restart the analog input
localStartAI(handles.ai);
