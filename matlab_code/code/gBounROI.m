function varargout = gBounROI(varargin)
% GBOUNROI MATLAB code for gBounROI.fig
%      GBOUNROI, by itself, creates a new GBOUNROI or raises the existing
%      singleton*.
%
%      H = GBOUNROI returns the handle to a new GBOUNROI or the handle to
%      the existing singleton*.
%
%      GBOUNROI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GBOUNROI.M with the given input arguments.
%
%      GBOUNROI('Property','Value',...) creates a new GBOUNROI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gBounROI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gBounROI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%   COPYRIGHT TBB 2012
%   For use only within ox.ac.uk domain
%   version 1.0
% Edit the above text to modify the response to help gBounROI

% Last Modified by GUIDE v2.5 27-Feb-2012 17:39:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gBounROI_OpeningFcn, ...
                   'gui_OutputFcn',  @gBounROI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1}); % converts str to function handle
end

if nargout %this checks if outout exists
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gBounROI is made visible.
function gBounROI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gBounROI (see VARARGIN)

% Choose default command line output for gBounROI
handles.output = hObject;

%load the image
I1=varargin{1};

%set up inputs
filter_set=varargin{2};
roisize=varargin{3};

%set the max slider values
gSetFilterSlideMax(roisize)

%set the filter values in the string boxes
gSetFilter(filter_set)


%create h window
hfilter=HWindow(roisize);



%plot the pattern
axes(handles.Pattern);
imagesc(I1)
axis off; axis image
colormap('gray')

%determine the starting boundary values
I1_size=size(I1);
boundary_par=[1,1,I1_size(2),I1_size(1)];
P1=PlotBoundary(boundary_par,handles.Pattern);

BoundaryData.P1=P1;
BoundaryData.BoundaryPar=boundary_par;

set(handles.PickBoundary,'UserData',BoundaryData);

%start with a ROI in the centre
roiloc=round([I1_size(2)/2,I1_size(1)/2]);

iData=struct;
iData.roisize=roisize;
iData.isize=I1_size;
iData.I1=I1;
iData.hfilter=hfilter;

%plot the ROI box
[roi_Xmax,roi_Ymax,roi_Xmin,roi_Ymin]=ROIpos(roiloc,iData.roisize);

%plot the ROI
axes(handles.Pattern);
hold on;
roi_bounds=[ roi_Xmin,roi_Xmax,roi_Xmax,roi_Xmin,roi_Xmin;
             roi_Ymax,roi_Ymax,roi_Ymin,roi_Ymin,roi_Ymax];
roibox_h=plot(roi_bounds(1,:),roi_bounds(2,:),'w','Linewidth',2);
hold off
set(handles.PickTestROI,'UserData',roibox_h);

set(handles.PickTestROI,'UserData',roibox_h);
set(handles.Pattern,'UserData',iData);

%set up the ROI data and plot in 1st image box
%   calc the raw FFT and plot in 2nd image box
[roi_data,fft_roi]=ROIExtract(I1,hfilter,roi_Xmin,roi_Ymin,roi_Xmax,roi_Ymax,handles.WROI,handles.FFT);

%set up filter
FFTfilter=FFT_Filter(filter_set,roisize);
%filter FFT and plot filtered FFT and BT
[fft_roi_filt,roi_bt]=FFTFilterCreate(fft_roi,FFTfilter,handles.FFT_Filt,handles.BTROI);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gBounROI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gBounROI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% fpassset=gGetFilters;
% BoundaryData=get(handles.PickBoundary,'UserData');
% boundary_par=BoundaryData.BoundaryPar;
fpassset=gGetFilters;
BoundaryData=get(handles.PickBoundary,'UserData');
varargout{1} = fpassset(:);
varargout{2} = BoundaryData.BoundaryPar(:);

delete(handles.figure1);


% --- Executes on slider movement.
function Slide_HPW_Callback(hObject, eventdata, handles)
% hObject    handle to Slide_HPW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%update the read out
sval=get(hObject,'Value');
HH = findobj(gcf,'Tag','Text_HPW');
set(HH,'String',int2str(sval));

%update the filter
fpassset=gGetFilters;

iData=get(handles.Pattern,'UserData');
FFTfilter=FFT_Filter(fpassset,iData.roisize);
fft_roi=get(handles.FFT,'UserData');
[~,~]=FFTFilterCreate(fft_roi,FFTfilter,handles.FFT_Filt,handles.BTROI);


% --- Executes during object creation, after setting all properties.
function Slide_HPW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slide_HPW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Slide_LP_Callback(hObject, eventdata, handles)
% hObject    handle to Slide_LP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%update the read out
sval=get(hObject,'Value');
HH = findobj(gcf,'Tag','Text_LP');
set(HH,'String',int2str(sval));

%update the filter
fpassset=gGetFilters;
iData=get(handles.Pattern,'UserData');
FFTfilter=FFT_Filter(fpassset,iData.roisize);
fft_roi=get(handles.FFT,'UserData');
[~,~]=FFTFilterCreate(fft_roi,FFTfilter,handles.FFT_Filt,handles.BTROI);


% --- Executes during object creation, after setting all properties.
function Slide_LP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slide_LP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Slide_LPW_Callback(hObject, eventdata, handles)
% hObject    handle to Slide_LPW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%update the read out
sval=get(hObject,'Value');
HH = findobj(gcf,'Tag','Text_LPW');
set(HH,'String',int2str(sval));

%update the filter
fpassset=gGetFilters;

iData=get(handles.Pattern,'UserData');
FFTfilter=FFT_Filter(fpassset,iData.roisize);
fft_roi=get(handles.FFT,'UserData');
[~,~]=FFTFilterCreate(fft_roi,FFTfilter,handles.FFT_Filt,handles.BTROI);


% --- Executes during object creation, after setting all properties.
function Slide_LPW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slide_LPW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Slide_HP_Callback(hObject, eventdata, handles)
% hObject    handle to Slide_HP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%update the read out
sval=get(hObject,'Value');
HH = findobj(gcf,'Tag','Text_HP');
set(HH,'String',int2str(sval));

%update the filter
fpassset=gGetFilters;

iData=get(handles.Pattern,'UserData');
FFTfilter=FFT_Filter(fpassset,iData.roisize);
fft_roi=get(handles.FFT,'UserData');
[~,~]=FFTFilterCreate(fft_roi,FFTfilter,handles.FFT_Filt,handles.BTROI);

 
% --- Executes during object creation, after setting all properties.
function Slide_HP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slide_HP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in PickBoundary.
function PickBoundary_Callback(hObject, eventdata, handles)
% hObject    handle to PickBoundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
BoundaryData=get(handles.PickBoundary,'UserData');


[boundary_par]=ClickBoundary(handles.Pattern);
                if exist('S1','var') == 1
            if ishandle(S1) == 1
                delete(S1);
            end
                end
        
if isempty(boundary_par) == 0
    BoundaryData.BoundaryPar=boundary_par;
end

BoundaryData.P1=PlotBoundary(BoundaryData.BoundaryPar,handles.Pattern,BoundaryData.P1);

set(handles.PickBoundary,'UserData',BoundaryData);


% --- Executes on button press in PickTestROI.
function PickTestROI_Callback(hObject, eventdata, handles)
% hObject    handle to PickTestROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get pattern data
iData=get(handles.Pattern,'UserData');
BoundaryData=get(handles.PickBoundary,'UserData');

%mouse input for new ROI loc
roiloc=ginput(1);
roiloc=round([roiloc(1,1),roiloc(1,2)]);

%determine the roi bounds

%determine the ROI box
[roi_Xmax,roi_Ymax,roi_Xmin,roi_Ymin]=ROIpos(roiloc,iData.roisize);

%check if these points (1) lie on the screen, (2) lie within the boundary
points_check=[roi_Xmax,roi_Ymax;roi_Xmax,roi_Ymin;roi_Xmin,roi_Ymax;roi_Xmin,roi_Ymin];

if all(points_check(:) > 0) && all(points_check(:) <  iData.isize(1))
        
        %delete the current ROI box
        roibox_h=get(handles.PickTestROI,'UserData');
        delete(roibox_h);
        
        %plot the ROI
        roi_bounds=[ roi_Xmin,roi_Xmax,roi_Xmax,roi_Xmin,roi_Xmin;
            roi_Ymax,roi_Ymax,roi_Ymin,roi_Ymin,roi_Ymax];
        axes(handles.Pattern);
        hold on;
        roibox_h=plot(roi_bounds(1,:),roi_bounds(2,:),'w','Linewidth',2);
        hold off
        set(handles.PickTestROI,'UserData',roibox_h);
        
        %set up the ROI data and plot in 1st image box
        %   calc the raw FFT and plot in 2nd image box
        [roi_data,fft_roi]=ROIExtract(iData.I1,iData.hfilter,roi_Xmin,roi_Ymin,roi_Xmax,roi_Ymax,handles.WROI,handles.FFT);
        
        %set up filter
        fpassset=gGetFilters;
        FFTfilter=FFT_Filter(fpassset,iData.roisize);
        %filter FFT and plot filtered FFT and BT
        [fft_roi_filt,roi_bt]=FFTFilterCreate(fft_roi,FFTfilter,handles.FFT_Filt,handles.BTROI);
end


% --- Executes on button press in SaveAndContinue.
function SaveAndContinue_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAndContinue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure1_CloseRequestFcn(hObject, eventdata, handles)


function [roi_Xmax,roi_Ymax,roi_Xmin,roi_Ymin]=ROIpos(roiloc,roisize)
%plot the box with the roi location on the EBSP
roi_Xmax=floor(roiloc(1,1)+roisize/2+0.5);
roi_Ymax=floor(roiloc(1,2)+roisize/2+0.5);
roi_Xmin=ceil(roiloc(1,1)-roisize/2+0.5);
roi_Ymin=ceil(roiloc(1,2)-roisize/2+0.5);


function [roi_data,fft_roi]=ROIExtract(I1,hfilter,roi_Xmin,roi_Ymin,roi_Xmax,roi_Ymax,handle_wroi,handle_fft)
%extract ROI data from image and plot
roi_data=I1(roi_Ymin:roi_Ymax,roi_Xmin:roi_Xmax);
roi_data=hfilter.*(roi_data-mean(roi_data(:)))/std(roi_data(:));

axes(handle_wroi);
imagesc(roi_data);
axis off; axis image; 

%calc 2D fft and plot
fft_roi=fftshift(fft2(roi_data));
axes(handle_fft);
imagesc(log10(abs(real(fft_roi))));
axis off; axis image; c=caxis; caxis([0 c(2)])

set(handle_fft,'UserData',fft_roi)


function [fft_roi_filt,roi_bt]=FFTFilterCreate(fft_roi,FFTfilter,handle_fftfilt,handle_bt)
%filter FFT
fft_roi_filt=fft_roi.*fftshift(FFTfilter);
axes(handle_fftfilt);
imagesc(log10(abs(real(fft_roi_filt))));
axis off; axis image; c=caxis; caxis([0 c(2)])

%back transform
roi_bt=ifft2(fftshift(fft_roi_filt));
axes(handle_bt);
imagesc(real(roi_bt));
axis off; axis image; 


function hfilter=HWindow(roisize)
%create window function
u=1:roisize;
[meshv,meshu]=meshgrid(u,u);
meshvf=meshv-roisize/2-0.5;
meshuf=meshu-roisize/2-0.5;

%create the Hfilter
hfilter=cos(pi.*(meshuf)/roisize).*cos(pi.*(meshvf)/roisize);


function [FFTfilter]=FFT_Filter(fpassset,roisize)
%create the Fourier domain filter
lcutoff=fpassset(3);
lwidth=fpassset(4)/2;
hcutoff=fpassset(1);
hwidth=fpassset(2)/2;

if lcutoff < hcutoff
    error('The low pass filter is smaller than the high pass filter');
end

%generate an x and y grid
u=1:roisize;
[meshv,meshu]=meshgrid(u,u);
meshvf=meshv-roisize/2-0.5;
meshuf=meshu-roisize/2-0.5;

%create the Hfilter
hfilter=cos(pi.*(meshuf)/roisize).*cos(pi.*(meshvf)/roisize);

%create the FFTfilter

distf=sqrt(meshvf.*meshvf+meshuf.*meshuf);

%lowpass
lFFTfilter=exp(-((distf-lcutoff)/(sqrt(2)*lwidth/2)).^2);
lFFTfilter(distf>(lcutoff+2*lwidth))=0;
lFFTfilter(distf<lcutoff)=1;
%highpass
hFFTfilter=exp(-((hcutoff-distf)/(sqrt(2)*hwidth/2)).^2);
hFFTfilter(distf<(hcutoff-2*hwidth))=0;
hFFTfilter(distf>hcutoff)=1;

%combine
FFTfilter=hFFTfilter.*lFFTfilter;
FFTfilter=fftshift(FFTfilter);


function gSetFilter(fpassset)
sl_objnames={'Slide_HP','Slide_HPW','Slide_LP','Slide_LPW'};
tx_objnames={'Text_HP','Text_HPW','Text_LP','Text_LPW'};

for n=1:4
    HH = findobj(gcf,'Tag',tx_objnames{n});
    set(HH,'String',int2str(fpassset(n)));
    HH = findobj(gcf,'Tag',sl_objnames{n});
    set(HH,'Value',fpassset(n));
end


function gSetFilterSlideMax(roisize)

flimit=floor(roisize/2);
objnames={'Slide_HP','Slide_HPW','Slide_LP','Slide_LPW'};
for n=1:4
    HH = findobj(gcf,'Tag',objnames{n});
    set(HH,'Max',flimit);
    set(HH,'SliderStep',[1/flimit 1/flimit]);
end


function fpassset=gGetFilters
fpassset=zeros(4,1);
objnames={'Slide_HP','Slide_HPW','Slide_LP','Slide_LPW'};
for n=1:4
    HH = findobj(gcf,'Tag',objnames{n});
fpassset(n)=get(HH,'Value');
end


function [boundary_par]=ClickBoundary(handle_pattern)
axes(handle_pattern)
hold on
    Boundary_corners = ginput(2);
    Boundary_corners = round(Boundary_corners);
    Boundary_height  = Boundary_corners(2,2)-Boundary_corners(1,2);
    Boundary_width   = Boundary_corners(2,1)-Boundary_corners(1,1);
    if Boundary_height < 0 && Boundary_width >0
        Rectangle_height = abs(Boundary_height);
        Rectangle_width  = Boundary_width;
        Rectangle_origin = [Boundary_corners(1,1),Boundary_corners(2,2)];     
    elseif Boundary_height > 0 && Boundary_width >0
        Rectangle_height = Boundary_height;
        Rectangle_width  = Boundary_width;
        Rectangle_origin = [Boundary_corners(1,1),Boundary_corners(1,2)];  
    elseif Boundary_height < 0 && Boundary_width <0
        Rectangle_height = -Boundary_height;
        Rectangle_width  = -Boundary_width;
        Rectangle_origin = [Boundary_corners(2,1),Boundary_corners(2,2)];
    elseif Boundary_height > 0 && Boundary_width <0
        Rectangle_height = -Boundary_height;
        Rectangle_width  = -Boundary_width;
        Rectangle_origin = [Boundary_corners(2,1),Boundary_corners(1,2)]; 
    end 
    
    boundary_par     = [Rectangle_origin(1,1),Rectangle_origin(1,2),Rectangle_width, Rectangle_height];
% the structure of boundary_par is: corner1_x,corner1_y,
% corner2_x,corner2_y, recentangle_width, rectangle_height
    
    hold off


function P1=PlotBoundary(boundary_par,image_handle,P1)
        if exist('P1','var') == 1
            if ishandle(P1) == 1
                delete(P1);
            end
        end
    axes(image_handle); hold on
    P1 = rectangle('position',[boundary_par(1) boundary_par(2),boundary_par(3), boundary_par(4)],'EdgeColor','r','LineWidth',2);    
    hold off


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

if isequal(get(handles.figure1,'WaitStatus'),'waiting')
    uiresume(handles.figure1);
else
    delete(handles.figure1);
end