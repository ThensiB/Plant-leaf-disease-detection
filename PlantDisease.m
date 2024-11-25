function varargout = PlantDisease(varargin)
% PLANTDISEASE M-file for PlantDisease.fig
%      PLANTDISEASE, by itself, creates a new PLANTDISEASE or raises the existing
%      singleton*.
%
%      H = PLANTDISEASE returns the handle to a new PLANTDISEASE or the handle to
%      the existing singleton*.
%
%      PLANTDISEASE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLANTDISEASE.M with the given input arguments.
%
%      PLANTDISEASE('Property','Value',...) creates a new PLANTDISEASE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PlantDisease_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PlantDisease_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PlantDisease

% Last Modified by GUIDE v2.5 04-Dec-2012 17:29:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PlantDisease_OpeningFcn, ...
                   'gui_OutputFcn',  @PlantDisease_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before PlantDisease is made visible.
function PlantDisease_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PlantDisease (see VARARGIN)

% Choose default command line output for PlantDisease
handles.output = hObject;

set(handles.color_transformation_structure_pushbutton,'enable','off');
set(handles.RGB_to_HSI_pushbutton,'enable','off');
set(handles.kmeans_pushbutton,'enable','off');
set(handles.masking_pushbutton,'enable','off');
set(handles.infected_cluster_pushbutton,'enable','off');
set(handles.Feature_Extraction_pushbutton,'enable','off');
set(handles.Neural_Network_pushbutton,'enable','off');


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PlantDisease wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PlantDisease_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in select_image_pushbutton.
function select_image_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to select_image_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rgbImage name
[filename pathname]=uigetfile('*.jpg','Select An Image');
rgbImage=imread([pathname filename]);
name=filename;
arrayfun(@cla,findall(0,'type','axes'));

axes(handles.axes1);
imshow(rgbImage);
axis equal;axis off;

set(handles.color_transformation_structure_pushbutton,'enable','off');
set(handles.RGB_to_HSI_pushbutton,'enable','off');
set(handles.kmeans_pushbutton,'enable','off');
set(handles.masking_pushbutton,'enable','off');
set(handles.infected_cluster_pushbutton,'enable','off');
set(handles.Feature_Extraction_pushbutton,'enable','off');
set(handles.Neural_Network_pushbutton,'enable','off');

preProcessing(rgbImage);

set(handles.color_transformation_structure_pushbutton,'enable','on');


% --- Executes on button press in color_transformation_structure_pushbutton.
function color_transformation_structure_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to color_transformation_structure_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Cmd=1;
colorTransformationStructure(Cmd)
Cmd=2;
colorTransformationStructure(Cmd)

set(handles.color_transformation_structure_pushbutton,'enable','off');
set(handles.RGB_to_HSI_pushbutton,'enable','on');


% --- Executes on button press in RGB_to_HSI_pushbutton.
function RGB_to_HSI_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to RGB_to_HSI_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


clc
global rgbImage
[F,hsi,C,H,S,I] = rgbtohsi(rgbImage);
fontSize = 15;
figure(3), 
subplot(2,3,1), imshow(F),title('RGB Image', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(2,3,2), imshow(hsi),title('HSI Image', 'FontSize', fontSize);
subplot(2,3,3), imshow(C),title('HSI Equalized', 'FontSize', fontSize);
subplot(2,3,4), imshow(H),title('Hue Image', 'FontSize', fontSize);
subplot(2,3,5), imshow(S),title('Saturation Image', 'FontSize', fontSize);
subplot(2,3,6), imshow(I),title('Intensity Image', 'FontSize', fontSize);

set(handles.RGB_to_HSI_pushbutton,'enable','off');
set(handles.kmeans_pushbutton,'enable','on');

% --- Executes on button press in kmeans_pushbutton.
function kmeans_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to kmeans_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rgbImage pixel_labels
cform = makecform('srgb2lab');
lab_he = applycform(rgbImage,cform);
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);
nColors = 4;
[cluster_idx cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean','Replicates',3);
pixel_labels = reshape(cluster_idx,nrows,ncols);
segmented_images = cell(1,3);
rgb_label = repmat(pixel_labels,[1 1 3]);

for k = 1:nColors
    color = rgbImage;
    color(rgb_label ~= k) = 0;
    segmented_images{k} = color;
end

fontSize = 15;
figure(4),
subplot(2, 3, 1);imshow(rgbImage), title('RGB image', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(2, 3, 2);imshow(segmented_images{1}), title('objects in cluster 1', 'FontSize', fontSize);
subplot(2, 3, 3);imshow(segmented_images{2}), title('objects in cluster 2', 'FontSize', fontSize);
subplot(2, 3, 4);imshow(segmented_images{3}), title('objects in cluster 3', 'FontSize', fontSize);
subplot(2, 3, 5);imshow(segmented_images{4}), title('objects in cluster 4', 'FontSize', fontSize);
subplot(2, 3, 6);imshow(pixel_labels,[]), title('Image labeled by cluster index', 'FontSize', fontSize);

set(handles.kmeans_pushbutton,'enable','off');
set(handles.masking_pushbutton,'enable','on');


% --- Executes on button press in masking_pushbutton.
function masking_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to masking_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rgbImage maskedRgbImage greenThreshold
% Get the dimensions of the image.  numberOfColorBands should be = 3.
[rows columns numberOfColorBands] = size(rgbImage);
fontSize = 15;
grayImage = rgb2gray(rgbImage);
I=grayImage;
level = graythresh(I); 
BW = im2bw(I,level);
figure(5),
subplot(3, 3, 1), imshow(rgbImage, []), title('Original Color Image', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(3, 3, 2), imshow(grayImage), title('Gray Image', 'FontSize', fontSize);
subplot(3, 3, 3), imshow(BW), title('Binarized Image', 'FontSize', fontSize);
% Display the original color image.
subplot(3, 3, 4), imshow(rgbImage, []);
title('Noise Free RGB Image', 'FontSize', fontSize);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

% Extract the individual red, green, and blue color channels.
redChannel = rgbImage(:, :, 1);
greenChannel = rgbImage(:, :, 2);
blueChannel = rgbImage(:, :, 3);
% Display
subplot(3, 3, 5), imshow(redChannel, []);
title('Red Channel Image', 'FontSize', fontSize);
subplot(3, 3, 6), imshow(greenChannel, []);
title('Green Channel Image', 'FontSize', fontSize);
subplot(3, 3, 7), imshow(blueChannel, []);
title('Blue Channel Image', 'FontSize', fontSize);

% Get a mask based on the green threshold
greenThreshold = level*255;
mask = greenChannel < greenThreshold;
subplot(3, 3, 8), imshow(mask, []);
title('Green Pixels Masked', 'FontSize', fontSize);

maskedRed = redChannel; % Initialize
maskedGreen = greenChannel; % Initialize
maskedBlue = blueChannel; % Initialize
% Mask 
maskedRed(mask) = 0; % Initialize
maskedGreen(mask) = 255; % Initialize
maskedBlue(mask) = 0; % Initialize

maskedRgbImage = cat(3, maskedRed, maskedGreen, maskedBlue);
subplot(3, 3, 9), imshow(maskedRgbImage);
title('Masked RGB Image', 'FontSize', fontSize);

set(handles.masking_pushbutton,'enable','off');
set(handles.infected_cluster_pushbutton,'enable','on');


% --- Executes on button press in infected_cluster_pushbutton.
function infected_cluster_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to infected_cluster_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global maskedRgbImage H S
[F,hsi,C,H,S,I] = rgbtohsi(maskedRgbImage);
fontSize = 15;
figure(6), 
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1), imshow(H),title('Hue Image of Infected Cluster', 'FontSize', fontSize);
subplot(1,3,2), imshow(S),title('Saturation Image of Infected Cluster', 'FontSize', fontSize);
subplot(1,3,3), imshow(I),title('Intensity Image of Infected Cluster', 'FontSize', fontSize);

set(handles.infected_cluster_pushbutton,'enable','off');
set(handles.Feature_Extraction_pushbutton,'enable','on');


% --- Executes on button press in Feature_Extraction_pushbutton.
function Feature_Extraction_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Feature_Extraction_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc
global H S
SGDM_H = graycomatrix(H,'NumLevels',8,'offset', [0 1], 'Symmetric', true);
display(SGDM_H);
GLCM_H = graycomatrix(H, 'offset', [0 1]);
stats_H = Texture_Features(GLCM_H,0);
display(stats_H);

SGDM_S = graycomatrix(S,'NumLevels',8,'offset', [0 1], 'Symmetric', true);
display(SGDM_S);
GLCM_S = graycomatrix(S, 'offset', [0 1]);
stats_S = Texture_Features(GLCM_S,0);
display(stats_S);

set(handles.Feature_Extraction_pushbutton,'enable','off');
set(handles.Neural_Network_pushbutton,'enable','on');


% --- Executes on button press in Neural_Network_pushbutton.
function Neural_Network_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Neural_Network_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global greenThreshold pixel_labels

create_fit_net();
figure(7);
P = [-3.0 +2.0];
T = [+0.4 +0.8];
wv = -4:0.4:4;
bv = -4:0.4:4;
es = errsurf(P,T,wv,bv,'logsig');
plotes(wv,bv,es,[60 30]);
warning('off','all');
[w,b] = initff(P,T,'logsig');
df = 5; 
me = 100; 
eg = 0.01; 
lr = 2; 
[w,b,ep,tr] = tbp1(w,b,'logsig',P,T,[df me eg lr],wv,bv,es,[60 30]);

load PlantDisease

IdentifyingFeatureValue = IdentifyingFeature;
DiseaseName = Disease;

n = length(IdentifyingFeatureValue);
temp = greenThreshold;


j=0;
value=0;
for i=1:1:n
    if(temp==IdentifyingFeatureValue(i))        
        value=1;
        j=i;
        break;
    end
end

if(value==1)
    DiseaseIdentified = DiseaseName(j);
    axes(handles.axes2);
    imshow(pixel_labels,[]);
    axis equal;axis off;
    message = sprintf('Disease identified as %s', DiseaseIdentified{1});
    uiwait(warndlg(message));     
else
    msgbox('Plant Disease Not Identified');
end

set(handles.Neural_Network_pushbutton,'enable','off');