function varargout = multifocus(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @multifocus_OpeningFcn, ...
                   'gui_OutputFcn',  @multifocus_OutputFcn, ...
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

function multifocus_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = multifocus_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function pushbutton1_Callback(hObject, eventdata, handles)
disp('Choose the first image:')
[filename, pathname]= uigetfile({'*.jpg;*.png;*.tif'},'Choose the first image');
path=fullfile(pathname, filename);
global image1
image1=imread(path);
disp('The first image is succesfully selected')
axes(handles.axes1);
imshow(image1);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Choose the second image:')
[filename, pathname]= uigetfile({'*.jpg;*.png;*.tif'},'Choose the second image');
path=fullfile(pathname, filename);
global image1
global image2
image2=imread(path);
disp('The second image is succesfully selected')


axes(handles.axes2);
imshow(image2);
if size(image1,3) == 3 %Does the image have 3 channel?
    image1 = rgb2gray(image1);
end
if size(image2,3) == 3
    image2 = rgb2gray(image2);
end

if size(image1) ~= size(image2)	%Are the images' size same?
    disp('Images have different sizes!')
end


function popupmenu3_Callback(hObject, eventdata, handles)
contents=cellstr(get(hObject,'String'));
popChoice= contents{get(hObject,'Value')};

switch popChoice
    
    
    case 'DCT+VARIANCE'
    uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal'));
    
    global image1
    global image2

 

%Unsharp Filter
%**************************************************
h = fspecial('unsharp');                                                                                                                                                                                                                                       
image1_sharp=imfilter(image1,h,'symmetric');
image2_sharp=imfilter(image2,h,'symmetric');
%**************************************************

% Specify the image's size
%**************************************************
[m,n] = size(image1);
FusedDCTSharp = zeros(m,n);
FusedDCTSharp_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	% for CV
%**************************************************

% Level shifting for DCT
%**************************************************
image1 = double(image1)-128;
image2 = double(image2)-128;
image1_sharp = double(image1_sharp)-128;
image2_sharp = double(image2_sharp)-128;
%**************************************************

% Dividing the image into 8*8 pieces
for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
        image1s = image1(8*i-7:8*i,8*j-7:8*j);
        image2s = image2(8*i-7:8*i,8*j-7:8*j);
        image1Sub = image1_sharp(8*i-7:8*i,8*j-7:8*j);
        image2Sub = image2_sharp(8*i-7:8*i,8*j-7:8*j);
        
        % DCT
%**************************************************
        image1sdct = dct2(image1s);
        image2sdct = dct2(image2s);
        image1SubDct = dct2(image1Sub);
        image2SubDct = dct2(image2Sub);
%**************************************************
        
        %NTC
%**************************************************
        image1Norm = image1SubDct ./ 8;
        image2Norm = image2SubDct ./ 8;
%**************************************************
        
        % Mean
%**************************************************
        image1Mean = image1Norm(1,1);
        image2Mean = image2Norm(1,1);
%**************************************************
        
        % Variance
%**************************************************
        image1Var = sum(sum(image1Norm.^2)) - image1Mean.^2;
        image2Var = sum(sum(image2Norm.^2)) - image2Mean.^2;
%**************************************************
        

         t=60;
        if image1Var > (image2Var+t)
            dctSub = image1sdct;
            
       
        end
        if image1Var < (image2Var-t)
            dctSub = image2sdct;
        end
% 
        if image1Var < (image2Var+t)&& image1Var > (image2Var-t)
            
            dctSub = (image1sdct+image2sdct)./2;

        end
        

        % inverse DCT
        FusedDCTSharp(8*i-7:8*i,8*j-7:8*j) = idct2(dctSub);	% DCT+Variance method
        
    end
end


    
% inverse shift 
%**************************************************
image1 = uint8(double(image1)+128);
image2 = uint8(double(image2)+128);
FusedDCTSharp = uint8(double(FusedDCTSharp)+128);
%**************************************************



figure, imshow(FusedDCTSharp), title('DCT+VARIANCE');

 case 'DCT+VARIANCE+CV'

 uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal'));
    
    global image1
    global image2

%Unsharp Filter
%**************************************************
h = fspecial('unsharp');                                                                                                                                                                                                                                       
image1_sharp=imfilter(image1,h,'symmetric');
image2_sharp=imfilter(image2,h,'symmetric');
%**************************************************

% Specify the image's size
%**************************************************
[m,n] = size(image1);
FusedDCTSharp = zeros(m,n);
FusedDCTSharp_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	% cv i�in
%**************************************************

% Level shifting for DCT
%**************************************************
image1 = double(image1)-128;
image2 = double(image2)-128;
image1_sharp = double(image1_sharp)-128;
image2_sharp = double(image2_sharp)-128;
%**************************************************

% Dividing the image into 8*8 pieces
for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
       image1s = image1(8*i-7:8*i,8*j-7:8*j);
        image2s = image2(8*i-7:8*i,8*j-7:8*j);
        image1Sub = image1_sharp(8*i-7:8*i,8*j-7:8*j);
        image2Sub = image2_sharp(8*i-7:8*i,8*j-7:8*j);
        
        
        % DCT
        %**************************************************
        image1sdct = dct2(image1s);
        image2sdct = dct2(image2s);
        image1SubDct = dct2(image1Sub);
        image2SubDct = dct2(image2Sub);
        %**************************************************
        
        % NTC 
        %**************************************************
        image1Norm = image1SubDct ./ 8;
        image2Norm = image2SubDct ./ 8;
        %**************************************************
        
        % Mean 
        %**************************************************
        image1Mean = image1Norm(1,1);
        image2Mean = image2Norm(1,1);
        %**************************************************
        
        % Variance
        %**************************************************
        image1Var = sum(sum(image1Norm.^2)) - image1Mean.^2;
        image2Var = sum(sum(image2Norm.^2)) - image2Mean.^2;
        %**************************************************
        

          t=60;
        if image1Var > (image2Var+t)
            dctSub = image1sdct;
            Map(i,j) =-1;	% CV 
       
        end
        if  image1Var < (image2Var-t)
            dctSub = image2sdct;
            Map(i,j) = +1;    % CV
        end
 
        if image1Var < (image2Var+t)&& image1Var > (image2Var-t)
            Map(i,j)=0;
         dctSub = (image1sdct+image2sdct)./2;

        end
        

        % Inverse DCT
        FusedDCTSharp(8*i-7:8*i,8*j-7:8*j) = idct2(dctSub);	% DCT+Variance 
    end
end

% CV
fi = ones(3)/9;	% Filter kernel 7*7

cvMapFiltered = imfilter(Map, fi,'symmetric');	% Filtered index map

for i = 1:m/8
    for j = 1:n/8
        % DCT+Variance+CV method
       

       
     FusedDCTSharp_CV(8*i-7:8*i,8*j-7:8*j) =(((1-cvMapFiltered(i,j))/2)*image1(8*i-7:8*i,8*j-7:8*j))+(((1+cvMapFiltered(i,j))/2)*image2(8*i-7:8*i,8*j-7:8*j));


   
    end
end

for i = 1:m/8
    for j = 1:n/8
        % DCT+Variance+CV method
        if cvMapFiltered(i,j) < -0.06
            FusedDCTSharp_CV(8*i-7:8*i,8*j-7:8*j) = image1(8*i-7:8*i,8*j-7:8*j);

        end
        if cvMapFiltered(i,j) > 0.06
            FusedDCTSharp_CV(8*i-7:8*i,8*j-7:8*j) = image2(8*i-7:8*i,8*j-7:8*j);
       
        end
        
    end
end

    
% inverse shift 
%**************************************************
image1 = uint8(double(image1)+128);
image2 = uint8(double(image2)+128);
%**************************************************

FusedDCTSharp_CV = uint8(double(FusedDCTSharp_CV)+128);
figure, imshow(FusedDCTSharp_CV), title('DCT+VARINACE+CV');


    case 'DCT+CORR'
     global image1
     global image2

    uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal'));
    
     


% low pass

%**********************************************************
image1low=imfilter(image1,fspecial('average',5),'symmetric');
image2low=imfilter(image2,fspecial('average',5),'symmetric');
%**********************************************************

% Get the image's size

[m,n] = size(image1);

% Two new images that contain zeros

FusedDCTCorr = zeros(m,n);
FusedDCTCorr_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	%for CV

% Level shifting
%**************************************************
image1 = double(image1)-128;
image2 = double(image2)-128;
image1low = double(image1low)-128;
image2low = double(image2low)-128;
%**************************************************


for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
        image1_Block = image1(8*i-7:8*i,8*j-7:8*j);
        image2_Block = image2(8*i-7:8*i,8*j-7:8*j);
        image1_Block_Blurred = image1low(8*i-7:8*i,8*j-7:8*j);
        image2_Block_Blurred = image2low(8*i-7:8*i,8*j-7:8*j);
        
        % DCT 
        image1_Block_DCT = dct2(image1_Block);
        image2_Block_DCT = dct2(image2_Block);
        image1_Block_Blurred_DCT = dct2(image1_Block_Blurred);
        image2_Block_Blurred_DCT = dct2(image2_Block_Blurred);
        
        % NTC
        image1Norm = image1_Block_DCT ./ 8;
        image2Norm = image2_Block_DCT ./ 8;
        image1Normlow = image1_Block_Blurred_DCT ./ 8;
        image2Normlow = image2_Block_Blurred_DCT ./ 8;
        
        % Mean 
        image1ave = mean(mean(image1Norm));
        image2ave = mean(mean(image2Norm));
        image1avelow =  mean(mean(image1Normlow));
        image2avelow = mean(mean(image2Normlow));

% CC     
%**************************************************
a = image1_Block_DCT - image1ave;
b = image1_Block_Blurred_DCT - image1avelow;
image1cor = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));
%**************************************************

% CC  
%**************************************************
a = image2_Block_DCT - image2ave;
b = image2_Block_Blurred_DCT - image2avelow;
image2cor = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));
%**************************************************

 % Fusion 
        if image1cor > image2cor
            dctCorrBlock = image2_Block_DCT;
             
        end
        if image1cor <= image2cor
            dctCorrBlock = image1_Block_DCT;
            
        end
        
        % Inverse DCT 
        % DCT+Corr Method
        FusedDCTCorr(8*i-7:8*i,8*j-7:8*j) = idct2(dctCorrBlock);
        
    end
end


% Inverse level shifting 
%**************************************************
image1 = uint8(double(image1)+128);
image2 = uint8(double(image2)+128);
%**************************************************
FusedDCTCorr = uint8(double(FusedDCTCorr)+128);
figure, imshow(FusedDCTCorr), title('"DCT+Corr"');

  
case 'DCT+CORR+CV'
     global image1
     global image2
 uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal'));

   
%**************************************************     
image1low=imfilter(image1,fspecial('average',5),'symmetric');
image2low=imfilter(image2,fspecial('average',5),'symmetric');
%**************************************************

% Sizes
[m,n] = size(image1);
% Two new images that contain zeros
FusedDCTCorr = zeros(m,n);
FusedDCTCorr_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	

% Level shifting
%**************************************************
image1 = double(image1)-128;
image2 = double(image2)-128;
image1low = double(image1low)-128;
image2low = double(image2low)-128;
%**************************************************

%  8*8 
for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
        image1_Block = image1(8*i-7:8*i,8*j-7:8*j);
        image2_Block = image2(8*i-7:8*i,8*j-7:8*j);
        image1_Block_Blurred = image1low(8*i-7:8*i,8*j-7:8*j);
        image2_Block_Blurred = image2low(8*i-7:8*i,8*j-7:8*j);
        
        % DCT
        %**************************************************
        image1_Block_DCT = dct2(image1_Block);
        image2_Block_DCT = dct2(image2_Block);
        image1_Block_Blurred_DCT = dct2(image1_Block_Blurred);
        image2_Block_Blurred_DCT = dct2(image2_Block_Blurred);
        %**************************************************
        
        % NTC
        %**************************************************
        image1Norm = image1_Block_DCT ./ 8;
        image2Norm = image2_Block_DCT ./ 8;
        image1Normlow = image1_Block_Blurred_DCT ./ 8;
        image2Normlow = image2_Block_Blurred_DCT ./ 8;
        %**************************************************
        
        % Mean 
        %**************************************************
        image1ave = mean(mean(image1Norm));
        image2ave = mean(mean(image2Norm));
        image1avelow =  mean(mean(image1Normlow));
        image2avelow = mean(mean(image2Normlow));
        %**************************************************

%CC
%**************************************************
a = image1_Block_DCT - image1ave;
b = image1_Block_Blurred_DCT - image1avelow;
image1cor = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));
%**************************************************

%CC   
%**************************************************
a = image2_Block_DCT - image2ave;
b = image2_Block_Blurred_DCT - image2avelow;
image2cor = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));
%**************************************************

 % Fusion 
        if image1cor > image2cor
            dctCorrBlock = image2_Block_DCT;
            Map(i,j) =+1;	%CV
        end
        if image1cor <= image2cor
            dctCorrBlock = image1_Block_DCT;
            Map(i,j) = -1;    % CV
        end
        
        % Inverse DCT
        % DCT+Corr 
        FusedDCTCorr(8*i-7:8*i,8*j-7:8*j) = idct2(dctCorrBlock);
        
    end
end

% (CV)
%**************************************************
Filter=fspecial('average',3);

Map_Filtered = imfilter(Map, Filter,'symmetric');
%**************************************************	

% CV process
for i = 1:m/8
    for j = 1:n/8
        
        if Map_Filtered(i,j) < 0
            FusedDCTCorr_CV(8*i-7:8*i,8*j-7:8*j) = image1(8*i-7:8*i,8*j-7:8*j);
        else
            FusedDCTCorr_CV(8*i-7:8*i,8*j-7:8*j) = image2(8*i-7:8*i,8*j-7:8*j);
        end
        
    end
end

% inverse shifting 
%**************************************************
image1 = uint8(double(image1)+128);
image2 = uint8(double(image2)+128);
%**************************************************

FusedDCTCorr_CV = uint8(double(FusedDCTCorr_CV)+128);

figure, imshow(FusedDCTCorr_CV), title('"DCT+Corr+CV" ');

     

     case 'DCT+CORR_ENG'
        global image1
    global image2
    uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal')); 


[m,n] = size(image1);
FusedDCT = zeros(m,n);
FusedDCT_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	


image1 = double(image1)-128;
image2 = double(image2)-128;

x=0.0751;
y=0.1238;
z=0.2042;
C=dctmtx(8);
t=[

     y     x     0     0     0     0     0     0
     x     y     x     0     0     0     0     0
     0     x     y     x     0     0     0     0
     0     0     x     y     x     0     0     0
     0     0     0     x     y     x     0     0
     0     0     0     0     x     y     x     0
     0     0     0     0     0     x     y     x
     0     0     0     0     0     0     x     y];
 s=[

     z     y     0     0     0     0     0     0
     y     z     y     0     0     0     0     0
     0     y     z     y     0     0     0     0
     0     0     y     z     y     0     0     0
     0     0     0     y     z     y     0     0
     0     0     0     0     y     z     y     0
     0     0     0     0     0     y     z     y
     0     0     0     0     0     0     y     z];
 u=[

     x+2*y     0     0     0     0     0     0     0
     0         0     0     0     0     0     0     0
     0         0     0     0     0     0     0     0
     0         0     0     0     0     0     0     0
     0         0     0     0     0     0     0     0
     0         0     0     0     0     0     0     0
     0         0     0     0     0     0     0     0
     0         0     0     0     0     0     0     x+2*y ];
 
v=[

     0     x     0     0     0     0     0     0
     x     y     x     0     0     0     0     0
     0     x     y     x     0     0     0     0
     0     0     x     y     x     0     0     0
     0     0     0     x     y     x     0     0
     0     0     0     0     x     y     x     0
     0     0     0     0     0     x     y     x
     0     0     0     0     0     0     x     0];
 
  lu =[

     0     1     0     0     0     0     0     0
     1     0     1     0     0     0     0     0
     0     1     0     1     0     0     0     0
     0     0     1     0     1     0     0     0
     0     0     0     1     0     1     0     0
     0     0     0     0     1     0     1     0
     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     1     0];
   q =[

     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1];
 
 LU=C*lu*C';
 T=C*t*C';
 S=C*s*C';
 U=C*u*C';
 V=C*v*C';
 Q=C*q*C';
 



threshold1=0;
% 8*8 
for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
        image1_Block  = image1(8*i-7:8*i,8*j-7:8*j);
        image2_Block  = image2(8*i-7:8*i,8*j-7:8*j);
        
        % Compute the 2-D DCT of 8*8 Blocks
        image1_Block_DCT = C*image1_Block*C';
        image2_Block_DCT = C*image2_Block*C';
        
        image1SubDct_LOW=(LU*image1_Block_DCT*T)+(image1_Block_DCT*S)+((Q*image1_Block_DCT*U)+(V*image1_Block_DCT*Q)+(Q*image1_Block_DCT*V));
        image2SubDct_LOW=(LU*image2_Block_DCT*T)+(image2_Block_DCT*S)+((Q*image2_Block_DCT*U)+(V*image2_Block_DCT*Q)+(Q*image2_Block_DCT*V));

        PimA=image1_Block_DCT-mean2(image1_Block_DCT);
        PimB=image2_Block_DCT-mean2(image2_Block_DCT);
        PimA_Low=image1SubDct_LOW-mean2(image1SubDct_LOW);
        PimB_Low=image2SubDct_LOW-mean2(image2SubDct_LOW);
        
        cor1= sum(sum(PimA.*PimA_Low))/sqrt(sum(sum(PimA.*PimA))*sum(sum(PimA_Low.*PimA_Low)));
        cor2= sum(sum(PimB.*PimB_Low))/sqrt(sum(sum(PimB.*PimB))*sum(sum(PimB_Low.*PimB_Low)));
       
        energy_A=sum(sum(image1_Block_DCT.^2));
        energy_B=sum(sum(image2_Block_DCT.^2));
        energy_A_Low=sum(sum(image1SubDct_LOW.^2));
        energy_B_Low=sum(sum(image2SubDct_LOW.^2));
        
        corr_eng_1=energy_A*(1-cor1)*energy_A_Low;
        corr_eng_2=energy_B*(1-cor2)*energy_B_Low;
        
        
        z=corr_eng_1;
        zz=corr_eng_2;


        if z>=zz
           dctBlock = image1_Block_DCT;
            Map(i,j) = -1;	% CV

        end
        if z<zz
            dctBlock = image2_Block_DCT;
            Map(i,j) = +1;    % CV



        end
        if z<zz+threshold1 && z>zz-threshold1
            dctBlock = (image2_Block_DCT+image2_Block_DCT)./2;
            Map(i,j) =0 ;
        end
        
        % Inverse  DCT
        FusedDCT(8*i-7:8*i,8*j-7:8*j) = C'*dctBlock*C;	% DCT+Corr_Eng method
       
    end
end



threshold2=0.00;

% Inverse level shifting 
image1 = uint8(double(image1)+128);
image2 = uint8(double(image2)+128);

FusedDCT = uint8(double(FusedDCT)+128);


figure, imshow(FusedDCT), title('"DCT+Corr_Eng"');



case 'DCT+CORR_ENG+CV'
    
    global image1
    global image2
    uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal')); 


[m,n] = size(image1);
FusedDCT = zeros(m,n);
FusedDCT_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	


image1 = double(image1)-128;
image2 = double(image2)-128;

x=0.0751;
y=0.1238;
z=0.2042;
C=dctmtx(8);
t=[

     y     x     0     0     0     0     0     0
     x     y     x     0     0     0     0     0
     0     x     y     x     0     0     0     0
     0     0     x     y     x     0     0     0
     0     0     0     x     y     x     0     0
     0     0     0     0     x     y     x     0
     0     0     0     0     0     x     y     x
     0     0     0     0     0     0     x     y];
 s=[

     z     y     0     0     0     0     0     0
     y     z     y     0     0     0     0     0
     0     y     z     y     0     0     0     0
     0     0     y     z     y     0     0     0
     0     0     0     y     z     y     0     0
     0     0     0     0     y     z     y     0
     0     0     0     0     0     y     z     y
     0     0     0     0     0     0     y     z];
 u=[

     x+2*y     0     0     0     0     0     0     0
     0         0     0     0     0     0     0     0
     0         0     0     0     0     0     0     0
     0         0     0     0     0     0     0     0
     0         0     0     0     0     0     0     0
     0         0     0     0     0     0     0     0
     0         0     0     0     0     0     0     0
     0         0     0     0     0     0     0     x+2*y ];
 
v=[

     0     x     0     0     0     0     0     0
     x     y     x     0     0     0     0     0
     0     x     y     x     0     0     0     0
     0     0     x     y     x     0     0     0
     0     0     0     x     y     x     0     0
     0     0     0     0     x     y     x     0
     0     0     0     0     0     x     y     x
     0     0     0     0     0     0     x     0];
 
  lu =[

     0     1     0     0     0     0     0     0
     1     0     1     0     0     0     0     0
     0     1     0     1     0     0     0     0
     0     0     1     0     1     0     0     0
     0     0     0     1     0     1     0     0
     0     0     0     0     1     0     1     0
     0     0     0     0     0     1     0     1
     0     0     0     0     0     0     1     0];
   q =[

     1     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1];
 
 LU=C*lu*C';
 T=C*t*C';
 S=C*s*C';
 U=C*u*C';
 V=C*v*C';
 Q=C*q*C';
 



threshold1=0;
% 8*8 
for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
        image1_Block  = image1(8*i-7:8*i,8*j-7:8*j);
        image2_Block  = image2(8*i-7:8*i,8*j-7:8*j);
        
        % Compute the 2-D DCT of 8*8 Blocks
        image1_Block_DCT = C*image1_Block*C';
        image2_Block_DCT = C*image2_Block*C';
        
        image1SubDct_LOW=(LU*image1_Block_DCT*T)+(image1_Block_DCT*S)+((Q*image1_Block_DCT*U)+(V*image1_Block_DCT*Q)+(Q*image1_Block_DCT*V));
        image2SubDct_LOW=(LU*image2_Block_DCT*T)+(image2_Block_DCT*S)+((Q*image2_Block_DCT*U)+(V*image2_Block_DCT*Q)+(Q*image2_Block_DCT*V));

        PimA=image1_Block_DCT-mean2(image1_Block_DCT);
        PimB=image2_Block_DCT-mean2(image2_Block_DCT);
        PimA_Low=image1SubDct_LOW-mean2(image1SubDct_LOW);
        PimB_Low=image2SubDct_LOW-mean2(image2SubDct_LOW);
        
        cor1= sum(sum(PimA.*PimA_Low))/sqrt(sum(sum(PimA.*PimA))*sum(sum(PimA_Low.*PimA_Low)));
        cor2= sum(sum(PimB.*PimB_Low))/sqrt(sum(sum(PimB.*PimB))*sum(sum(PimB_Low.*PimB_Low)));
       
        energy_A=sum(sum(image1_Block_DCT.^2));
        energy_B=sum(sum(image2_Block_DCT.^2));
        energy_A_Low=sum(sum(image1SubDct_LOW.^2));
        energy_B_Low=sum(sum(image2SubDct_LOW.^2));
        
        corr_eng_1=energy_A*(1-cor1)*energy_A_Low;
        corr_eng_2=energy_B*(1-cor2)*energy_B_Low;
        
        
        z=corr_eng_1;
        zz=corr_eng_2;


        if z>=zz
           dctBlock = image1_Block_DCT;
            Map(i,j) = -1;	% CV

        end
        if z<zz
            dctBlock = image2_Block_DCT;
            Map(i,j) = +1;    % CV



        end
        if z<zz+threshold1 && z>zz-threshold1
            dctBlock = (image2_Block_DCT+image2_Block_DCT)./2;
            Map(i,j) =0 ;
        end
        
        % Inverse  DCT
        FusedDCT(8*i-7:8*i,8*j-7:8*j) = C'*dctBlock*C;	% DCT+Corr_Eng method
       
    end
end

% (CV) 
fi=fspecial('average',3);
Map_Filtered = imfilter(Map, fi,'symmetric');	% Filtered index map

threshold2=0.00;
for i = 1:m/8
    for j = 1:n/8
        % DCT+Variance+CV method
        if Map_Filtered(i,j) <= -threshold2
            FusedDCT_CV(8*i-7:8*i,8*j-7:8*j) = image1(8*i-7:8*i,8*j-7:8*j);
   
        end
        if Map_Filtered(i,j) > threshold2
            FusedDCT_CV(8*i-7:8*i,8*j-7:8*j) = image2(8*i-7:8*i,8*j-7:8*j);
      
        end
        if Map_Filtered(i,j) > -threshold2 &&  Map_Filtered(i,j) < threshold2
             FusedDCT_CV(8*i-7:8*i,8*j-7:8*j) = (image1(8*i-7:8*i,8*j-7:8*j)+image2(8*i-7:8*i,8*j-7:8*j))./2;
        end
        end
end

% Inverse level shifting 
image1 = uint8(double(image1)+128);
image2 = uint8(double(image2)+128);

FusedDCT_CV = uint8(double(FusedDCT_CV)+128);


figure, imshow(FusedDCT_CV), title('"DCT+Corr_Eng+CV"');

%DCT+SVD
    case 'DCT+SVD'
        
uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal')); 
%**************************************************
global image1
global image2
[m,n] = size(image1);
FusedDCT = zeros(m,n);
FusedDCT_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	
%**************************************************

% Level shifting
%**************************************************
image1 = double(image1)-128;
image2 = double(image2)-128;
%**************************************************

% 8*8 
for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
        image1_Block = image1(8*i-7:8*i,8*j-7:8*j);
        image2_Block = image2(8*i-7:8*i,8*j-7:8*j);
        %  DCT 
        image1_Block_DCT = dct2(image1_Block);
        image2_Block_DCT = dct2(image2_Block);
        sigma1=svd(image1_Block_DCT);
        sigma2=svd(image2_Block_DCT);
        
        x1=sigma1(1)*sigma1(2)*sigma1(3)*sigma1(4)*sigma1(5);
        x2=sigma2(1)*sigma2(2)*sigma2(3)*sigma2(4)*sigma2(5);
        
         % Fusion 
        if x1 > x2
            dctBlock = image1_Block_DCT;
             
        else
            dctBlock = image2_Block_DCT;
             
        end
        
        % Compute the 2-D inverse DCT of 8*8 Blocks and construct fused image
        % DCT+SVD Method
        FusedDCT(8*i-7:8*i,8*j-7:8*j) = idct2(dctBlock);
        
    end
end
% inverse shifting 
%**************************************************
image1 = uint8(double(image1)+128);
image2 = uint8(double(image2)+128);
FusedDCT = uint8(double(FusedDCT)+128);
%**************************************************




figure, imshow(FusedDCT), title('"DCT+SVD" fusion result');

%DCT+SVD+CV    
    case 'DCT+SVD+CV'
%uyar�
uiwait(msgbox('Bu i�lem biraz zaman alabilir.','Uyar�!','modal')); 
        
%**************************************************
global image1
global image2
        
[m,n] = size(image1);
FusedDCT = zeros(m,n);
FusedDCT_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	
%**************************************************

% Level shifting
%**************************************************
image1 = double(image1)-128;
image2 = double(image2)-128;
%**************************************************

% 8*8 
for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
        image1_Block = image1(8*i-7:8*i,8*j-7:8*j);
        image2_Block = image2(8*i-7:8*i,8*j-7:8*j);
        %  DCT  
        image1_Block_DCT = dct2(image1_Block);
        image2_Block_DCT = dct2(image2_Block);
        sigma1=svd(image1_Block_DCT);
        sigma2=svd(image2_Block_DCT);
        
        x1=sigma1(1)*sigma1(2)*sigma1(3)*sigma1(4)*sigma1(5);
        x2=sigma2(1)*sigma2(2)*sigma2(3)*sigma2(4)*sigma2(5);
        
         % Fusion
        if x1 > x2
            dctBlock = image1_Block_DCT;
            Map(i,j) =+1;	% CV
        else
            dctBlock = image2_Block_DCT;
            Map(i,j) = -1;    % CV 
        end
        
        % Inverse DCT 
        % DCT+SVD Method
        FusedDCT(8*i-7:8*i,8*j-7:8*j) = idct2(dctBlock);
        
    end
end

% (CV)(3x3 Averaging Filter)
%**************************************************
Filter=fspecial('average',3);

Map_Filtered = imfilter(Map, Filter,'symmetric');
%**************************************************	

%CV
%**************************************************
for i = 1:m/8
    for j = 1:n/8
        
        if Map_Filtered(i,j) > 0
            FusedDCT_CV(8*i-7:8*i,8*j-7:8*j) = image1(8*i-7:8*i,8*j-7:8*j);
        else
            FusedDCT_CV(8*i-7:8*i,8*j-7:8*j) = image2(8*i-7:8*i,8*j-7:8*j);
        end
        
    end
end
%**************************************************


% inverse shifting 
%**************************************************
image1 = uint8(double(image1)+128);
image2 = uint8(double(image2)+128);
%**************************************************

FusedDCT_CV = uint8(double(FusedDCT_CV)+128);
figure, imshow(FusedDCT_CV), title('"DCT+SVD+CV" ');
end
        

function popupmenu3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
