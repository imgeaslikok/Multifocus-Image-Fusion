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
global resim1
resim1=imread(path);
disp('The first image is succesfully selected')
axes(handles.axes1);
imshow(resim1);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Choose the second image:')
[filename, pathname]= uigetfile({'*.jpg;*.png;*.tif'},'Choose the second image');
path=fullfile(pathname, filename);
global resim1
global resim2
resim2=imread(path);
disp('The second image is succesfully selected')


axes(handles.axes2);
imshow(resim2);
if size(resim1,3) == 3 %Does the image have 3 channel?
    resim1 = rgb2gray(resim1);
end
if size(resim2,3) == 3
    resim2 = rgb2gray(resim2);
end

if size(resim1) ~= size(resim2)	%Are the images' size same?
    disp('Images have different sizes!')
end


function popupmenu3_Callback(hObject, eventdata, handles)
contents=cellstr(get(hObject,'String'));
popChoice= contents{get(hObject,'Value')};

switch popChoice
    
    
    case 'DCT+VARIANCE'
    uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal'));
    
    global resim1
    global resim2

 

%Unsharp Filter
%**************************************************
h = fspecial('unsharp');                                                                                                                                                                                                                                       
resim1_sharp=imfilter(resim1,h,'symmetric');
resim2_sharp=imfilter(resim2,h,'symmetric');
%**************************************************

% Specify the image's size
%**************************************************
[m,n] = size(resim1);
FusedDCTSharp = zeros(m,n);
FusedDCTSharp_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	% for CV
%**************************************************

% Level shifting for DCT
%**************************************************
resim1 = double(resim1)-128;
resim2 = double(resim2)-128;
resim1_sharp = double(resim1_sharp)-128;
resim2_sharp = double(resim2_sharp)-128;
%**************************************************

% Dividing the image into 8*8 pieces
for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
        resim1s = resim1(8*i-7:8*i,8*j-7:8*j);
        resim2s = resim2(8*i-7:8*i,8*j-7:8*j);
        resim1Sub = resim1_sharp(8*i-7:8*i,8*j-7:8*j);
        resim2Sub = resim2_sharp(8*i-7:8*i,8*j-7:8*j);
        
        % DCT
%**************************************************
        resim1sdct = dct2(resim1s);
        resim2sdct = dct2(resim2s);
        resim1SubDct = dct2(resim1Sub);
        resim2SubDct = dct2(resim2Sub);
%**************************************************
        
        %NTC
%**************************************************
        resim1Norm = resim1SubDct ./ 8;
        resim2Norm = resim2SubDct ./ 8;
%**************************************************
        
        % Mean
%**************************************************
        resim1Mean = resim1Norm(1,1);
        resim2Mean = resim2Norm(1,1);
%**************************************************
        
        % Variance
%**************************************************
        resim1Var = sum(sum(resim1Norm.^2)) - resim1Mean.^2;
        resim2Var = sum(sum(resim2Norm.^2)) - resim2Mean.^2;
%**************************************************
        

         t=60;
        if resim1Var > (resim2Var+t)
            dctSub = resim1sdct;
            
       
        end
        if resim1Var < (resim2Var-t)
            dctSub = resim2sdct;
        end
% 
        if resim1Var < (resim2Var+t)&& resim1Var > (resim2Var-t)
            
            dctSub = (resim1sdct+resim2sdct)./2;

        end
        

        % inverse DCT
        FusedDCTSharp(8*i-7:8*i,8*j-7:8*j) = idct2(dctSub);	% DCT+Variance method
        
    end
end


    
% inverse shift 
%**************************************************
resim1 = uint8(double(resim1)+128);
resim2 = uint8(double(resim2)+128);
FusedDCTSharp = uint8(double(FusedDCTSharp)+128);
%**************************************************



figure, imshow(FusedDCTSharp), title('DCT+VARIANCE');

 case 'DCT+VARIANCE+CV'

 uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal'));
    
    global resim1
    global resim2

%Unsharp Filter
%**************************************************
h = fspecial('unsharp');                                                                                                                                                                                                                                       
resim1_sharp=imfilter(resim1,h,'symmetric');
resim2_sharp=imfilter(resim2,h,'symmetric');
%**************************************************

% Specify the image's size
%**************************************************
[m,n] = size(resim1);
FusedDCTSharp = zeros(m,n);
FusedDCTSharp_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	% cv için
%**************************************************

% Level shifting for DCT
%**************************************************
resim1 = double(resim1)-128;
resim2 = double(resim2)-128;
resim1_sharp = double(resim1_sharp)-128;
resim2_sharp = double(resim2_sharp)-128;
%**************************************************

% Dividing the image into 8*8 pieces
for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
       resim1s = resim1(8*i-7:8*i,8*j-7:8*j);
        resim2s = resim2(8*i-7:8*i,8*j-7:8*j);
        resim1Sub = resim1_sharp(8*i-7:8*i,8*j-7:8*j);
        resim2Sub = resim2_sharp(8*i-7:8*i,8*j-7:8*j);
        
        
        % DCT
        %**************************************************
        resim1sdct = dct2(resim1s);
        resim2sdct = dct2(resim2s);
        resim1SubDct = dct2(resim1Sub);
        resim2SubDct = dct2(resim2Sub);
        %**************************************************
        
        % NTC 
        %**************************************************
        resim1Norm = resim1SubDct ./ 8;
        resim2Norm = resim2SubDct ./ 8;
        %**************************************************
        
        % Mean 
        %**************************************************
        resim1Mean = resim1Norm(1,1);
        resim2Mean = resim2Norm(1,1);
        %**************************************************
        
        % Variance
        %**************************************************
        resim1Var = sum(sum(resim1Norm.^2)) - resim1Mean.^2;
        resim2Var = sum(sum(resim2Norm.^2)) - resim2Mean.^2;
        %**************************************************
        

          t=60;
        if resim1Var > (resim2Var+t)
            dctSub = resim1sdct;
            Map(i,j) =-1;	% CV 
       
        end
        if  resim1Var < (resim2Var-t)
            dctSub = resim2sdct;
            Map(i,j) = +1;    % CV
        end
 
        if resim1Var < (resim2Var+t)&& resim1Var > (resim2Var-t)
            Map(i,j)=0;
         dctSub = (resim1sdct+resim2sdct)./2;

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
       

       
     FusedDCTSharp_CV(8*i-7:8*i,8*j-7:8*j) =(((1-cvMapFiltered(i,j))/2)*resim1(8*i-7:8*i,8*j-7:8*j))+(((1+cvMapFiltered(i,j))/2)*resim2(8*i-7:8*i,8*j-7:8*j));


   
    end
end

for i = 1:m/8
    for j = 1:n/8
        % DCT+Variance+CV method
        if cvMapFiltered(i,j) < -0.06
            FusedDCTSharp_CV(8*i-7:8*i,8*j-7:8*j) = resim1(8*i-7:8*i,8*j-7:8*j);

        end
        if cvMapFiltered(i,j) > 0.06
            FusedDCTSharp_CV(8*i-7:8*i,8*j-7:8*j) = resim2(8*i-7:8*i,8*j-7:8*j);
       
        end
        
    end
end

    
% inverse shift 
%**************************************************
resim1 = uint8(double(resim1)+128);
resim2 = uint8(double(resim2)+128);
%**************************************************

FusedDCTSharp_CV = uint8(double(FusedDCTSharp_CV)+128);
figure, imshow(FusedDCTSharp_CV), title('DCT+VARINACE+CV');


    case 'DCT+CORR'
     global resim1
     global resim2

    uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal'));
    
     


% low pass

%**********************************************************
resim1low=imfilter(resim1,fspecial('average',5),'symmetric');
resim2low=imfilter(resim2,fspecial('average',5),'symmetric');
%**********************************************************

% Get the image's size

[m,n] = size(resim1);

% Two new images that contain zeros

FusedDCTCorr = zeros(m,n);
FusedDCTCorr_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	%for CV

% Level shifting
%**************************************************
resim1 = double(resim1)-128;
resim2 = double(resim2)-128;
resim1low = double(resim1low)-128;
resim2low = double(resim2low)-128;
%**************************************************


for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
        resim1_Blok = resim1(8*i-7:8*i,8*j-7:8*j);
        resim2_Blok = resim2(8*i-7:8*i,8*j-7:8*j);
        resim1_Blok_Blurred = resim1low(8*i-7:8*i,8*j-7:8*j);
        resim2_Blok_Blurred = resim2low(8*i-7:8*i,8*j-7:8*j);
        
        % DCT 
        resim1_Blok_DCT = dct2(resim1_Blok);
        resim2_Blok_DCT = dct2(resim2_Blok);
        resim1_Blok_Blurred_DCT = dct2(resim1_Blok_Blurred);
        resim2_Blok_Blurred_DCT = dct2(resim2_Blok_Blurred);
        
        % NTC
        resim1Norm = resim1_Blok_DCT ./ 8;
        resim2Norm = resim2_Blok_DCT ./ 8;
        resim1Normlow = resim1_Blok_Blurred_DCT ./ 8;
        resim2Normlow = resim2_Blok_Blurred_DCT ./ 8;
        
        % Mean 
        resim1ave = mean(mean(resim1Norm));
        resim2ave = mean(mean(resim2Norm));
        resim1avelow =  mean(mean(resim1Normlow));
        resim2avelow = mean(mean(resim2Normlow));

% CC     
%**************************************************
a = resim1_Blok_DCT - resim1ave;
b = resim1_Blok_Blurred_DCT - resim1avelow;
resim1cor = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));
%**************************************************

% CC  
%**************************************************
a = resim2_Blok_DCT - resim2ave;
b = resim2_Blok_Blurred_DCT - resim2avelow;
resim2cor = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));
%**************************************************

 % Fusion 
        if resim1cor > resim2cor
            dctCorrBlok = resim2_Blok_DCT;
             
        end
        if resim1cor <= resim2cor
            dctCorrBlok = resim1_Blok_DCT;
            
        end
        
        % Inverse DCT 
        % DCT+Corr Method
        FusedDCTCorr(8*i-7:8*i,8*j-7:8*j) = idct2(dctCorrBlok);
        
    end
end


% Inverse level shifting 
%**************************************************
resim1 = uint8(double(resim1)+128);
resim2 = uint8(double(resim2)+128);
%**************************************************
FusedDCTCorr = uint8(double(FusedDCTCorr)+128);
figure, imshow(FusedDCTCorr), title('"DCT+Corr"');

  
case 'DCT+CORR+CV'
     global resim1
     global resim2
 uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal'));

   
%**************************************************     
resim1low=imfilter(resim1,fspecial('average',5),'symmetric');
resim2low=imfilter(resim2,fspecial('average',5),'symmetric');
%**************************************************

% Sizes
[m,n] = size(resim1);
% Two new images that contain zeros
FusedDCTCorr = zeros(m,n);
FusedDCTCorr_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	

% Level shifting
%**************************************************
resim1 = double(resim1)-128;
resim2 = double(resim2)-128;
resim1low = double(resim1low)-128;
resim2low = double(resim2low)-128;
%**************************************************

%  8*8 
for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
        resim1_Blok = resim1(8*i-7:8*i,8*j-7:8*j);
        resim2_Blok = resim2(8*i-7:8*i,8*j-7:8*j);
        resim1_Blok_Blurred = resim1low(8*i-7:8*i,8*j-7:8*j);
        resim2_Blok_Blurred = resim2low(8*i-7:8*i,8*j-7:8*j);
        
        % DCT
        %**************************************************
        resim1_Blok_DCT = dct2(resim1_Blok);
        resim2_Blok_DCT = dct2(resim2_Blok);
        resim1_Blok_Blurred_DCT = dct2(resim1_Blok_Blurred);
        resim2_Blok_Blurred_DCT = dct2(resim2_Blok_Blurred);
        %**************************************************
        
        % NTC
        %**************************************************
        resim1Norm = resim1_Blok_DCT ./ 8;
        resim2Norm = resim2_Blok_DCT ./ 8;
        resim1Normlow = resim1_Blok_Blurred_DCT ./ 8;
        resim2Normlow = resim2_Blok_Blurred_DCT ./ 8;
        %**************************************************
        
        % Mean 
        %**************************************************
        resim1ave = mean(mean(resim1Norm));
        resim2ave = mean(mean(resim2Norm));
        resim1avelow =  mean(mean(resim1Normlow));
        resim2avelow = mean(mean(resim2Normlow));
        %**************************************************

%CC
%**************************************************
a = resim1_Blok_DCT - resim1ave;
b = resim1_Blok_Blurred_DCT - resim1avelow;
resim1cor = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));
%**************************************************

%CC   
%**************************************************
a = resim2_Blok_DCT - resim2ave;
b = resim2_Blok_Blurred_DCT - resim2avelow;
resim2cor = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));
%**************************************************

 % Fusion 
        if resim1cor > resim2cor
            dctCorrBlok = resim2_Blok_DCT;
            Map(i,j) =+1;	%CV
        end
        if resim1cor <= resim2cor
            dctCorrBlok = resim1_Blok_DCT;
            Map(i,j) = -1;    % CV
        end
        
        % Inverse DCT
        % DCT+Corr 
        FusedDCTCorr(8*i-7:8*i,8*j-7:8*j) = idct2(dctCorrBlok);
        
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
            FusedDCTCorr_CV(8*i-7:8*i,8*j-7:8*j) = resim1(8*i-7:8*i,8*j-7:8*j);
        else
            FusedDCTCorr_CV(8*i-7:8*i,8*j-7:8*j) = resim2(8*i-7:8*i,8*j-7:8*j);
        end
        
    end
end

% inverse shifting 
%**************************************************
resim1 = uint8(double(resim1)+128);
resim2 = uint8(double(resim2)+128);
%**************************************************

FusedDCTCorr_CV = uint8(double(FusedDCTCorr_CV)+128);

figure, imshow(FusedDCTCorr_CV), title('"DCT+Corr+CV" ');

     

     case 'DCT+CORR_ENG'
        global resim1
    global resim2
    uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal')); 


[m,n] = size(resim1);
FusedDCT = zeros(m,n);
FusedDCT_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	


resim1 = double(resim1)-128;
resim2 = double(resim2)-128;

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
        
        resim1_Blok  = resim1(8*i-7:8*i,8*j-7:8*j);
        resim2_Blok  = resim2(8*i-7:8*i,8*j-7:8*j);
        
        % Compute the 2-D DCT of 8*8 Bloks
        resim1_Blok_DCT = C*resim1_Blok*C';
        resim2_Blok_DCT = C*resim2_Blok*C';
        
        resim1SubDct_LOW=(LU*resim1_Blok_DCT*T)+(resim1_Blok_DCT*S)+((Q*resim1_Blok_DCT*U)+(V*resim1_Blok_DCT*Q)+(Q*resim1_Blok_DCT*V));
        resim2SubDct_LOW=(LU*resim2_Blok_DCT*T)+(resim2_Blok_DCT*S)+((Q*resim2_Blok_DCT*U)+(V*resim2_Blok_DCT*Q)+(Q*resim2_Blok_DCT*V));

        PimA=resim1_Blok_DCT-mean2(resim1_Blok_DCT);
        PimB=resim2_Blok_DCT-mean2(resim2_Blok_DCT);
        PimA_Low=resim1SubDct_LOW-mean2(resim1SubDct_LOW);
        PimB_Low=resim2SubDct_LOW-mean2(resim2SubDct_LOW);
        
        cor1= sum(sum(PimA.*PimA_Low))/sqrt(sum(sum(PimA.*PimA))*sum(sum(PimA_Low.*PimA_Low)));
        cor2= sum(sum(PimB.*PimB_Low))/sqrt(sum(sum(PimB.*PimB))*sum(sum(PimB_Low.*PimB_Low)));
       
        energy_A=sum(sum(resim1_Blok_DCT.^2));
        energy_B=sum(sum(resim2_Blok_DCT.^2));
        energy_A_Low=sum(sum(resim1SubDct_LOW.^2));
        energy_B_Low=sum(sum(resim2SubDct_LOW.^2));
        
        corr_eng_1=energy_A*(1-cor1)*energy_A_Low;
        corr_eng_2=energy_B*(1-cor2)*energy_B_Low;
        
        
        z=corr_eng_1;
        zz=corr_eng_2;


        if z>=zz
           dctBlok = resim1_Blok_DCT;
            Map(i,j) = -1;	% CV

        end
        if z<zz
            dctBlok = resim2_Blok_DCT;
            Map(i,j) = +1;    % CV



        end
        if z<zz+threshold1 && z>zz-threshold1
            dctBlok = (resim2_Blok_DCT+resim2_Blok_DCT)./2;
            Map(i,j) =0 ;
        end
        
        % Inverse  DCT
        FusedDCT(8*i-7:8*i,8*j-7:8*j) = C'*dctBlok*C;	% DCT+Corr_Eng method
       
    end
end



threshold2=0.00;

% Inverse level shifting 
resim1 = uint8(double(resim1)+128);
resim2 = uint8(double(resim2)+128);

FusedDCT = uint8(double(FusedDCT)+128);


figure, imshow(FusedDCT), title('"DCT+Corr_Eng"');



case 'DCT+CORR_ENG+CV'
    
    global resim1
    global resim2
    uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal')); 


[m,n] = size(resim1);
FusedDCT = zeros(m,n);
FusedDCT_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	


resim1 = double(resim1)-128;
resim2 = double(resim2)-128;

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
        
        resim1_Blok  = resim1(8*i-7:8*i,8*j-7:8*j);
        resim2_Blok  = resim2(8*i-7:8*i,8*j-7:8*j);
        
        % Compute the 2-D DCT of 8*8 Bloks
        resim1_Blok_DCT = C*resim1_Blok*C';
        resim2_Blok_DCT = C*resim2_Blok*C';
        
        resim1SubDct_LOW=(LU*resim1_Blok_DCT*T)+(resim1_Blok_DCT*S)+((Q*resim1_Blok_DCT*U)+(V*resim1_Blok_DCT*Q)+(Q*resim1_Blok_DCT*V));
        resim2SubDct_LOW=(LU*resim2_Blok_DCT*T)+(resim2_Blok_DCT*S)+((Q*resim2_Blok_DCT*U)+(V*resim2_Blok_DCT*Q)+(Q*resim2_Blok_DCT*V));

        PimA=resim1_Blok_DCT-mean2(resim1_Blok_DCT);
        PimB=resim2_Blok_DCT-mean2(resim2_Blok_DCT);
        PimA_Low=resim1SubDct_LOW-mean2(resim1SubDct_LOW);
        PimB_Low=resim2SubDct_LOW-mean2(resim2SubDct_LOW);
        
        cor1= sum(sum(PimA.*PimA_Low))/sqrt(sum(sum(PimA.*PimA))*sum(sum(PimA_Low.*PimA_Low)));
        cor2= sum(sum(PimB.*PimB_Low))/sqrt(sum(sum(PimB.*PimB))*sum(sum(PimB_Low.*PimB_Low)));
       
        energy_A=sum(sum(resim1_Blok_DCT.^2));
        energy_B=sum(sum(resim2_Blok_DCT.^2));
        energy_A_Low=sum(sum(resim1SubDct_LOW.^2));
        energy_B_Low=sum(sum(resim2SubDct_LOW.^2));
        
        corr_eng_1=energy_A*(1-cor1)*energy_A_Low;
        corr_eng_2=energy_B*(1-cor2)*energy_B_Low;
        
        
        z=corr_eng_1;
        zz=corr_eng_2;


        if z>=zz
           dctBlok = resim1_Blok_DCT;
            Map(i,j) = -1;	% CV

        end
        if z<zz
            dctBlok = resim2_Blok_DCT;
            Map(i,j) = +1;    % CV



        end
        if z<zz+threshold1 && z>zz-threshold1
            dctBlok = (resim2_Blok_DCT+resim2_Blok_DCT)./2;
            Map(i,j) =0 ;
        end
        
        % Inverse  DCT
        FusedDCT(8*i-7:8*i,8*j-7:8*j) = C'*dctBlok*C;	% DCT+Corr_Eng method
       
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
            FusedDCT_CV(8*i-7:8*i,8*j-7:8*j) = resim1(8*i-7:8*i,8*j-7:8*j);
   
        end
        if Map_Filtered(i,j) > threshold2
            FusedDCT_CV(8*i-7:8*i,8*j-7:8*j) = resim2(8*i-7:8*i,8*j-7:8*j);
      
        end
        if Map_Filtered(i,j) > -threshold2 &&  Map_Filtered(i,j) < threshold2
             FusedDCT_CV(8*i-7:8*i,8*j-7:8*j) = (resim1(8*i-7:8*i,8*j-7:8*j)+resim2(8*i-7:8*i,8*j-7:8*j))./2;
        end
        end
end

% Inverse level shifting 
resim1 = uint8(double(resim1)+128);
resim2 = uint8(double(resim2)+128);

FusedDCT_CV = uint8(double(FusedDCT_CV)+128);


figure, imshow(FusedDCT_CV), title('"DCT+Corr_Eng+CV"');

%DCT+SVD
    case 'DCT+SVD'
        
uiwait(msgbox('This process may take some time, please do not close this window.','Warning','modal')); 
%**************************************************
global resim1
global resim2
[m,n] = size(resim1);
FusedDCT = zeros(m,n);
FusedDCT_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	
%**************************************************

% Level shifting
%**************************************************
resim1 = double(resim1)-128;
resim2 = double(resim2)-128;
%**************************************************

% 8*8 
for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
        resim1_Blok = resim1(8*i-7:8*i,8*j-7:8*j);
        resim2_Blok = resim2(8*i-7:8*i,8*j-7:8*j);
        %  DCT 
        resim1_Blok_DCT = dct2(resim1_Blok);
        resim2_Blok_DCT = dct2(resim2_Blok);
        sigma1=svd(resim1_Blok_DCT);
        sigma2=svd(resim2_Blok_DCT);
        
        x1=sigma1(1)*sigma1(2)*sigma1(3)*sigma1(4)*sigma1(5);
        x2=sigma2(1)*sigma2(2)*sigma2(3)*sigma2(4)*sigma2(5);
        
         % Fusion 
        if x1 > x2
            dctBlok = resim1_Blok_DCT;
             
        else
            dctBlok = resim2_Blok_DCT;
             
        end
        
        % Compute the 2-D inverse DCT of 8*8 Bloks and construct fused image
        % DCT+SVD Method
        FusedDCT(8*i-7:8*i,8*j-7:8*j) = idct2(dctBlok);
        
    end
end
% inverse shifting 
%**************************************************
resim1 = uint8(double(resim1)+128);
resim2 = uint8(double(resim2)+128);
FusedDCT = uint8(double(FusedDCT)+128);
%**************************************************




figure, imshow(FusedDCT), title('"DCT+SVD" fusion result');

%DCT+SVD+CV    
    case 'DCT+SVD+CV'
%uyarý
uiwait(msgbox('Bu iþlem biraz zaman alabilir.','Uyarý!','modal')); 
        
%**************************************************
global resim1
global resim2
        
[m,n] = size(resim1);
FusedDCT = zeros(m,n);
FusedDCT_CV = zeros(m,n);
Map = zeros(floor(m/8),floor(n/8));	
%**************************************************

% Level shifting
%**************************************************
resim1 = double(resim1)-128;
resim2 = double(resim2)-128;
%**************************************************

% 8*8 
for i = 1:floor(m/8)
    for j = 1:floor(n/8)
        
        resim1_Blok = resim1(8*i-7:8*i,8*j-7:8*j);
        resim2_Blok = resim2(8*i-7:8*i,8*j-7:8*j);
        %  DCT  
        resim1_Blok_DCT = dct2(resim1_Blok);
        resim2_Blok_DCT = dct2(resim2_Blok);
        sigma1=svd(resim1_Blok_DCT);
        sigma2=svd(resim2_Blok_DCT);
        
        x1=sigma1(1)*sigma1(2)*sigma1(3)*sigma1(4)*sigma1(5);
        x2=sigma2(1)*sigma2(2)*sigma2(3)*sigma2(4)*sigma2(5);
        
         % Fusion
        if x1 > x2
            dctBlok = resim1_Blok_DCT;
            Map(i,j) =+1;	% CV
        else
            dctBlok = resim2_Blok_DCT;
            Map(i,j) = -1;    % CV 
        end
        
        % Inverse DCT 
        % DCT+SVD Method
        FusedDCT(8*i-7:8*i,8*j-7:8*j) = idct2(dctBlok);
        
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
            FusedDCT_CV(8*i-7:8*i,8*j-7:8*j) = resim1(8*i-7:8*i,8*j-7:8*j);
        else
            FusedDCT_CV(8*i-7:8*i,8*j-7:8*j) = resim2(8*i-7:8*i,8*j-7:8*j);
        end
        
    end
end
%**************************************************


% inverse shifting 
%**************************************************
resim1 = uint8(double(resim1)+128);
resim2 = uint8(double(resim2)+128);
%**************************************************

FusedDCT_CV = uint8(double(FusedDCT_CV)+128);
figure, imshow(FusedDCT_CV), title('"DCT+SVD+CV" ');
end
        

function popupmenu3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
