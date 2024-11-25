function varargout = colorspace(Conversion,varargin)

if nargin < 2, error('Not enough input arguments.'); end
[SrcSpace,DestSpace] = parse(Conversion);

if nargin == 2
   Image = varargin{1};
elseif nargin >= 3
   Image = cat(3,varargin{:});
else
   error('Invalid number of input arguments.');
end

FlipDims = (size(Image,3) == 1);

if FlipDims, Image = permute(Image,[1,3,2]); end
if ~isa(Image,'double'), Image = double(Image)/255; end
if size(Image,3) ~= 3, error('Invalid input size.'); end

SrcT = gettransform(SrcSpace);
DestT = gettransform(DestSpace);

if ~ischar(SrcT) && ~ischar(DestT)
    T = [DestT(:,1:3)*SrcT(:,1:3),DestT(:,1:3)*SrcT(:,4)+DestT(:,4)];      
    Temp = zeros(size(Image));
    Temp(:,:,1) = T(1)*Image(:,:,1) + T(4)*Image(:,:,2) + T(7)*Image(:,:,3) + T(10);
    Temp(:,:,2) = T(2)*Image(:,:,1) + T(5)*Image(:,:,2) + T(8)*Image(:,:,3) + T(11);
    Temp(:,:,3) = T(3)*Image(:,:,1) + T(6)*Image(:,:,2) + T(9)*Image(:,:,3) + T(12);
    Image = Temp;
elseif ~ischar(DestT)
    Image = rgb(Image,SrcSpace);
    Temp = zeros(size(Image));
    Temp(:,:,1) = DestT(1)*Image(:,:,1) + DestT(4)*Image(:,:,2) + DestT(7)*Image(:,:,3) + DestT(10);
    Temp(:,:,2) = DestT(2)*Image(:,:,1) + DestT(5)*Image(:,:,2) + DestT(8)*Image(:,:,3) + DestT(11);
    Temp(:,:,3) = DestT(3)*Image(:,:,1) + DestT(6)*Image(:,:,2) + DestT(9)*Image(:,:,3) + DestT(12);
    Image = Temp;
else
	Image = feval(DestT,Image,SrcSpace);
end

%%% Output format %%%
if nargout > 1
    varargout = {Image(:,:,1),Image(:,:,2),Image(:,:,3)};
else
    if FlipDims, Image = permute(Image,[1,3,2]); end
    varargout = {Image};
end

return;


function [SrcSpace,DestSpace] = parse(Str)

if ischar(Str)
    Str = lower(strrep(strrep(Str,'-',''),'=',''));
    k = find(Str == '>');

    if length(k) == 1        
        SrcSpace = Str(1:k-1);
        DestSpace = Str(k+1:end);
    else
        k = find(Str == '<');
        
        if length(k) == 1     
            DestSpace = Str(1:k-1);
            SrcSpace = Str(k+1:end);
        else
            error(['Invalid conversion, ''',Str,'''.']);
        end   
    end

    SrcSpace = alias(SrcSpace);
    DestSpace = alias(DestSpace);
else
    SrcSpace = 1;             
    DestSpace = Conversion;
    if any(size(Conversion) ~= 3), error('Transformation matrix must be 3x3.'); end
end
return;


function Space = alias(Space)

Space = strrep(strrep(Space,'cie',''),' ','');

if isempty(Space)
    Space = 'rgb';
end

switch Space
case {'hsi'}
    Space = 'hsi';
case {'rgb'}
    return;
end
return;


function T = gettransform(Space)

switch Space
case {'rgb','hsi'}
    T = Space;
end
return;


function Image = rgb(Image,SrcSpace)

switch SrcSpace
case 'rgb'
    return;
case 'hsv'
    Image = huetorgb((1 - Image(:,:,2)).*Image(:,:,3),Image(:,:,3),Image(:,:,1));
case 'hsi'
    L = Image(:,:,3);
    Delta = Image(:,:,2).*min(L,1-L);
    Image = huetorgb(L-Delta,L+Delta,Image(:,:,1));
end

Image = min(max(Image,0),1);
return;


function Image = hsi(Image,SrcSpace)

switch SrcSpace
case 'hsv'
    MaxVal = Image(:,:,3);
    MinVal = (1 - Image(:,:,2)).*MaxVal;
    L = 0.5*(MaxVal + MinVal);
    temp = min(L,1-L);
    Image(:,:,2) = 0.5*(MaxVal - MinVal)./(temp + (temp == 0));
    Image(:,:,3) = L;
otherwise
    Image = rgb(Image,SrcSpace);  
    MinVal = min(Image,[],3);
    MaxVal = max(Image,[],3);
    L = 0.5*(MaxVal + MinVal);
    temp = min(L,1-L);
    S = 0.5*(MaxVal - MinVal)./(temp + (temp == 0));
    Image(:,:,1) = rgbtohue(Image);
    Image(:,:,2) = S;
    Image(:,:,3) = L;
end
return;


function Image = huetorgb(m0,m2,H)
% Convert hsi hue to RGB
N = size(H);
H = min(max(H(:),0),360)/60;
m0 = m0(:);
m2 = m2(:);
F = H - round(H/2)*2;
M = [m0, m0 + (m2-m0).*abs(F), m2];
Num = length(m0);
j = [2 1 0;1 2 0;0 2 1;0 1 2;1 0 2;2 0 1;2 1 0]*Num;
k = floor(H) + 1;
Image = reshape([M(j(k,1)+(1:Num).'),M(j(k,2)+(1:Num).'),M(j(k,3)+(1:Num).')],[N,3]);
return;


function H = rgbtohue(Image)
% Convert RGB to HSI hue
[M,i] = sort(Image,3);
i = i(:,:,3);
Delta = M(:,:,3) - M(:,:,1);
Delta = Delta + (Delta == 0);
R = Image(:,:,1);
G = Image(:,:,2);
B = Image(:,:,3);
H = zeros(size(R));
k = (i == 1);
H(k) = (G(k) - B(k))./Delta(k);
k = (i == 2);
H(k) = 2 + (B(k) - R(k))./Delta(k);
k = (i == 3);
H(k) = 4 + (R(k) - G(k))./Delta(k);
H = 60*H + 360*(H < 0);
H(Delta == 0) = nan;
return;