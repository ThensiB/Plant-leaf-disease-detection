function [out] = Texture_Features(glcmin,pairs)

if ((nargin > 2) || (nargin == 0))
   error('Too many or too few input arguments. Enter GLCM and pairs.');
elseif ( (nargin == 2) ) 
    if ((size(glcmin,1) <= 1) || (size(glcmin,2) <= 1))
       error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcmin,1) ~= size(glcmin,2) )
        error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end    
elseif (nargin == 1) % only GLCM is entered
    pairs = 0; % default is numbers and input 1 for percentage
    if ((size(glcmin,1) <= 1) || (size(glcmin,2) <= 1))
       error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcmin,1) ~= size(glcmin,2) )
       error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end    
end


format long e
if (pairs == 1)
    newn = 1;
    for nglcm = 1:2:size(glcmin,3)
        glcm(:,:,newn)  = glcmin(:,:,nglcm) + glcmin(:,:,nglcm+1);
        newn = newn + 1;
    end
elseif (pairs == 0)
    glcm = glcmin;
end

size_glcm_1 = size(glcm,1);
size_glcm_2 = size(glcm,2);
size_glcm_3 = size(glcm,3);

% display 
out.E = zeros(1,size_glcm_3); % Energy
out.cov = zeros(1,size_glcm_3); % Correlation
out.se = zeros(1,size_glcm_3); % Sum entropy
out.de = zeros(1,size_glcm_3); % Difference entropy
out.e = zeros(1,size_glcm_3); % Entropy
out.inf2h = zeros(1,size_glcm_3); % Informaiton measure of correlation
out.id = zeros(1,size_glcm_3); % Contrast
out.c = zeros(1,size_glcm_3); % Correlation

glcm_sum  = zeros(size_glcm_3,1);
glcm_mean = zeros(size_glcm_3,1);
glcm_var  = zeros(size_glcm_3,1);

u_x = zeros(size_glcm_3,1);
u_y = zeros(size_glcm_3,1);
s_x = zeros(size_glcm_3,1);
s_y = zeros(size_glcm_3,1);

p_x = zeros(size_glcm_1,size_glcm_3);   
p_y = zeros(size_glcm_2,size_glcm_3);
p_xplusy = zeros((size_glcm_1*2 - 1),size_glcm_3); 
p_xminusy = zeros((size_glcm_1),size_glcm_3); 
% hxy hxy1 hxy2 hx hy
hxy  = zeros(size_glcm_3,1);
hxy1 = zeros(size_glcm_3,1);
hx   = zeros(size_glcm_3,1);
hy   = zeros(size_glcm_3,1);
hxy2 = zeros(size_glcm_3,1);

for k = 1:size_glcm_3 
    glcm_sum(k) = sum(sum(glcm(:,:,k)));
    glcm(:,:,k) = glcm(:,:,k)./glcm_sum(k); % Normalize each glcm
    glcm_mean(k) = mean2(glcm(:,:,k)); % compute mean after norm
    glcm_var(k) = (std2(glcm(:,:,k)))^2;
    
    for i = 1:size_glcm_1
        for j = 1:size_glcm_2
            out.id(k) = out.id(k) + (glcm(i,j,k)/(1+(abs(i-j)).^2));
            out.E(k) = out.E(k) + (glcm(i,j,k).^2);
            out.e(k) = out.e(k) - (glcm(i,j,k)*log(glcm(i,j,k) + eps));
            u_x(k) = u_x(k) + (i)*glcm(i,j,k);
            u_y(k) = u_y(k) + (j)*glcm(i,j,k);
        end
    end
end

for k = 1:size_glcm_3
    for i = 1:size_glcm_1
        for j = 1:size_glcm_2
            p_x(i,k) = p_x(i,k) + glcm(i,j,k); 
            p_y(i,k) = p_y(i,k) + glcm(j,i,k);
            if (ismember((i + j),[2:2*size_glcm_1])) 
                p_xplusy((i+j)-1,k) = p_xplusy((i+j)-1,k) + glcm(i,j,k);
            end
            if (ismember(abs(i-j),[0:(size_glcm_1-1)])) 
                p_xminusy((abs(i-j))+1,k) = p_xminusy((abs(i-j))+1,k) +...
                    glcm(i,j,k);
            end
        end
    end   
end

% computing sum and difference entropy
for k = 1:(size_glcm_3)
    for i = 1:(2*(size_glcm_1)-1)
        out.se(k) = out.se(k) - (p_xplusy(i,k)*log(p_xplusy(i,k) + eps));
    end
end

for k = 1:size_glcm_3
    for i = 0:(size_glcm_1-1)
        out.de(k) = out.de(k) - (p_xminusy(i+1,k)*log(p_xminusy(i+1,k) + eps));
    end
end

% compute information measure of correlation
for k = 1:size_glcm_3
    hxy(k) = out.e(k);
    for i = 1:size_glcm_1
        for j = 1:size_glcm_2
            hxy2(k) = hxy2(k) - (p_x(i,k)*p_y(j,k)*log(p_x(i,k)*p_y(j,k) + eps));
        end
    end
    out.inf2h(k) = ( 1 - exp( -2*( hxy2(k) - hxy(k) ) ) )^0.5;
end

corm = zeros(size_glcm_3,1);
corp = zeros(size_glcm_3,1);
for k = 1:size_glcm_3
    for i = 1:size_glcm_1
        for j = 1:size_glcm_2
            s_x(k)  = s_x(k)  + (((i) - u_x(k))^2)*glcm(i,j,k);
            s_y(k)  = s_y(k)  + (((j) - u_y(k))^2)*glcm(i,j,k);
            corp(k) = corp(k) + ((i)*(j)*glcm(i,j,k));
            corm(k) = corm(k) + (((i) - u_x(k))*((j) - u_y(k))*glcm(i,j,k));
        end
    end
    s_x(k) = s_x(k) ^ 0.5;
    s_y(k) = s_y(k) ^ 0.5;
    out.c(k) = (corp(k) - u_x(k)*u_y(k))/(s_x(k)*s_y(k));
    out.cov(k) = corm(k);
end