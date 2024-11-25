function level = graythresh(I) 
 
if ~isempty(I) 
    I = im2uint8(I(:)); 
    num_bins = 256; 
    counts = imhist(I,num_bins); 

    p = counts / sum(counts); 
    omega = cumsum(p); 
    mu = cumsum(p .* (1:num_bins)'); 
    mu_t = mu(end); 

    previous_state = warning('off', 'MATLAB:divideByZero'); 
    sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega)); 
    warning(previous_state); 

    maxval = max(sigma_b_squared);
    if isfinite(maxval) 
        idx = mean(find(sigma_b_squared == maxval)); 
     
        % Normalize the threshold to the range [0, 1]. 
        level = (idx - 1) / (num_bins - 1);
    else 
        level = 0.0;
    end 
else 
    level = 0.0; 
end 
