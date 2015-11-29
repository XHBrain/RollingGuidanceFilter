function result = RollingGuidanceFilter(I, sigma_s, sigma_r, ...
     iteration,GaussianPrecision)
%function result = RollingGuidanceFilter(I, sigma_s, sigma_r, iteration,GaussianPrecision)
%is executed a Rolling Guidance Filter.
%site:http://www.cse.cuhk.edu.hk/leojia/projects/rollguidance
%input    I:the origic image; sigma_s: the Gaussian standard deviation; sigma_r: control the spatial and range weights respectively
          iteration: the number of iteration; GaussianPrecision: contorl the size of Block.
%output   result: generate a new image

if ~exist('iteration','var')
    iteration = 4;
end

if ~exist('sigma_s','var')
    sigma_s = 3;
end

if ~exist('sigma_r','var')
    sigma_r = 0.1;
end

if ~exist('GaussianPrecision','var')
    GaussianPrecision = 0.08;
end
 
[a b c] = size(I);

GaussianWindow = WindowBlock(sigma_s, GaussianPrecision);

N = ( size(GaussianWindow,1) -1)/2;

J_plus = zeros(size(I));
I = ExpandBorder(I, N);

J = ones(size(I));

for l = 1 : iteration
    for k = 1 : c
        for j = 1+N : b+N
            for i = 1+N : a+N
                H = GaussianWindow .* exp( -(J(i-N:i+N,j-N:j+N,k) - J(i,j,k)).^2/2/sigma_r^2 );
                K_p = sum(sum(H));
                J_plus(i-N,j-N,k) = 1/K_p*sum(sum(H .* I(i-N:i+N,j-N:j+N,k)));                
            end
        end
    end
    %figure;   image(J_plus);
    J = ExpandBorder(J_plus, N);
    
end
    result = J_plus;
    
end

%==================================================================

function GaussianWindow = WindowBlock(sigma_s, GaussianPrecision)
%function GaussianWindow = WindowBlock(sigma_s, GaussianPrecision)
%product a Block of 2D Gaussian.
%input    sigma_s: the Gaussian standard deviation; GaussianPrecision: contorl the size of Block.
%output   GaussianWindow: a square matrix of 2D Gaussian.

%right below
pq = bsxfun(@plus, ([0:sigma_s*3].^2)', [0:sigma_s*3].^2); 

% gaussian distribution
pqrb = exp(-pq/2/sigma_s^2); 

% element that is less than GaussianPrecision ar equal zero
pqrb = pqrb .* (pqrb>GaussianPrecision); 

% remove all zero column
pqrb(:, pqrb(1,:)==0) = []; 

% remove all zero row
pqrb(pqrb(:,1)==0, :) = []; 

%left below
pqlb = fliplr(pqrb); 

%right upper
pqru = flipud(pqrb); 

%left upper
pqlu = fliplr(pqru); 

GaussianWindow = [pqlu(:, 1:end-1)     pqru;
                  pqlb(2:end, 1:end-1) pqrb(2:end, :)];
  
end

%==================================================================

function imageExpand = ExpandBorder(image, N)
%function imageExpand = ExpandBorder(image, N)
%is On the edge of the pixels extend outward N pixels.
%input    image: the origic image; N: the number of expansion need.
%output   imageExpand: generate a new image

imageExpand = [repmat(image(1,:,:), [N,1,1]) ; image ; repmat(image(end,:,:), [N,1,1])];
imageExpand = [repmat(imageExpand(:,1,:), [1,N,1])  imageExpand  repmat(imageExpand(:,end,:), [1,N,1])];

end
