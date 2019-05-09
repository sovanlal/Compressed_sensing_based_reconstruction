function A = get_system_matrix(N,theta)
%GET_SYSTEM_MATRIX Set up tomography system matrix using parallel rays
% The N times N image is assumed to be centered around the origo. The side
% length of a pixel is 1, so the inter-pixel borders are the lines x = k,
% k=-N/2,-N/2+1,...,N/2+1,N/2 and y = k, k=-N/2,-N/2+1,...,N/2+1,N/2. theta
% is a vector with angles. The k'th entry of theta specifies the angle
% between the x-rays in the k'th view and the positive x-axis.
%
% A = get_system_matrix(N) returns the system matrix A, using the default
% value of theta (0-179 with 1 degree interval).
%
% A = get_system_matrix(N,theta) returns the system matrix A, using the
% specified values of theta.
%
% See also GET_TEST_IMAGE, ADD_NOISE, CIMMINO, CGLS,
% DISPLAY_RECONSTRUCTION, DEMO1, DEMO2, DEMO3
%
% Based on a function created by Jakob Heide Joergensen, Per Christian
% Hansen and Maria Saxild-Hansen, IMM, DTU, 2010.
% Modified by Sovan Mukherjee, April 2015

if nargin < 1
 errordlg('Not enough input arguments. You must at least specify the number of pixels N')
 return
elseif nargin < 2
 theta = 0:179;
end
if ~isscalar(N) || N < 1 || ~ceil(N)==floor(N)
 errordlg('N must be a positive integer')
 return
elseif ~isrow(theta) || size(theta,2)<1
 errordlg('Theta must be a row vector containing at least two elements')
 return
end
% Set the number of rays in a view
% p = ceil(sqrt(2)*N); p = p + ~mod(p,2);
p = N; 
% p = p + ~mod(p,2);

% Get the number of angles.
nA = length(theta);

% The starting values both the x and the y coordinates for theta=0 degrees,
% which corresponds to vertical rays. The spacing between each ray is 1.
x0 = (-(p-1)/2:(p-1)/2)';
y0 = zeros(p,1);
% The inter-pixel lines.
x = (-N/2:N/2)';
y = x;
% Preallocate the various indices for the sparse matrix to be created later
M=ceil(2*N*sqrt(2)*p*nA); % total amount of nonzeros in A matrix
I=zeros(1,M);
J=zeros(1,M);
X=zeros(1,M);
count = 1; % Define a counter for the indices
% Loop over the angles
for i = 1:nA

 % Rotate the points (x0,y0) by theta(i) around the origo to obtain the
 % starting points for the current angle
 x0theta = cosd(theta(i))*x0-sind(theta(i))*y0;
 y0theta = sind(theta(i))*x0+cosd(theta(i))*y0;

 % The direction vector for all the rays corresponding to the current
 % angle.
 a = -sind(theta(i));
 b = cosd(theta(i));

 % Loop over the single rays at the given angle
 for j = 1:p
 % Determine which pixels are hit by the current ray and the path
 % length inside each pixel. col holds the pixel index to the hit
 % pixels and d the corresponding path lengths
 [col,d] = trace_single_ray(N,x,x0theta(j),y,y0theta(j),a,b);

 % Create indices for creating the sparse matrix in the end (if col is not empty)
 if ~isempty(col)
 count_end = count+length(col)-1;
 I(count:count_end) = (i-1)*p+j;
 J(count:count_end) = col;
 X(count:count_end) = d;
 count = count_end+1;
 end
 end
end
% Create A as a sparse matrix
A = sparse(I(1:count-1),J(1:count-1),X(1:count-1),p*nA,N*N);
function [col,d] = trace_single_ray(N,x,x0theta,y,y0theta,a,b)
% Subfunction to compute which pixels are hit by a ray. The image is N by
% N, with inter-pixel lines at x and y. The ray goes through
% (x0theta,y0theta) and has direction (a,b).
% The idea is to determine the intersection points of the ray with each of
% the inter-pixel lines. From these points the path lengths and the index
% of the pixels can be determined.
% The ray is parametrized as [x(t);y(t)] = [a;b]*t + [x0;y0].
% The intersection points are found by determining the intersection time
% with all the x-inter-pixel lines and evaluating the parametrization for
% this time. And similarly for the y-inter-pixel lines.
% The distances between consecutive intersection points are then the path
% lengths inside the pixels.
% Use the parametrization of the ray to get the y-coordinates of
% intersection points with the x-inter-pixel lines x = x(k),
tx = (x - x0theta)/a;
yx = b*tx + y0theta;

% Use the parametrization of the ray to get the x-coordinates of
% intersection points with the y-inter-pixel lines y = y(k),
ty = (y - y0theta)/b;
xy = a*ty + x0theta;
% Collect the intersection times and coordinates.
t = [tx; ty];
xxy = [x; xy];
yxy = [yx; y];
% Sort the coordinates according to intersection time.
[t I] = sort(t);
xxy = xxy(I);
yxy = yxy(I);
% Skip the points outside the box, using find.
I = find(xxy >= -N/2 & xxy <= N/2 & yxy >= -N/2 & yxy <= N/2);
xxy = xxy(I);
yxy = yxy(I);
% If the ray goes from one pixel to another exactly through the corner,
% it intersects both an x- and a y-inter-pixel line at the same time, which
% makes this corner point occur twice in (xxy,yxy). These double points
% must be skipped.
I = find(abs(diff([xxy;inf])) <= 1e-10 & abs(diff([yxy;inf])) <= 1e-10);
xxy(I) = [];
yxy(I) = [];

% Calculate the length within each pixel and determine the number of
% pixels that are hit.
d = sqrt(diff(xxy).^2 + diff(yxy).^2);
numvals = numel(d);
col = [];
% Store the values inside the box.
if numvals > 0

 % Rays exactly along the inter-pixel lines must be assigned to either
 % the pixels one on side or the other. When using the floor function
 % this is automatically handled, except for at the rightmost and
 % topmost boundary of the image. These rays are considered outside, and
 % their path lengths are excluded by the if statement
 if ~((b == 0 && abs(y0theta - N/2) < 1e-15) || ...
 (a == 0 && abs(x0theta - N/2) < 1e-15))

 % Calculates the midpoints of the line within the cells
 xm = 0.5*(xxy(1:end-1)+xxy(2:end));
 ym = 0.5*(yxy(1:end-1)+yxy(2:end));

 % To determine the index of the hit pixels, translate the midpoint
 % coordinates by N/2, take the integer part and combine to single
 % index.
 col = floor(xm+N/2)*N + (N - floor(ym+N/2));
 end
end