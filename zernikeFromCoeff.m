function zernikeFromCoeff(coef, c, k, type, titleSag, titleSurf)
% 直接使用CodeV或Zemax优化得到的Zernike系数仿真镜面
% coef: Code V
% 求出的参数；c:曲面顶点处的曲率；k:圆锥曲面系数；titleSag：平面图的标题；titleSurf：三维图的标题；type：一般zernike或XY
% zernike
% example:
% coef = [3.608, -0.00519, 0.89, -0.00050];
% c = 0.001994;
% k = -5.6933;
% titleSag = 'sag';
% titleSurf = 'surf';
% References: 
% [1] 刘军. 自由曲面在成像光学系统中的研究[D]. 中国科学院研究生院(长春光学精密机械与物理研究所), 2016.
% [2] 王超. 自由曲面表征函数及其应用研究[D]. 中国科学院研究生院(长春光学精密机械与物理研究所), 2014.
% zernikeFromCoeff(coef, c, k, 0, titleSag, titleSurf);
    coef = [1.0,coef];
    len = length(coef);
    %nm = [0,0; 1,-1; 1,1; 2,-2; 2,0; 2,2; 3,-3; 3,-1; 3,1; 3,3; 4,-4; 4,-2; 4,0; 4,2; 4,4];
    nm = [0,0; 1,1; 1,-1; 2,2; 2,0; 2,-2; 3,3; 3,1; 3,-1; 3,-3; 4,4; 4,2; 4,0; 4,-2; 4,-4];
    xy = cell(15,1);
    xy{1,1} = [1,1];
    xy{2,1} = [2,1];
    xy{3,1} = [3,1];
    xy{4,1} = [1,4;4,2;5,4];
    xy{5,1} = [6,2];
    xy{6,1} = [1,4;4,-2;5,4];
    xy{7,1} = [2,2;7,4;8,4];
    xy{8,1} = [10,4;3,6;9,12];
    xy{9,1} = [2,6;7,-4;8,12];
    xy{10,1} = [3,2;10,-4;9,4];
    xy{11,1} = [1,8;11,8;12,8;13,16;4,2.667;5,5.333];
    xy{12,1} = [14,16;15,8;6,5.333];
    xy{13,1} = [1,24;11,-8;13,48;5,16];
    xy{14,1} = [14,16;15,-8;6,5.333];
    xy{15,1} = [1,8;11,8;12,-8;13,16;4,-2.667;5,5.333];

    if len>length(nm)
        error('the coef is too long!');
    end

    if type==0
        x = -1:0.01:1;
        [X,Y] = meshgrid(x,x);
        [theta,r] = cart2pol(X,Y);
        idx = r<=1;
        % z = cell(len,1);
        z = calcZer(X, idx, len, r, theta, nm);
        for j=1:len
            %z{j} = nan(size(X));
            z{j}(idx) = z{j}(idx)*coef(j);
        end
    else
        x = -1:0.01:1;
        [X,Y] = meshgrid(x,x);
        [theta,r] = cart2pol(X,Y);
        idx = r<=1;
        % z = cell(len,1);
        zz = calcZer(X, idx, 15, r, theta, nm);
        z = cell(len,1);
        for q=1:len
            z{q} = nan(size(X));
            for p = 1:size(xy{q},1)
                if p==1
                    z{q} = zz{xy{q}(p,1)}/xy{q}(p,2);
                else
                    z{q} = z{q}+zz{xy{q}(p,1)}/xy{q}(p,2);
                end
            end
            z{q} = z{q}*coef(q);
        end
        
        % for i=1:len
        %     z{i} = nan(size(X));
        %     z{i}(idx) = zernfun(nm(i,1),nm(i,2),r(idx),theta(idx))*coef(i);
        % end
    end
    rr = r.^2+theta.^2;
    res = c*rr./(1+(1-(1+k)*c^2*rr).^0.5);
    for i=1:len
        res = res+z{i,1};
    end

    figure,
    pcolor(x,x,res), shading interp;
    axis square, colorbar, colormap jet;
    title(titleSag);

    figure,
    surf(X, Y, res), shading interp, colormap jet;
    axis square, colorbar;
    title(titleSurf);
end

function z = calcZer(X, idx, len, r, theta, nm)
z = cell(len,1);
    for j=1:len
        z{j} = nan(size(X));
        z{j}(idx) = zernfun(nm(j,1),nm(j,2),r(idx),theta(idx));
    end
end


function z = zernfun(n,m,r,theta,nflag)
%ZERNFUN Zernike functions of order N and frequency M on the unit circle.
%   Z = ZERNFUN(N,M,R,THETA) returns the Zernike functions of order N
%   and angular frequency M, evaluated at positions (R,THETA) on the
%   unit circle.  N is a vector of positive integers (including 0), and
%   M is a vector with the same number of elements as N.  Each element
%   k of M must be a positive integer, with possible values M(k) = -N(k)
%   to +N(k) in steps of 2.  R is a vector of numbers between 0 and 1,
%   and THETA is a vector of angles.  R and THETA must have the same
%   length.  The output Z is a matrix with one column for every (N,M)
%   pair, and one row for every (R,THETA) pair.
%
%   Z = ZERNFUN(N,M,R,THETA,'norm') returns the normalized Zernike
%   functions.  The normalization factor sqrt((2-delta(m,0))*(n+1)/pi),
%   with delta(m,0) the Kronecker delta, is chosen so that the integral
%   of (r * [Znm(r,theta)]^2) over the unit circle (from r=0 to r=1,
%   and theta=0 to theta=2*pi) is unity.  For the non-normalized
%   polynomials, max(Znm(r=1,theta))=1 for all [n,m].
%
%   The Zernike functions are an orthogonal basis on the unit circle.
%   They are used in disciplines such as astronomy, optics, and
%   optometry to describe functions on a circular domain.
%
%   The following table lists the first 15 Zernike functions.
%
%       n    m    Zernike function             Normalization
%       ----------------------------------------------------
%       0    0    1                              1/sqrt(pi)
%       1    1    r * cos(theta)                 2/sqrt(pi)
%       1   -1    r * sin(theta)                 2/sqrt(pi)
%       2    2    r^2 * cos(2*theta)             sqrt(6/pi)
%       2    0    (2*r^2 - 1)                    sqrt(3/pi)
%       2   -2    r^2 * sin(2*theta)             sqrt(6/pi)
%       3    3    r^3 * cos(3*theta)             sqrt(8/pi)
%       3    1    (3*r^3 - 2*r) * cos(theta)     sqrt(8/pi)
%       3   -1    (3*r^3 - 2*r) * sin(theta)     sqrt(8/pi)
%       3   -3    r^3 * sin(3*theta)             sqrt(8/pi)
%       4    4    r^4 * cos(4*theta)             sqrt(10/pi)
%       4    2    (4*r^4 - 3*r^2) * cos(2*theta) sqrt(10/pi)
%       4    0    6*r^4 - 6*r^2 + 1              sqrt(5/pi)
%       4   -2    (4*r^4 - 3*r^2) * sin(2*theta) sqrt(10/pi)
%       4   -4    r^4 * sin(4*theta)             sqrt(10/pi)
%       ----------------------------------------------------
%
%   Example 1:
%
%       % Display the Zernike function Z(n=5,m=1)
%       x = -1:0.01:1;
%       [X,Y] = meshgrid(x,x);
%       [theta,r] = cart2pol(X,Y);
%       idx = r<=1;
%       z = nan(size(X));
%       z(idx) = zernfun(5,1,r(idx),theta(idx));
%       figure
%       pcolor(x,x,z), shading interp
%       axis square, colorbar
%       title('Zernike function Z_5^1(r,\theta)')
%
%   Example 2:
%
%       % Display the first 10 Zernike functions
%       x = -1:0.01:1;
%       [X,Y] = meshgrid(x,x);
%       [theta,r] = cart2pol(X,Y);
%       idx = r<=1;
%       z = nan(size(X));
%       n = [0  1  1  2  2  2  3  3  3  3];
%       m = [0 -1  1 -2  0  2 -3 -1  1  3];
%       Nplot = [4 10 12 16 18 20 22 24 26 28];
%       y = zernfun(n,m,r(idx),theta(idx));
%       figure('Units','normalized')
%       for k = 1:10
%           z(idx) = y(:,k);
%           subplot(4,7,Nplot(k))
%           pcolor(x,x,z), shading interp
%           set(gca,'XTick',[],'YTick',[])
%           axis square
%           title(['Z_{' num2str(n(k)) '}^{' num2str(m(k)) '}'])
%       end
%
%   See also ZERNPOL, ZERNFUN2.

%   Paul Fricker 2/28/2012

% Check and prepare the inputs:
% -----------------------------
if ( ~any(size(n)==1) ) || ( ~any(size(m)==1) )
    error('zernfun:NMvectors','N and M must be vectors.')
end

if length(n)~=length(m)
    error('zernfun:NMlength','N and M must be the same length.')
end

n = n(:);
m = m(:);
if any(mod(n-m,2))
    error('zernfun:NMmultiplesof2', ...
          'All N and M must differ by multiples of 2 (including 0).')
end

if any(m>n)
    error('zernfun:MlessthanN', ...
          'Each M must be less than or equal to its corresponding N.')
end

if any( r>1 | r<0 )
    error('zernfun:Rlessthan1','All R must be between 0 and 1.')
end

if ( ~any(size(r)==1) ) || ( ~any(size(theta)==1) )
    error('zernfun:RTHvector','R and THETA must be vectors.')
end

r = r(:);
theta = theta(:);
length_r = length(r);
if length_r~=length(theta)
    error('zernfun:RTHlength', ...
          'The number of R- and THETA-values must be equal.')
end

% Check normalization:
% --------------------
if nargin==5 && ischar(nflag)
    isnorm = strcmpi(nflag,'norm');
    if ~isnorm
        error('zernfun:normalization','Unrecognized normalization flag.')
    end
else
    isnorm = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Zernike Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the required powers of r:
% -----------------------------------
m_abs = abs(m);
rpowers = [];
for j = 1:length(n)
    rpowers = [rpowers m_abs(j):2:n(j)];
end
rpowers = unique(rpowers);

% Pre-compute the values of r raised to the required powers,
% and compile them in a matrix:
% -----------------------------
if rpowers(1)==0
    rpowern = arrayfun(@(p)r.^p,rpowers(2:end),'UniformOutput',false);
    rpowern = cat(2,rpowern{:});
    rpowern = [ones(length_r,1) rpowern];
else
    rpowern = arrayfun(@(p)r.^p,rpowers,'UniformOutput',false);
    rpowern = cat(2,rpowern{:});
end

% Compute the values of the polynomials:
% --------------------------------------
z = zeros(length_r,length(n));
for j = 1:length(n)
    s = 0:(n(j)-m_abs(j))/2;
    pows = n(j):-2:m_abs(j);
    for k = length(s):-1:1
        p = (1-2*mod(s(k),2))* ...
                   prod(2:(n(j)-s(k)))/              ...
                   prod(2:s(k))/                     ...
                   prod(2:((n(j)-m_abs(j))/2-s(k)))/ ...
                   prod(2:((n(j)+m_abs(j))/2-s(k)));
        idx = (pows(k)==rpowers);
        z(:,j) = z(:,j) + p*rpowern(:,idx);
    end
    
    if isnorm
        z(:,j) = z(:,j)*sqrt((1+(m(j)~=0))*(n(j)+1)/pi);
    end
end
% END: Compute the Zernike Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the Zernike functions:
% ------------------------------
idx_pos = m>0;
idx_neg = m<0;

if any(idx_pos)
    z(:,idx_pos) = z(:,idx_pos).*cos(theta*m_abs(idx_pos)');
end
if any(idx_neg)
    z(:,idx_neg) = z(:,idx_neg).*sin(theta*m_abs(idx_neg)');
end
end
% EOF zernfun