%% ---- Path to Siemens raw (.dat)
fn = '/Users/ekintaskin/Desktop/pulseq/pypulseq/fat_sep/data/meas_MID00086_FID319290_256_v2_dcm.dat';

%% ---- Load with mapVBVD
twix = mapVBVD(fn);
if iscell(twix), twix = twix{end}; end     % keep imaging dataset
twix.image.flagRemoveOS = false;           % keep 256 RO samples

%% ---- Pull k-space and make dims explicit
k    = twix.image();                       % complex raw data
k    = squeeze(k);                         % drop length-1 dims so sqzDims matches
dims = twix.image.sqzDims;                 % names of each dimension in 'k'
fprintf('Raw (squeezed) k-space size: [%s]\n', num2str(size(k)));

% Find named dimensions in current array
ixCol = find(strcmpi(dims,'Col'));  assert(~isempty(ixCol),'No Col dim.');
ixLin = find(strcmpi(dims,'Lin'));  assert(~isempty(ixLin),'No Lin dim.');
ixCha = find(strcmpi(dims,'Cha'));  % may be empty if already combined on scanner
ixEco = find(strcmpi(dims,'Eco'));  assert(~isempty(ixEco),'No Eco dim.');

% Reorder to [Col, Lin, (Cha), Eco, ...rest]
perm = [ixCol ixLin];
if ~isempty(ixCha), perm(end+1) = ixCha; end
perm(end+1) = ixEco;
perm = [perm setdiff(1:ndims(k), perm, 'stable')];
k = permute(k, perm);

% Sizes
sz = size(k);
Nx = sz(1); Ny = sz(2);
hasCha = ~isempty(ixCha) && (numel(sz)>=3) && (sz(3) > 1);
Nc = hasCha*sz(3) + (~hasCha)*1;
Ne = sz(2 + 1 + hasCha);   % echo dim index after the permute above

fprintf('k-space size after reorder: [%s]\n', num2str(size(k)));
assert(Nx==256 && Ny==256, 'kx/ky are not 256x256; set flagRemoveOS=false to keep 256.');

%% ---- 2D centered FFT (Col & Lin only) -> image space [Nx Ny (Nc) Ne ...]
k = fftshift(ifft(ifftshift(k,1),[],1),1);
k = fftshift(ifft(ifftshift(k,2),[],2),2);

%% ---- Complex coil combination (create img_echo [Nx Ny Ne], COMPLEX!)
if hasCha && Nc > 1
    % --- Adaptive/Walsh combine using echo 1 to compute weights
    ref_all = k(:,:,1:Nc,1);                 % [Nx Ny Nc]
    w = zeros(size(ref_all), 'like', ref_all);
    win = 5;                                  % small local window
    for y = 1:Ny
        y1 = max(1,y-win); y2 = min(Ny,y+win);
        for x = 1:Nx
            x1 = max(1,x-win); x2 = min(Nx,x+win);
            P = reshape(ref_all(x1:x2,y1:y2,:), [], Nc);   % [samples x Nc]
            C = (P' * conj(P)) + 1e-6*eye(Nc);             % coil covariance
            [V,~] = eigs(C,1,'largestreal');               % principal eigenvector
            w(x,y,:) = V;
        end
    end
    img_echo = zeros(Nx,Ny,Ne,'like',ref_all(:,:,1));
    for e = 1:Ne
        S = k(:,:,1:Nc,e);
        img_echo(:,:,e) = sum(conj(w).*S, 3);              % complex combine
    end

    % ---- (Optional) sanity check with Sum-of-Squares magnitude:
    % img_sos = sqrt(sum(abs(k(:,:,1:Nc,:)).^2, 3)); % [Nx Ny Ne] magnitude only
else
    % Single coil or already combined by scanner
    img_echo = squeeze(k(:,:,1,1:Ne));       % -> [Nx Ny Ne] complex
end

% Ensure echoes are along 3rd dim
img_echo = squeeze(img_echo);
assert(ndims(img_echo)==3 && size(img_echo,3)>=2, ...
    'Echo dimension not found or <2; check dims and Ne.');

%% ---- Choose echoes (set to your in/opposed pair)
ip = 3; op = 4;   % <- adjust if needed (or use the auto-pick block below)

% Optional: auto-pick using header TEs and B0
%{
B0 = 3.0; gamma = 42.577478518e6; delta_ppm = 3.5e-6;
TEs = [];
try
    if isfield(twix.hdr,'Meas') && isfield(twix.hdr.Meas,'alTE')
        TEs = double([twix.hdr.Meas.alTE{:}]) * 1e-6;  % seconds
        disp('TEs (s):'); disp(TEs(:).');
    end
catch, end
if ~isempty(TEs) && numel(TEs) >= 2
    df  = gamma*B0*delta_ppm;
    phi = mod(2*pi*df*TEs, 2*pi);
    [~,ip] = min(abs(angdiff(phi, 0)));
    [~,op] = min(abs(angdiff(phi, pi)));
    fprintf('Auto-picked echoes: IP=%d, OP=%d\n', ip, op);
end
%}

S1 = img_echo(:,:,ip);
S2 = img_echo(:,:,op);

%% ---- Phase alignment between echoes (global + optional linear plane)
mask = abs(S1) > 0.2*max(abs(S1(:)));
phi0 = angle(sum(S1(mask) .* conj(S2(mask))));  % global phase
S2a  = S2 .* exp(1i*phi0);

% (Optional) linear phase plane (helps with bipolar readouts)
%{
dphi = angle(S1 .* conj(S2a)); dphi(~mask)=0;
[Nx,Ny] = size(S1);
[xg, yg] = ndgrid((0:Nx-1)/Nx - 0.5, (0:Ny-1)/Ny - 0.5);
A = [xg(mask), yg(mask), ones(nnz(mask),1)];
p = A \ wrapToPi(dphi(mask));                 % [ax; by; c]
phi_map = p(1)*xg + p(2)*yg + p(3);
S2a = S2a .* exp(-1i*phi_map);
%}

%% ---- 2-echo Dixon in IMAGE space
W  = 0.5*(S1 + S2a);
F  = 0.5*(S1 - S2a);
FF = abs(F) ./ (abs(F) + abs(W) + eps);
%% ---- Show
figure('Name','2-Echo Dixon','NumberTitle','off'); clf;
subplot(2,3,1), imagesc(abs(S1)), axis image off, colormap gray, title(sprintf('Echo %d |S|',ip))
subplot(2,3,2), imagesc(abs(S2)), axis image off, title(sprintf('Echo %d |S|',op))
subplot(2,3,3), imagesc(angle(S1.*conj(S2))), axis image off, title('Phase diff (rad)')
subplot(2,3,4), imagesc(abs(W)), axis image off, title('|Water|')
subplot(2,3,5), imagesc(abs(F)), axis image off, title('|Fat|')
subplot(2,3,6), imagesc(FF,[0 1]), axis image off, title('Fat Fraction')
drawnow; shg;
