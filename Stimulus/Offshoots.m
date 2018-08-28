%% Bits of code, no longer needed but might be useful
%% AG.Mitchell 12.07.18

%% Creating a noise mask
% Noise mask variables
mask = struct;
mask.size.pix = 600;
mask.size.text = mask.size.pix/2;
mask.contrast = 0.8; 
mask.freq.sec = 1; %freq cycles per second
mask.freq.pix = 0.01; %freq cycles per pixel
mask.freq.pixper = ceil(1/mask.freq.pix); %pixels per cycle rounded to nearest pixel
mask.freq.rad = mask.freq.pix*2*pi; %frequency in radians
mask.size.vis = 2*mask.size.text + 1; %visible size of the grating

% Define the grating
mask.m = meshgrid(-mask.size.vis:mask.size.vis + mask.freq.pixper, 1);
mask.grating = grey * cos(mask.freq.rad * mask.m) + grey;
mask.layer = ones(1, numel(mask.m), 2) * grey; %make 2 layer mask filled with the background colour
% Place grating in 'alpha' channel of mask
mask.layer(:,:,2) = mask.grating .* mask.contrast;
% Rect
mask.rect.coord = [0 0 mask.size.vis mask.size.vis];
mask.rect.d = CenterRect(mask.rect.coord, rect);

% Draw mask
mask.tex = Screen('MakeTexture', window, mask.layer); %making texture
% Make black and white noise mask
mask.noise = rand(round(mask.size.vis/2)) .* white;
mask.noisetex = Screen('MakeTexture', window, mask.noise);
% Drawing mask
Screen('DrawTexture', window, mask.noisetex, [], mask.rect.d, []);
Screen('Flip', window);
time.mask.start(i) = GetSecs;
WaitSecs(stim.time.mask);
Screen('Flip', window);
time.mask.end(i) = GetSecs;
WaitSecs(0.2); %wait a further 200ms
