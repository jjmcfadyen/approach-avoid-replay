function get_corticalMesh()

maxvertices = 8200;

directory = 'D:\2020_RiskyReplay\data\mri';
addpath D:\Toolboxes\fieldtrip-20191119\external\iso2mesh

for h = 1:2 % for each hemisphere separately...
    
    if h==1
        htag = 'l';
    elseif h==2
        htag = 'r';
    end
    
    % Load gifti
    G = gifti(fullfile(directory,[htag 'h.orig.gii']));
    
    % Downsample
    dsfactor = maxvertices/size(G.vertices,1);
    [v,f] = meshresample(G.vertices,G.faces,dsfactor);
    
    orig = G;
    G.faces = int32(f);
    G.vertices = single(v);
    
    % Save
    save(G,fullfile(directory,[htag 'h.orig.downsampled.gii']));

    % Plot
    R=-1;
    if h==2
        R=1;
    end
    figure
    
    subplot(1,2,1)
    patch('vertices',orig.vertices,'faces',orig.faces,'facealpha',0.1,'edgealpha',0.1);
    axis equal
    view(70*R,10)
    title('Original')
    
    subplot(1,2,2)
    patch('vertices',G.vertices,'faces',G.faces,'facealpha',0.1,'edgealpha',0.3);
    axis equal
    view(70*R,10)
    title('Downsampled')
    
    if h==1
        sgtitle('Left Hemisphere')
    elseif h==2
        sgtitle('Right Hemisphere')
    end
    
end

% Combine left and right hemispheres
L = gifti(fullfile(directory,'lh.orig.downsampled.gii'));
R = gifti(fullfile(directory,'rh.orig.downsampled.gii'));

[v,f] = mergemesh(L.vertices,L.faces,R.vertices,R.faces);

[v,f] = meshcheckrepair(v,f);

M = L;
M.faces = f(:,1:3);
M.vertices = v;

% Save
save(M,fullfile(directory,'merged.orig.downsampled.gii'));

M.vertices = spm_eeg_inv_transform_points(M.mat, M.vertices);
save(M,fullfile(directory,'merged.orig.downsampled.warped.gii'));

% Plot
figure
patch('vertices',M.vertices,'faces',M.faces,'facealpha',0.1,'edgealpha',0.3);
axis equal
view(-70,10)
title('Merged')

%% Compare to SPM mesh

% Load giftis
M = gifti(fullfile(directory,'merged.orig.downsampled.gii'));
S = gifti('D:\Toolboxes\spm12\canonical\orig\cortex_8196.surf.gii');

% Plot
figure
hold on
patch('vertices',M.vertices,'faces',M.faces,...
    'facecolor','r','edgealpha',0,'facealpha',.5);
patch('vertices',S.vertices,'faces',S.faces,...
    'facecolor','b','edgealpha',0,'facealpha',.5);
lightangle(gca,-70,10)
lighting gouraud
axis equal
view(-70,10)
title('Merged')

end