function [h_ctx,h_iskl,h_oskl,h_slp] = checkMesh(varargin)
% Adapted from spm_eeg_inv_checkmeshes.m
% Use as checkMesh(D,val) where val is the inversion number

% OPTS
innerskullOn = false;
outerskullOn = false;
brainOn = true;
scalpOn = false;

%-Initialisation
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});
try
    %disp(D.inv{val}.mesh);
    mesh  = spm_eeg_inv_transform_mesh(eye(4), D.inv{val}.mesh);
    
    Mctx  = mesh.tess_ctx;
    Miskl = mesh.tess_iskull;
    Moskl = mesh.tess_oskull;
    Mslp  = mesh.tess_scalp;
    
    N = nifti(mesh.sMRI);
    M = N.mat;
    Mctx.vert = spm_eeg_inv_transform_points(M,Mctx.vert);

catch
    spm('alert!','Please create meshes',mfilename);
    return
end

if spm('CmdLine'), [h_ctx,h_iskl,h_oskl,h_slp] = deal([]); return; end

%-SPM Graphics figure
%--------------------------------------------------------------------------
Fgraph  = spm_figure('GetWin','Graphics'); spm_figure('Clear',Fgraph)

for a = 2%1:3
    
%     subplot(3,1,a)
    
    %-Cortex Mesh Display
    %--------------------------------------------------------------------------
    if brainOn
        h_ctx   = patch('vertices',Mctx.vert,'faces',Mctx.face,...
                    'EdgeColor','b','FaceColor','b');
    end
    hold on

    %-Inner-skull Mesh Display
    %--------------------------------------------------------------------------
    if innerskullOn
        h_iskl  = patch('vertices',Miskl.vert,'faces',Miskl.face,...
                    'EdgeColor','r','FaceColor','none');
    end

    %-Outer-skull Mesh Display
    %--------------------------------------------------------------------------
    if outerskullOn
        h_oskl  = patch('vertices',Moskl.vert,'faces',Moskl.face,...
                    'EdgeColor',[1 .5 .35],'FaceColor','none');
    end

    %-Inner Scalp Mesh Display
    %--------------------------------------------------------------------------
    if scalpOn
        h_slp   = patch('vertices',Mslp.vert,'faces',Mslp.face,...
                    'EdgeColor',[1 .7 .55],'FaceColor','none');
    end

    %-Slices
    %--------------------------------------------------------------------------
    pls = 0.05:0.1:0.9;
    N   = nifti(mesh.sMRI);
    d   = size(N.dat);
    pls = round(pls.*d(a));
    for i=1:numel(pls)
        if a==1
            [x,y,z] = ndgrid(pls(i),1:d(2),1:d(3));
            f1 = squeeze(N.dat(pls(i),:,:));
        elseif a==2
            [x,y,z] = ndgrid(1:d(1),pls(i),1:d(3));
            f1 = squeeze(N.dat(:,pls(i),:));
        elseif a==3
            [x,y,z] = ndgrid(1:d(1),1:d(2),pls(i));
            f1 = squeeze(N.dat(:,:,pls(i)));
        end
        x = squeeze(x);
        y = squeeze(y);
        z = squeeze(z);
        
        M  = N.mat;
        x1 = M(1,1)*x+M(1,2)*y+M(1,3)*z+M(1,4);
        y1 = M(2,1)*x+M(2,2)*y+M(2,3)*z+M(2,4);
        z1 = M(3,1)*x+M(3,2)*y+M(3,3)*z+M(3,4);

        if ~all(f1(:)==0)
            s  = surf(x1,y1,z1,f1,'facealpha',1);
            set(s,'EdgeColor','none')
        end
    end

    %-Fiducials
    %--------------------------------------------------------------------------
    pnt   = mesh.fid.fid.pnt;
    h_fid = plot3(pnt(:,1), pnt(:,2), pnt(:,3), '.c', 'MarkerSize', 30);
    h_fid_txt = text(pnt(:,1), pnt(:,2), pnt(:,3), mesh.fid.fid.label);

    axis image off;
    colormap('gray');
    view(-135,45);
    rotate3d on
    drawnow
    hold off
end
end