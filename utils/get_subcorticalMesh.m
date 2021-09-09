function get_subcorticalMesh()

maxvertices = 8200;

directory = 'D:\2020_RiskyReplay\data\mri';
addpath D:\Toolboxes\fieldtrip-20191119\external\iso2mesh

cd(directory)
meshes = {'leftHipp.gii','rightHipp.gii','leftAmg.gii','rightAmg.gii'};
nMesh = length(meshes);

%% Downsample subcortical meshes

figure
hold on

mesh = gifti('D:\Toolboxes\spm12\canonical\orig\cortex_8196.surf.gii');
patch('vertices',mesh.vertices,'faces',mesh.faces,...
        'FaceColor',[.8 .8 .8],'edgealpha',0);

for m = 1:nMesh

    % Load
    mesh = gifti(meshes{m});

    % Downsample
    dsfactor = maxvertices/size(mesh.vertices,1)/10;
    [v,f] = meshresample(mesh.vertices,mesh.faces,dsfactor);
    
    orig = mesh;
    mesh.faces = int32(f);
    mesh.vertices = single(v);
    
    % Transform
    mesh.vertices = spm_eeg_inv_transform_points(mesh.mat, mesh.vertices);

    if m==1 || m==2
        col = 'b';
    elseif m==3 || m==4
        col = 'r';
    end

    patch('vertices',mesh.vertices,'faces',mesh.faces,...
        'EdgeColor',col,'FaceColor',col,'edgealpha',.3,'facealpha',.1);
    
    % Save
    tmp = strsplit(meshes{m},'.gii');
    save(mesh,[tmp{1} '.trans.gii']);

end
axis equal
view(-70,10)

lightangle(gca,-70,10)
lighting gouraud

%% Merge

% Merge hemispheres
meshes = {'Hipp','Amg'};
for m = 1:length(meshes)
    
    L = gifti(['left' meshes{m} '.trans.gii']);
    R = gifti(['right' meshes{m} '.trans.gii']);
    
    [v,f] = mergemesh(L.vertices,L.faces,R.vertices,R.faces);
    [v,f] = meshcheckrepair(v,f);
    
    M = L;
    M.faces = f(:,1:3);
    M.vertices = v;
    
    save(M,['both' meshes{m} '.trans.gii']);
    
end

% Add hippocampus to canonical cortical mesh
M1 = gifti('D:\Toolboxes\spm12\canonical\orig\cortex_8196.surf.gii');
M2 = gifti('bothHipp.trans.gii');

[v,f] = mergemesh(M1.vertices,M1.faces,M2.vertices,M2.faces);
[v,f] = meshcheckrepair(v,f);

M = L;
M.faces = f(:,1:3);
M.vertices = v;

save(M,'canonicalHipp.gii');

% Add amygdala to canonical mesh
M1 = gifti('D:\Toolboxes\spm12\canonical\orig\cortex_8196.surf.gii');
M2 = gifti('bothAmg.trans.gii');

[v,f] = mergemesh(M1.vertices,M1.faces,M2.vertices,M2.faces);
[v,f] = meshcheckrepair(v,f);

M = L;
M.faces = f(:,1:3);
M.vertices = v;

save(M,'canonicalAmg.gii');

% Add hippocampus AND amygdala to canonical mesh
M1 = gifti('canonicalHipp.gii');
M2 = gifti('bothAmg.trans.gii');

[v,f] = mergemesh(M1.vertices,M1.faces,M2.vertices,M2.faces);
[v,f] = meshcheckrepair(v,f);

M = L;
M.faces = f(:,1:3);
M.vertices = v;

save(M,'canonicalHippAmg.gii');

figure
patch('vertices',M.vertices,'faces',M.faces,...
        'FaceColor',[.8 .8 .8],'edgealpha',0);
axis equal
view(-155,-30)
lightangle(gca,-155,-30)
lighting gouraud

%% Check against template MRI
% 
% for a = 1:3
%     
%     figure
% 
%     % Get slices
%     V = cell(1,nMesh+1);
%     V{1} = nifti('D:\Toolboxes\spm12\canonical\orig\single_subj_T1.nii');
%     
%     [BB1,vx1] = spm_get_bbox('D:\Toolboxes\spm12\canonical\orig\single_subj_T1.nii'); % canonical dimensions
%     [BB2,vx2] = spm_get_bbox(fullfile(directory,'colin','mri','orig.nii')); % high-res dimensions
%     
%     % match bounding box
%     d = size(V{1}.dat);
%     M = V{1}.mat;
% 
%     for i = 1:3
%         tmp = ones(4,d(i));
%         tmp(1,:) = 1:d(i);
%         tmp = M(1:3,:) * tmp;
%         tmp = tmp(1,:);
%         if i==1
%             origx=tmp;
%         elseif i==2
%             origy=tmp;
%         elseif i==3
%             origz=tmp;
%         end
%     end
% 
%     for v = 1:nMesh
%         
%         tmp = strsplit(meshes{v},'.gii');
%         V{v+1} = nifti(fullfile(directory,[tmp{1} '.nii']));
%         
%         % match bounding box
%         d = size(V{v+1}.dat);
%         M = V{v+1}.mat;
%         
%         if v==1
%             for i = 1:3
%                 tmp = ones(4,d(i));
%                 tmp(1,:) = 1:d(i);
%                 tmp = M(1:3,:) * tmp;
%                 tmp = tmp(1,:);
%                 if i==1
%                     x=tmp;
%                 elseif i==2
%                     y=tmp;
%                 elseif i==3
%                     z=tmp;
%                 end
%             end
% 
%             xbounds = [findMin(origx(1),x) findMin(origx(end),x)];
%             ybounds = [findMin(origy(1),y) findMin(origy(end),y)];
%             zbounds = [findMin(origz(1),z) findMin(origz(end),z)];
%         end
%     end
%     
%     % Move through slices
%     D = size(V{1}.dat);
%     for i = 1:D(a)
%     
%         % Set up view
%         clf
%         hold on
%         axis equal
%         if a==1
%             view(90,0) % saggital (side view of head)
%         elseif a==2
%             
%         elseif a==3
%             
%         end
%         
% %         % Show cortex
% %         mesh = gifti('D:\Toolboxes\spm12\canonical\orig\cortex_8196.surf.gii');
% %         patch('vertices',mesh.vertices,'faces',mesh.faces,...
% %             'EdgeColor','c','FaceColor','c','edgealpha',.1,'facealpha',.05);
%         
%         % Show subcortical regions
%         for m = 1:length(meshes)
%             mesh = gifti(meshes{m});
%             if m==1 || m==2
%                 col = 'b';
%             elseif m==3 || m==4
%                 col = 'r';
%             end
%             patch('vertices',mesh.vertices,'faces',mesh.faces,...
%                 'EdgeColor',col,'FaceColor',col,'edgealpha',.2,'facealpha',.05);
%         end
%         
%         % Show MRI slice (T1 and subcortical regions)
%         for v = 1:length(V)
%             
%             N = V{v};
%             d = size(V{v}.dat);
%             M  = N.mat;
%             
%             if v==1
%                 I = i;
%             else
%                 coord = ones(4,1);
%                 coord(a) = i;
%                 coord = V{1}.mat(1:3,:) * coord;
%                 coord = coord(a);
%                 
%                 matchcoord = ones(4,d(a));
%                 matchcoord(a,:) = 1:d(a);
%                 matchcoord = V{v}.mat(1:3,:) * matchcoord;
%                 I = findMin(matchcoord(a,:),coord);
%             end
%             
%             if a==1
%                 [x,y,z] = ndgrid(I,1:d(2),1:d(3));
%                 f1 = squeeze(N.dat(I,:,:));
%             elseif a==2
%                 [x,y,z] = ndgrid(1:d(1),I,1:d(3));
%                 f1 = squeeze(N.dat(:,I,:));
%             elseif a==3
%                 [x,y,z] = ndgrid(1:d(1),1:d(2),I);
%                 f1 = squeeze(N.dat(:,:,I));
%             end
%             x = squeeze(x);
%             y = squeeze(y);
%             z = squeeze(z);
%             
%             x1 = M(1,1)*x+M(1,2)*y+M(1,3)*z+M(1,4);
%             y1 = M(2,1)*x+M(2,2)*y+M(2,3)*z+M(2,4);
%             z1 = M(3,1)*x+M(3,2)*y+M(3,3)*z+M(3,4);
% 
%             if ~all(f1(:)==0)
%                 if v>1
%                     x1(f1==0) = NaN;
%                     y1(f1==0) = NaN;
%                     z1(f1==0) = NaN;
%                     f1(f1==0) = NaN;
%                 end
%                 s  = surf(x1,y1,z1,f1,'facealpha',1);
%                 set(s,'EdgeColor','none')
% %                 colormap('gray')
%                 hold on
%             end
%         end
%         
%         input('Press ENTER for next slice ')
%         
%     end
% 
% end

end
    