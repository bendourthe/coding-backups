function [V_reg, mnorm] = Muscle_AlignRegion(F, V, F_reg, V_reg, wbar, method)
%Muscle_AlignRegion  Align database muscle region to new bone mesh.
%
%   Syntax:
%    [V_reg, mnorm] = Muscle_AlignRegion(F, V, F_reg, V_reg, wbar, method)
%    [V_regs, mnorms] = Muscle_AlignRegion(F, V, F_regs, V_regs, wbar, ...
%                                          method)
%
%   Input:
%    F:      J-by-3 array defining the new bone mesh. The rows correspond
%            to different triangles and the columns correspond to the three
%            vertices that make up each triangle. The elements are row
%            indices into V.
%    V:      K-by-3 array defining vertices in the new bone mesh. The rows
%            correspond to different vertices and the columns correspond to
%            X-, Y- and Z-coordinates. The elements are coordinate values.
%    F_reg:  L-by-3 array defining the muscle attachment region. The rows
%            correspond to different triangles and the columns correspond
%            to the three vertices that make up each triangle. The elements
%            are row indices into V_reg.
%    V_reg:  M-by-3 array defining vertices in the muscle attachment
%            region. The rows correspond to different vertices and the
%            columns correspond to X-, Y- and Z-coordinates. The elements
%            are coordinate values.
%    wbar:   Logical indicating whether or not a wait bar should be shown.
%            Optional, defaults to true.
%    method: String indicating the type of transformation to be used when
%            aligning the regions to the bone. Will be used as input for
%            ICP_RegisterMatchedDatasets.m; see that function's help text
%            for possible transformations. Optional, defaults to 'rigid'.
%    F_regs: N-element column cell concatenation of F_reg. The elements
%            correspond to different muscle attachment regions.
%    V_regs: N-element column cell concatenation of V_reg. The elements
%            correspond to different muscle attachment regions.
%
%   Output:
%    V_reg:  M-by-3 array defining vertices in the muscle attachment
%            region. The rows correspond to different vertices and the
%            columns correspond to X-, Y- and Z-coordinates. The elements
%            are coordinate values.
%    mnorm:  3-element row vector defining the mean triangle normal of the
%            muscle attachment region, weighted by triangle areas. The
%            elements are X-, Y- and Z-coordinates.
%    V_regs: N-element column cell concatenation of V_reg. The elements
%            correspond to different muscle attachment regions.
%    mnorms: N-by-3 concatenation of mnorm. The rows correspond to
%            different muscle attachment regions.
%
%   Effect: This function will align a database muscle attachment region to
%   a new bone mesh. If the region lies fully outside the bone, it will
%   first be moved closer to the bone and rotated so that its mean normal
%   is parallel to the local normal of the bone surface. Next, the region
%   will be registered onto the bone mesh by use of the ICP algorithm. For
%   the ICP registration, bone vertices will only be considered if their
%   normal points in the same direction as the region mean normal.
%
%   Dependencies: TRI_Normals.m
%                 TRI_VertexNormals.m
%                 VectorNorms.m
%                 TRI_VertexAreas.m
%                 TRI_MeanNormal.m
%                 ICP_MatchDatasets.m
%                 VectorAligningRotation.m
%                 ICP_RegisterDatasets.m
%
%   Known parents: Muscle_ProjectDatabaseRegions.m

%Created on 21/08/2007 by Ward Bartels.
%WB, 19/07/2009: Added optional wait bar.
%WB, 23/07/2009: Added support for multiple muscle regions.
%WB, 15/04/2010: Changed pre-alignment procedure. It is now made of win.
%Stabile, fully functional.


%Parameter: minimum fraction of points to be used for pre-aligning region
minfrac = .2;

%Set defaults for inputs
if nargin<5, wbar = true; end
if nargin<6, method = 'rigid'; end

%Handle multiple regions if necessary
if iscell(V_reg)
    
    %Create wait bar if necessary <<GUI_waitbar.m>>
    numreg = numel(V_reg);
    if wbar
        verts = cumsum(cellfun(@(x) size(x, 1), V_reg)); %vertices "done"
        wb = waitbar(0, sprintf('Aligning muscle attachment group 1 of %d', numreg), ...
                         'Name', 'Aligning attachments...');
    end
    
    %Align all regions
    mnorm = zeros(numreg, 3);
    for ind = 1:numreg;
        
        %Align region <<Muscle_SelectByRegionBorder.m>>
        [V_reg{ind}, mnorm(ind,:)] = ...
            Muscle_AlignRegion(F, V, F_reg{ind}, V_reg{ind}, false, method);
        
        %Update wait bar if necessary
        if wbar
            waitbar(verts(ind)/verts(end), wb, ...
                    sprintf('Aligning muscle attachment group %d of %d', ...
                            min(ind+1, numreg), numreg));
        end
    end
    
    %Close wait bar if necessary
    if wbar, delete(wb); end
    
    return
end

%Calculate bone vertex normals and vertex areas <<TRI_Normals.m>>
%<<TRI_VertexNormals.m>> <<VectorNorms.m>> <<TRI_VertexAreas.m>>
tnorms = TRI_Normals(F, V, false);
normals = TRI_VertexNormals(F, V, true, tnorms);
W_mod = TRI_VertexAreas(F, VectorNorms(tnorms), true);

%Calculate region mean normal and vertex areas <<TRI_Normals.m>>
%<<TRI_MeanNormal.m>> <<VectorNorms.m>> <<TRI_VertexAreas.m>>
tnorms_reg = TRI_Normals(F_reg, V_reg, false);
mnorm = TRI_MeanNormal(tnorms_reg);
W_reg = TRI_VertexAreas(F_reg, VectorNorms(tnorms_reg), true);

%Use only bone vertices with normal facing the same way as region
dotprod = normals*mnorm.';
use = dotprod>0;
V = V(use,:);

%Take normal direction into account for weights
W_mod = W_mod(use,:).*dotprod(use,:);

%Project bone vertices on plane perpendicular to mnorm
PR = null(mnorm);
V_pr = V*PR;

%Get attachment region border edge loops
%<<TRI_DetermineBorderEdges.m>> <<Contour_Loops.m>>
loops = Contour_Loops(TRI_DetermineBorderEdges(F_reg));

%Isolate vertices projected inside region border edge loops
in = cellfun(@(l) inpolygon(V_pr(:,1), V_pr(:,2), ...
                             V_reg(l,:)*PR(:,1), V_reg(l,:)*PR(:,2)), ...
              loops, 'UniformOutput', false);
in = any([in{:}], 2);

%Only pre-align region if enough points were intercepted
if sum(in)>ceil(minfrac*size(V_reg, 1))
    
    %Translate region along mnorm so its peak coincides with that of
    %isolated vertices
    V_reg = V_reg+ones(size(V_reg, 1), 1)*((max(V(in,:)*mnorm.')-...
                                            max(V_reg*mnorm.'))*mnorm);
end

%Register region onto bone <<ICP_RegisterDatasets.m>>
[R, ignoble, V_reg] = ...
    ICP_RegisterDatasets(V_reg, V, W_reg, W_mod, ...
                         wbar, true, zeros(0, 1), method);
mnorm = mnorm*R;


% %Original code (wb20100415), begin:
% 
% %Consider only bone vertices facing the same way as region mean normal
% vertices = find(normals*mnorm.'>0);
% 
% %Find closest bone vertex for each region vertex <<ICP_MatchDatasets.m>>
% ind = vertices(ICP_MatchDatasets(V_reg, V(vertices,:)));
% displ = V_reg-V(ind,:);
% 
% %If region lies fully outside bone, move it closer and rotate it
% if all(dot(displ, normals(ind,:), 2)>0)
%     
%     %Determine new normal direction parallel to bone normals
%     ind = unique(ind);
%     newnorm = W_mod(ind,:).'*normals(ind,:);
%     newnorm = newnorm/norm(newnorm);
% 
%     %Create orthonormal rotation matrix <<VectorAligningRotation.m>>
%     R = VectorAligningRotation(mnorm, newnorm);
%     
%     %Rotate region to new normal direction, and displace to bone mesh
%     centroid = W_reg.'*V_reg/sum(W_reg);
%     T = centroid-centroid*R-W_reg.'*displ/sum(W_reg);
%     V_reg = V_reg*R+T(ones(size(V_reg, 1), 1),:);
%     mnorm = newnorm;
%     
%     %Update list of vertices considered for registration
%     vertices = normals*mnorm.'>0;
% end
% 
% %Register region onto bone <<ICP_RegisterDatasets.m>>
% [R, ignoble, V_reg] = ...
%     ICP_RegisterDatasets(V_reg, V(vertices,:), W_reg, W_mod(vertices), ...
%                          wbar, true, zeros(0, 1), method);
% mnorm = mnorm*R;
% 
% %Original code (wb20100415), end.