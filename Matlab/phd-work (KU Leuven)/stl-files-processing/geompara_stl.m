function [Vtotal,Mtotal,J] = geompara_stl(Face,rho,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geopara_stl: calculates different geometrical parameters from a stl
% file:
%   - surf_area: total area of the surface of the stl
%   - vol: total volume of the stl
%   - cm: center of mass of the stl
%   - J: tensor of inertia of the stl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization

  filename = stl2m();
  handles.filename = filename;
  set(handles.file,'String',filename(1:end-2));
  set(handles.pb_create,'Enable','on');
  % Update handles structure
  guidata(hObject, handles);

contents = cellstr(get(handles.pu_material,'String'));
  material = contents{get(handles.pu_material,'Value')};
  density = handles.d(material);
  
  [M,Xg,J] = stl2sm(handles.filename(1:end-2),density);

  numtriangle = numel(Face);
  Ccan=1/60*(1/2*eye(3)+1/2*ones(3,3))*rho;
  CM = zeros(3,1);
  V=zeros(numtriangle,1);
  M=zeros(numtriangle,1);
  C=zeros(3,3,numtriangle);
  Xg=zeros(3,1,numtriangle);
  
  %% Loop (origin at [0;0;0]
  for i = 1:numtriangle
    A = reshape(struct2array(Face(i)),3,3);
    V(i) = 1/6 * det(A); %Volume
    M(i) = rho * V(i);   %Mass
    Xg(:,:,i) = mean([zeros(3,1) A],2); %Center of Mass
    C(:,:,i) = det(A)*A*Ccan*A';        %Covariance
  end
  
  Vtotal = sum(V);
  Mtotal = sum(M);
  Ctotal = sum(C,3);
  %Vtotal = sum(V);
  for i=1:numtriangle
    CM = CM + Xg(:,:,i)*M(i);
  end
  CM = CM/Mtotal;
  
  % Translate to the center of mass
  dx = zeros(3,1) - CM;
  Ctotal = translateC(Ctotal,dx,CM,Mtotal);
  
  % Tensor of Inertia
  J=eye(3)*trace(Ctotal)-Ctotal;
end