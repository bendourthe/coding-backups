function [Mtotal,Xgtotal,J] = stl2sm(filename,rho)
% Copyright 2011 The MathWorks, Inc.
  
  % stl2sm
  % Built SimMechanics block based on stl file.
  % The function calculates the body mass, center of mass and inertia tensor 
  % and creates a SimMechanics body.
  
  %%
  run(filename);
  
  %% Init
  numtriangle = numel(triangle);
  Ccan=1/60*(1/2*eye(3)+1/2*ones(3,3))*rho;
  Xgtotal = zeros(3,1);
  V=zeros(numtriangle,1);
  M=zeros(numtriangle,1);
  C=zeros(3,3,numtriangle);
  Xg=zeros(3,1,numtriangle);
  
  %% Loop (origin at [0;0;0]
  for i = 1:numtriangle
    A = reshape(struct2array(triangle(i)),3,3);
    V(i) = 1/6 * det(A); %Volume
    M(i) = rho * V(i);   %Mass
    Xg(:,:,i) = mean([zeros(3,1) A],2); %Center of Mass
    C(:,:,i) = det(A)*A*Ccan*A';        %Covariance
  end
  
  Mtotal = sum(M);
  Ctotal = sum(C,3);
  %Vtotal = sum(V);
  for i=1:numtriangle
    Xgtotal = Xgtotal + Xg(:,:,i)*M(i);
  end
  Xgtotal = Xgtotal/Mtotal;
  
  %Translate to Xgtotal (center of mass)
  dx = zeros(3,1) - Xgtotal;
  Ctotal = translateC(Ctotal,dx,Xgtotal,Mtotal);
  
  %Tensor of Inertia
  J=eye(3)*trace(Ctotal)-Ctotal;
end