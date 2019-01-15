%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo of "Fuzzy region-based active contours driven by hybrid fitted energy 
% with local and global information for image segmentation"
% Jiangxiong Fang
% East China University of Technology&&Nanchang University, Nanchang, China
% 6th, Jan., 2019
% Email: fangchj2002@163.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,e,deltaF]= FRAGL_v1(Img,u0,Ksigma,lambda1,lambda2,alpha,belta,m1,m2,g)

  epsilon = 1;
  timestep = 0.01;
  mu = 1.5;% level set regularization term, please refer to "Chunming Li 
  nu = 5;%length term


  u1 = u0.^2;
  u2 = (1-u0).^2;
  
  Iu1 = Img.*u1;
  Iu2 = Img.*u2;
  
  c1 = sum(sum(Iu1))/sum(sum(u1));
  c2 = sum(sum(Iu2))/sum(sum(u2));  
  
  Ku1 = imfilter(u1,Ksigma,'replicate'); 
  Ku2 = imfilter(u2,Ksigma,'replicate'); 
         
  KI1 = imfilter(Iu1,Ksigma,'replicate');
  KI2 = imfilter(Iu2,Ksigma,'replicate');
         
  fo = sum(sum(KI1))/sum(sum(Ku1));
  fb = sum(sum(KI2))/sum(sum(Ku2));

  F1_old = (Img-belta*c1-alpha*fo).^2.*u1;
  F2_old = (Img-belta*c2-alpha*fb).^2.*u2;
  
  un= 1./(1+(lambda1*((Img-alpha*fo-belta*c1).^2))./(lambda2*((Img-alpha*fb-belta*c2).^2)));
  
  un1 = un.^2;
  un2 = (1-un).^2;
  
  delta_u1 = un1-u1;
  delta_u2 = un2-u2;
 
  delta_Ku1 = imfilter(delta_u1,Ksigma,'replicate'); 
  delta_Ku2 = imfilter(delta_u2,Ksigma,'replicate'); 
  
  fs1 = alpha*(delta_Ku1./(delta_Ku1+Ku1)).*(Img-fo)+belta*(delta_u1./(delta_u1+u1)).*(Img-c1);
  fun1 = fs1.^2.*u1.*g;
  
  fs2 = alpha*(Ku1./(delta_Ku1+Ku1)).*(Img-fo)+belta*(u1./(delta_u1+u1)).*(Img-c1);
  fun2 = fs2.^2.*delta_u1.*g;
  
  fs3 = alpha*(delta_Ku2./(delta_Ku1+Ku2)).*(Img-fb)+belta*(delta_u2./(delta_u2+u2)).*(Img-c2);
  fun3 = fs3.^2.*u2.*g;
  
  fs4 = alpha*(Ku2./(delta_Ku2+Ku2)).*(Img-fb)+belta*(u2./(delta_u2+u2)).*(Img-c2);
  fun4 = fs4.^2.*delta_u2.*g;

  deltaF = fun1+fun2+fun3+fun4;
  
  idx = find(deltaF<0);
  u0(idx)=un(idx); 
  
  e =sum(sum(F1_old+F2_old));
  deltaF = sum(sum(deltaF));
  u = u0;   
  u = imfilter(u,Ksigma,'replicate'); 
  Delta = Dirac(u,epsilon);
  K=curvature_central(u);
  P=mu*(4*del2(u) - K);
  L=nu.*Delta.*K;%length term   
  u = u+timestep*(m1*L+m2*P);  
  u = imfilter(u,Ksigma,'replicate'); 
end

function K = curvature_central(u)
   [bdx,bdy]=gradient(u);
   mag_bg=sqrt(bdx.^2+bdy.^2)+1e-10;
   nx=bdx./mag_bg;
   ny=bdy./mag_bg;
   [nxx,nxy]=gradient(nx);
   [nyx,nyy]=gradient(ny);
   K=nxx+nyy;
end

function f = Dirac(x, epsilon)
   f=(epsilon/pi)./(epsilon^2.+x.^2);
end
