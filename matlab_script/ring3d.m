%% imports
import kde2d.*

%% functions
nearzero=@(x,tol) abs(x)<tol;

%% control panel
coeff=.3;
skip=1;

%% points generator
N1 = 20;                                                           
R1 = 100. + 2e-10*(rand(1, N1)-0.5);                                
A1 = linspace(0, 2*pi, N1);                                         

Ring1 = R1.*cos(A1);
Ring2 = R1.*sin(A1);
ring= [Ring1',Ring2'];

%% kde
[bandwidth,density,X,Y]=kde2d(ring);
% fid = fopen('Mymatrix.txt','wt');
% 
% for ii = 1:size(density,1)
%     fprintf(fid,'%g\t',density(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid);

%% threshold and data size
quantum=mean(density(:))/10;
threshold=quantum*15;
sizeA=size(density);
dim=ndims(density);
%% parameters
[gx, gy] = gradient(double(density));
[gxx, gxy] = gradient(gx);
[gyx, gyy] = gradient(gy);

tolx=mean(abs(gx(:)));
toly=mean(abs(gy(:)));
tol=coeff*hypot(tolx,toly);
result=false(sizeA);

grad=cat(3,gx,gy);
Hx=cat(3,gxx,gxy);
Hy=cat(3,gyx,gyy);
H=cat(4,Hx,Hy);
eig_val=zeros(size(H));
eig_vec=zeros(size(H));
grad_p=zeros(sizeA(1),sizeA(2),dim);

%% first check
for i=1:sizeA(1)
    for j=1:sizeA(2)
        [V,D]=eig(squeeze(H(i,j,:,:)));
        eig_val(i,j,:,:)=D;
        eig_vec(i,j,:,:)=V;
        grad_p(i,j,:)=squeeze(eig_vec(i,j,:,1:end-1)).'*squeeze(grad(i,j,:));
        G_test=nearzero(max(grad_p(i,j,:)),tol);
        result(i,j)=G_test && eig_val(i,j,dim-1,dim-1)<0 && density(i,j)>threshold;
    end 
end

%% double check
for i=1:sizeA
    for j=1:sizeA
        if result(i,j)
            gxp=grad_p(i,j,1);
            gyp=grad_p(i,j,2);
            doublecheck=false;
            
%             for ii=min(i,i-sign(gxp)):max(i,i-sign(gxp))
%                 for jj=min(j,j-sign(gyp)):max(j,j-sign(gyp))
%                     if ii>0 && jj>0 && ii<sizeA(1) && jj<sizeA(2)
%                         if grad_p(i,j)*grad_p(ii,jj)<0 && result(ii,jj)~=true || grad_p(i,j)==0
%                             doublecheck=true;
%                         end
%                     end
%                 end
%             end
%             result(i,j)=doublecheck;

            ang=atan(abs(gyp/gxp))/pi*180;
            k=zeros(2);
            if ang<60 && ang>30
                k=[sign(gxp),sign(gyp)];
            else
                [~,idx]=max([abs(gxp),abs(gyp)]);
                k(idx)=sign(gxp)*(idx==1)+sign(gyp)*(idx==2);
            end
            for ii=min(i,i-k(1)):max(i,i-k(1))
                for jj=min(j,j-k(2)):max(j,j-k(2))
                    if ii>0 && jj>0 && ii<sizeA(1) && jj<sizeA(2)
                        if grad_p(i,j)*grad_p(ii,jj)<0 && result(ii,jj)~=true || grad_p(i,j)==0
                        	doublecheck=true;
                        end
                    end
                end
            end
            result(i,j)=doublecheck;
        end
    end
end

%% figures
% figure(1)
% mesh(X(1:skip:end,1:skip:end),Y(1:skip:end,1:skip:end),density(1:skip:end,1:skip:end))
% hold on
% plot3(X(result),Y(result),density(result),'.','MarkerSize',15)
% hold off
figure
contourf(X,Y,density,50)
hold on
scatter(X(result),Y(result))%,'.','MarkerSize',15)
axis equal
axis tight
hold off