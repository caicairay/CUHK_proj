function result=flmtdtct(density)
%% functions
nearzero=@(x,tol) abs(x)<tol;
%% control panel
coeff=0.5;
%% threshold and data size
quantum=mean(density(:))/10;
threshold=quantum*10;
sizeA=size(density);
dim=ndims(density);
%% pre-define
if dim==2
    [gx, gy] = gradient(double(density));
    [gxx, gxy] = gradient(gx);
    [gyx, gyy] = gradient(gy);
    tolx=mean(abs(gx(:)));
    toly=mean(abs(gy(:)));
    tol=coeff*hypot(tolx,toly);
    grad=cat(3,gx,gy);
    Hx=cat(3,gxx,gxy);
    Hy=cat(3,gyx,gyy);
    H=cat(4,Hx,Hy);
    grad_p=zeros(sizeA(1),sizeA(2),dim);
else
    [gx,gy,gz] = gradient(double(density));
    [gxx, gxy,gxz] = gradient(gx);
    [gyx, gyy,gyz] = gradient(gy);
    [gzx, gzy,gzz] = gradient(gz);
    tolx=mean(abs(gx(:)));
    toly=mean(abs(gy(:)));
    tolz=mean(abs(gy(:)));
    tol=coeff*sqrt(tolx^2+toly^2+tolz^2);
    grad=cat(4,gx,gy,gz);
    Hx=cat(4,gxx,gxy,gxz);
    Hy=cat(4,gyx,gyy,gyz);
    Hz=cat(4,gzx,gzy,gzz);
    H=cat(5,Hx,Hy,Hz);
    grad_p=zeros(sizeA(1),sizeA(2),sizeA(3),dim);
end
eig_val=zeros(size(H));
eig_vec=zeros(size(H));
result=false(sizeA);
%% detection
if dim==2
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
else
    for i=1:sizeA(1)
        for j=1:sizeA(2)
            for k=1:sizeA(3)
                [V,D]=eig(squeeze(H(i,j,k,:,:)));
                eig_val(i,j,k,:,:)=D;
                eig_vec(i,j,k,:,:)=V;
                grad_p(i,j,k,1)=squeeze(eig_vec(i,j,k,:,1))'*squeeze(grad(i,j,k,:));
                grad_p(i,j,k,2)=squeeze(eig_vec(i,j,k,:,2))'*squeeze(grad(i,j,k,:));
                
                G_test=nearzero(max(grad_p(i,j,k,:)),tol);
                result(i,j,k)=G_test && eig_val(i,j,k,dim-1,dim-1)<0 && density(i,j,k)>threshold;
            end
        end
    end
end
%% double check
%result=true(sizeA);
if dim==2
    for i=1:sizeA(1)
        for j=1:sizeA(2)
            if result(i,j)
                gxp=grad_p(i,j,1);
                gyp=grad_p(i,j,2);
                doublecheck=true;
                ang=atan(abs(gyp/gxp))/pi*180;
                step=zeros(2);
                if ang<60 && ang>30
                    step=[sign(gxp),sign(gyp)];
                else
                    [~,idx]=max([abs(gxp),abs(gyp)]);
                    step(idx)=sign(gxp)*(idx==1)+sign(gyp)*(idx==2);
                end
                for ii=min(i-step(1),i):max(i-step(1),i)
                    for jj=min(j-step(2),j):max(j-step(2),j)
                        if ii>0 && jj>0 && ii<sizeA(1) && jj<sizeA(2)
                            if abs(grad_p(i,j))>abs(grad_p(ii,jj))
                                if result(ii,jj)==true
                                doublecheck=false;
                                else
                                result(ii,jj)=true;
                                end
                            end
                        end
                    end
                end
                
%                 ang=atan(abs(gyp/gxp))/pi*180;
%                 step=zeros(2);
%                 if ang<60 && ang>30
%                     step=[sign(gxp),sign(gyp)];
%                 else
%                     [~,idx]=max([abs(gxp),abs(gyp)]);
%                     step(idx)=sign(gxp)*(idx==1)+sign(gyp)*(idx==2);
%                 end
%                 for ii=min(i,i-step(1)):max(i,i-step(1))
%                     for jj=min(j,j-step(2)):max(j,j-step(2))
%                         if ii>0 && jj>0 && ii<sizeA(1) && jj<sizeA(2)
%                             if abs(grad_p(i,j))<abs(grad_p(ii,jj)) || grad_p(i,j)==0 % && result(ii,jj)~=true
%                                 doublecheck=true;
%                             end
%                         end
%                     end
%                 end
                result(i,j)=doublecheck;
            end
        end
    end
else
    for i=1:sizeA(1)
        for j=1:sizeA(2)
            for k=1:size(3)
                if result(i,j,k)
                    gxp=grad_p(i,j,k,1);
                    gyp=grad_p(i,j,k,2);
                    gzp=grad_p(i,j,k,3);
                    gp=[gxp,gyp,gzp];
                    step=zeros(3);
                    doublecheck=false;
                    [~,basis]=max([abs(gxp),abs(gyp),abs(gzp)]);
                    step(basis)=sign(gp(basis));
                    for n=1:3
                        if n~=basis
                            ang=atan(abs(gp(n)/gp(basis)))/pi*180;
                            if ang<60 && ang>30
                                step(n)=sign(gp(n));
                            else
                                step(n)=0;
                            end
                        end
                    end
                    for ii=min(i,i-step(1)):max(i,i-step(1))
                        for jj=min(j,j-step(2)):max(j,j-step(2))
                            for kk=min(k,k-step(3)):max(k,k-step(3))
                                if ii>0 && jj>0 && kk>0 && ii<sizeA(1) && jj<sizeA(2) && kk<sizeA(3)
                                    if grad_p(i,j,k)*grad_p(ii,jj,kk)<0 && result(ii,jj,kk)~=true || grad_p(i,j,k)==0
                                        doublecheck=true;
                                    end
                                end
                            end
                        end
                    end
                    result(i,j,k)=doublecheck;
                end
            end
        end
    end
end
end