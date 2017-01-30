% color colordeconv2 on 8-bit RGB image: returns H,E,R channels from an HE image
% inputs: -path to he image
%         -imagetype

function [H,E,R] = colordeconv2(RGB,imagetype)

mx=zeros(3,1);
my=zeros(3,1);
mz=zeros(3,1);

if strcmp(imagetype,'PSL'),
    
    %H
    mx(1,1)= 0.8249;
    my(1,1)= 0.8040;
    mz(1,1)= 0.5969;
    
    %E
    mx(2,1)= 0.0825;
    my(2,1)= 0.5454;
    mz(2,1)= 0.3680;
    
    mx(3,1)= 0.0;
    my(3,1)= 0.0;
    mz(3,1)= 0.0;
    
elseif strcmp(imagetype,'NUH'),
    
    %H
    mx(1,1)= 0.650;
    my(1,1)= 0.704;
    mz(1,1)= 0.286;
    
    %E
    mx(2,1)= 0.072;
    my(2,1)= 0.990;
    mz(2,1)= 0.105;
    
    mx(3,1)= 0.0;
    my(3,1)= 0.0;
    mz(3,1)= 0.0;
    
else
    
    error(['Unknown image type : ' imagetype]);
    
end

%convert stain vectors to 'q' vector

cx=zeros(3,1);
cy=zeros(3,1);
cz=zeros(3,1);

l=zeros(3,1);

%normalize length
for i=1:3,
    
    l(i,1)=sqrt(mx(i,1)^2+my(i,1)^2+mz(i,1)^2);%0 for i=3
    
    if(l(i,1)>0),
    
        cx(i,1)=mx(i,1)/l(i,1);
        cy(i,1)=my(i,1)/l(i,1);
        cz(i,1)=mz(i,1)/l(i,1);
        
    end
    
end

%first and second vectors never empty, need to compute third vector

if cx(1,1)^2+cx(2,1)^2 <= 1,
    
    cx(3,1) = sqrt(1 - cx(1,1)^2 - cx(2,1)^2);

end

if cy(1,1)^2+cy(2,1)^2 <= 1,
    
    cy(3,1) = sqrt(1 - cy(1,1)^2 - cy(2,1)^2);

end

if cz(1,1)^2+cz(2,1)^2 <= 1,
    
    cz(3,1) = sqrt(1 - cz(1,1)^2 - cz(2,1)^2);

end

l3=sqrt(cx(3,1)^2+cy(3,1)^2+cz(3,1)^2);
cx(3,1)=cx(3,1)/l3;
cy(3,1)=cy(3,1)/l3;
cz(3,1)=cz(3,1)/l3;

%matrix invertion

a=cy(2,1) - cx(2,1) * cy(1,1) / cx(1,1);
v=cz(2,1) - cx(2,1) * cz(1,1) / cx(1,1);
c=cz(3,1) - cy(3,1) * (v/a) + cx(3,1) * ((v/a) * (cy(1,1) / cx(1,1)) - (cz(1,1) / cx(1,1)));

q=zeros(9,1);

q(3,1)=(-cx(3,1)/cx(1,1)-cx(3,1)/a*cx(2,1)/cx(1,1)*cy(1,1)/cx(1,1)+cy(3,1)/a*cx(2,1)/cx(1,1))/c;
q(2,1)=-q(3,1)*v/a-cx(2,1)/(cx(1,1)*a);
q(1,1)=1/cx(1,1)-q(2,1)*cy(1,1)/cx(1,1)-q(3,1)*cz(1,1)/cx(1,1);
q(6,1)=(-cy(3,1)/a+cx(3,1)/a*cy(1,1)/cx(1,1))/c;
q(5,1)=-q(6,1)*v/a+1/a;
q(4,1)=-q(5,1)*cy(1,1)/cx(1,1)-q(6,1)*cz(1,1)/cx(1,1);
q(9,1)=1/c;
q(8,1)=-q(9,1)*v/a;
q(7,1)=-q(8,1)*cy(1,1)/cx(1,1)-q(9,1)*cz(1,1)/cx(1,1);

%log transform the 8-bit RGB image

I=double(RGB);

lI=-((255*log((I+1)/255))/log(255));

%H channel
H=exp(-((q(1,1)*lI(:,:,1) + q(2,1)*lI(:,:,2) + q(3,1)*lI(:,:,3))-255)*log(255)/255);
H=min(H,255);
H=floor(H+0.5);

%E channel 
E=exp(-((q(4,1)*lI(:,:,1) + q(5,1)*lI(:,:,2) + q(6,1)*lI(:,:,3))-255)*log(255)/255);
E=min(E,255);
E=floor(E+0.5);

%R channel
R=exp(-((q(7,1)*lI(:,:,1) + q(8,1)*lI(:,:,2) + q(9,1)*lI(:,:,3))-255)*log(255)/255);
R=min(R,255);
R=floor(R+0.5);

end

