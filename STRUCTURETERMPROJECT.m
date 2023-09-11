clear all
clc
NumNodes = 8
NumElements = 8
NumSupports = 3
NumLoadJoints = 3

XY = [0,0;0,3;0,6;6,6;6,3;11.5,3;6,0;11.5,0]
M= [0.1,7*10^-4,2*10^8]
C= [1,2,1;2,3,2;3,4,3;4,5,4;2,5,5;5,6,6;5,7,7;6,8,8]
S=[1,1,1,0;7,1,1,0;8,1,1,1]
La=[1,0,0,0;2,50,0,0;3,70,0,0;4,0,-30,-22.5;5,0,0,0;6,0,0,0;7,0,0,0;8,0,0,0]
Lf=[1,0,0,0;2,0,60,60;3,0,60,60;4,0,60,-60;5,0,110.42,-9.58;6,0,50.42,-50.42;7,0,0,0;8,0,0,0]
L=La-Lf
E=zeros(NumNodes,3)
for i=1:size(S,1)
    E(S(i,1),:)=-1
    if S(i,2)==0
        E(S(i,1),1)=0
    end
   if S(i,3)==0
       E(S(i,1),2)= 0
   end
   
    if S(i,4)==0
        E(S(i,1),3)= 0
    end
end

c=1;
for i=1:NumNodes
    for j=1:3
        if E(i,j)~= -1
            E(i,j)=c
            c=c+1
             
        end
    end
    
end
   E(E(:,:)==-1)=0
NumEq= max(max(E))
KG= zeros(NumEq,NumEq);
for i=1:NumElements
     x1= XY(C(i,1),1) 
     x2= XY(C(i,2),1)
     y1= XY(C(i,1),2)
     y2= XY(C(i,2),2)
    %teta = atan((XY(C(i,1),1)-XY(C(i,2),1))/(XY(C(i,1),2)-XY(C(i,2),2)))
    Le=sqrt((x2-x1)^2+((y2-y1)^2));
    cos=(x2-x1)/Le;
    sin=(y2-y1)/Le ;
    R=[cos sin 0 0 0 0;
      -sin cos 0 0 0 0;
       0 0 1 0 0 0;  
      0   0  0 cos sin 0;
      0 0 0 -sin cos 0;
        0   0  0 0 0 1]
    k=[M(3)*M(1)/Le   0                 0         -M(3)*M(2)/Le    0                  0;
       0       12*M(3)*M(2)/Le^3 6*M(3)*M(2)/Le^2      0        -12*M(3)*M(2)/Le^3 6*M(3)*M(2)/Le^2;
       0       6*M(3)*M(2)/Le^2  4*M(3)*M(2)/Le        0        -6*M(3)*M(2)/Le^2  2*M(3)*M(2)/Le;
      -M(3)*M(2)/Le   0                 0          M(3)*M(2)/Le       0                  0; 
       0      -12*M(3)*M(2)/Le^3 -6*M(3)*M(2)/Le^2     0        12*M(3)*M(2)/Le^3 -6*M(3)*M(2)/Le^2;
       0       6*M(3)*M(2)/Le^2   2*M(3)*M(2)/Le       0        -6*M(3)*M(2)/Le^2  4*M(3)*M(2)/Le]
    K=R'*k*R;
    
    
   
    G=[E(C(i,1),:),E(C(i,2),:)]
      
        for i3=1:6 
            for j=1:6
                if G(i3)==0 || G(j)==0
                    continue
                else
                    KG(G(i3),G(j))=KG(G(i3),G(j))+K(i3,j)
                end
            end
        end

end

FNumEq=zeros(NumEq,1)

for i=1:NumNodes
   for j=1:3
    if E(i,j)==0 
    continue
    else
        FNumEq(E(i,j))=L(i,j+1) 
    end
   end
end
D=KG^-1*FNumEq

for i=1:NumElements
    d=zeros(6,1);
    %for m=1:3
        for n=1:3
            if E(C(i,1),n)==0
                continue
            else
            d(n,1)=D(E(C(i,1),n))         
            end
        end
    %end
    
    %for m=1:3
        for n=1:3
            if E(C(i,2),n)==0
                continue
            else
            d(n+3,1)=D(E(C(i,2),n))         
            end
        end
    %end
     x1= XY(C(i,1),1) 
     x2= XY(C(i,2),1)
     y1= XY(C(i,1),2)
     y2= XY(C(i,2),2)
    Le=sqrt((x2-x1)^2+((y2-y1)^2));
    cos=(x2-x1)/Le;
    sin=(y2-y1)/Le ;
    
    R=[cos sin 0 0 0 0;
      -sin cos 0 0 0 0;
       0 0 1 0 0 0;  
      0   0  0 cos sin 0;
      0 0 0 -sin cos 0;
        0   0  0 0 0 1]
    
    k=[M(3)*M(1)/Le   0                 0         -M(3)*M(2)/Le    0                  0;
       0       12*M(3)*M(2)/Le^3 6*M(3)*M(2)/Le^2      0        -12*M(3)*M(2)/Le^3 6*M(3)*M(2)/Le^2;
       0       6*M(3)*M(2)/Le^2  4*M(3)*M(2)/Le        0        -6*M(3)*M(2)/Le^2  2*M(3)*M(2)/Le;
      -M(3)*M(2)/Le   0                 0          M(3)*M(2)/Le       0                  0; 
       0      -12*M(3)*M(2)/Le^3 -6*M(3)*M(2)/Le^2     0        12*M(3)*M(2)/Le^3 -6*M(3)*M(2)/Le^2;
       0       6*M(3)*M(2)/Le^2   2*M(3)*M(2)/Le       0        -6*M(3)*M(2)/Le^2  4*M(3)*M(2)/Le]
   d=R*d;
   f(:,i)=k*d;
end

