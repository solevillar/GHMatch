function [G,r] =tri2graph(s)
n=size(s.X,1);
G=zeros(n,n);
m=size(s.TRIV,1);
for i=1:m
    i1=s.TRIV(i,1);
    i2=s.TRIV(i,2);
    i3=s.TRIV(i,3);
    if G(i1,i2)==0
        G(i1,i2)=norm([s.X(i1), s.Y(i1), s.Z(i1)]- [s.X(i2), s.Y(i2), s.Z(i2)],2);
    end
    if G(i1,i3)==0
        G(i1,i3)=norm([s.X(i1), s.Y(i1), s.Z(i1)]- [s.X(i3), s.Y(i3), s.Z(i3)],2);
    end
    if G(i2,i3)==0
        G(i2,i3)=norm([s.X(i2), s.Y(i2), s.Z(i2)]- [s.X(i3), s.Y(i3), s.Z(i3)],2);
    end
end
G=(G+G')/2;
G=sparse(G);
end
