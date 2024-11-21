function [t,x,v,a,at] = TimeHistory(ag,Tn,zeta,g)

wn = 2*pi / Tn ;
[j,~] = size(ag) ; 
R = zeros(j , 3 ) ;
R (1,1) = -ag(1,2)*g ; 

for i=2:j
     dt = ag(i,1) - ag (i-1,1);
     A1 = 1 + zeta * wn * dt + (1/4)* wn^2 * dt^2 ;
     B1 = A1 - 1 ;
     C2 = 2 * zeta * wn +  wn^2 * dt ;
     
     R(i,1) = ( -(ag(i,2))*g - B1*R(i-1,1) - C2 * R(i-1,2) - wn^2*R(i-1,3))/A1 ;
     
     R(i,2) = R(i-1,2) + 0.5*(R(i-1)+ R(i,1))*dt ;
     
     R(i,3) = R(i-1,3)+ R(i-1,2)*dt  + 0.25*(R(i-1)+ R(i,1))*dt^2 ;  
end

   t = ag(:,1);
   x = R(:,3);
   v = R(:,2);
   a = R(:,1);
   at = a+ag(:,2)*g;

end




