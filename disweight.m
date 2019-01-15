function wght = disweight(rad)
   if rad == 0
       wght =0;
   else
       dist = ones(2*rad+1,2*rad+1);
       
       for m=1:2*rad+1
           for n=1:2*rad+1
               dist(m,n) = 1/(1+sqrt((m-rad)^2+(n-rad)^2));
           end
       end
       wght = dist/sum(sum(dist));
   end    
end