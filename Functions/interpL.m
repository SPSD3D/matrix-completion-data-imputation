function L = interpL( x )

indsz=find(x==0);
   L=eye(length(x));
   w={};
   for i=1:length(indsz)    
   
       count1=0;
       w{i}=count1;
       docontinue=true;
       pointer=1;
       while docontinue
            if (((i+pointer)<=length(indsz)))
                if (indsz(i+pointer)==indsz(i+pointer-1)+1)
                    count1=count1+1;  
                    pointer=pointer+1;
                    docontinue=true;
                else
                    docontinue=false; 
                end
            else
                docontinue=false; 
            end
       end

       count2=0;
       docontinue=true;
       pointer=1;
       while docontinue
            if (((i-pointer)>=1))
                if (indsz(i-pointer+1)==indsz(i-pointer)+1)
                    count2=count2+1;  
                    pointer=pointer+1;
                    docontinue=true;
                else
                    docontinue=false; 
                end
            else
                docontinue=false; 
            end
       end
       w{i}=[count2 count1];
   end
    

    for i=1:length(w)
        j=indsz(i);
        L(j,j)=0.0;
         L(j,j-1- w{i}(1))=(w{i}(2)+1)/(w{i}(1)+w{i}(2)+2);
         L(j,j+1+ w{i}(2))=(w{i}(1)+1)/(w{i}(1)+w{i}(2)+2);
    end

end

