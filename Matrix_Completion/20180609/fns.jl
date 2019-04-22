function myhist(data, min, max, nbins)
  N = length(data)             # How many elements in the input vector 'data' ?
  delta = (max-min)/nbins      # Bin size is inferred here from the maximal, minimal, and bin number
  out = zeros(nbins)           # Let's initialize the output data structures for the bin count
  bin = zeros(nbins)           # and for the bin centres...

  start = min                  # Left edge
  for k=1:nbins
    stop   = start + delta   # Right edge
    out[k] = length(find((data .>= start) & (data .< stop))) # Count how many elements are between left and right
    bin[k] = start + delta/2. # Centre of the bin
    start  = stop            # New left edge
   end
   return out, bin
  end

  function myhist(data,bin)
    N = length(data)             # How many elements in the input vector 'data' ?
    delta = (max-min)/nbins      # Bin size is inferred here from the maximal, minimal, and bin number
    out = zeros(nbins)           # Let's initialize the output data structures for the bin count
    bin = zeros(nbins)           # and for the bin centres...

    start = min                  # Left edge
    for k=1:nbins
      stop   = start + delta   # Right edge
      out[k] = length(find((data .>= start) & (data .< stop))) # Count how many elements are between left and right
      bin[k] = start + delta/2. # Centre of the bin
      start  = stop            # New left edge
     end
     return out, bin
    end


function interpL(x,indsz=find(x->x==0,x))
    #indsz=find(x->x==0,x)
    #indsz=find( x->x == "NA",x)
    emptyLen=length(indsz);
    L=sparse([],[],Float64[],length(x),length(x))
    w=zeros(emptyLen,2)
    #w=eye(emptyLen,2)
    for i=1:emptyLen
       count1=0;
       #w{i}=count1;
       docontinue=true;
       pointer=1;
       while docontinue
            if (((i+pointer)<=length(indsz)))
                if (indsz[i+pointer]==indsz[i+pointer-1]+1)
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
                if (indsz[i-pointer+1]==indsz[i-pointer]+1)
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
       w[i,1]=count2;
       w[i,2]=count1;
       end
     for i=1:emptyLen
          j=indsz[i];
          L[j,j]=1.0;
          ind1=trunc(Int,j-1- w[i,1])
          ind2=trunc(Int,j+1+ w[i,2])
          L[j,ind1]=-(w[i,2]+1)/(w[i,1]+w[i,2]+2);
          L[j,ind2]=-(w[i,1]+1)/(w[i,1]+w[i,2]+2);
     end
return L
end


function regionGrowLeft(x,i)
    z=0;
    while ((x[i]-1==x[i]) &&  i>0)
        i=i-1;
        z=i;
    end
    return z
end





function regionGrowRight(x,i)
    z=0;
    while ((x[i]+1==x[i]) && i<=length(x))
        i=i+1;
        z=i
    end
    return z
end



function findPositionLeftEqualTo(x,i,value)
    z=0;
    while ((x[i]-1!=value) &&  i>0)
        i=i-1;
        z=i;
    end
    return z
end

function findPositionRightEqualTo(x,i,value)
    z=0;
    while ((x[i]+1!=value) && i<=length(x))
        i=i+1;
        z=i;
    end
    return z
end



function interpLC(x,indsz=find(x->x==0,x))


    x=Omega
    indsz=find(x->x==0,x)
    A=zeros(length(x))
    A[indsz]=1
    for i=1:length(x)

        if i!=1
            k=regionGrowLeft(x,i)
        else
            k=i;
        end
        if i!=1
            k=regionGrowRight(x,i)
        else
            k=i;
        end


    end




    L=sparse([],[],Float64[],length(x),length(x))
    w=zeros(length(indsz),2)

    for i=1:length(indsz)
       count1=0;
       docontinue=true;
       pointer=1;
       while docontinue
            if (((i+pointer)<=length(indsz)))
                if (indsz[i+pointer]==indsz[i+pointer-1]+1)
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
                if (indsz[i-pointer+1]==indsz[i-pointer]+1)
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
       w[i,1]=count2;
       w[i,2]=count1;
       end
     for i=1:length(indsz)
          j=indsz[i];
          L[j,j]=1.0;
          ind1=trunc(Int,j-1- w[i,1])
          ind2=trunc(Int,j+1+ w[i,2])
          L[j,ind1]=-(w[i,2]+1)/(w[i,1]+w[i,2]+2);
          L[j,ind2]=-(w[i,1]+1)/(w[i,1]+w[i,2]+2);
     end
return L
end







function interpA(x)
    indsz=find(x->x==0,x)
    L=eye(length(x),length(x));
    emptyLen=length(indsz);
    #w=zeros(emptyLen,2)
    w=eye(emptyLen,2)
    for i=1:emptyLen
       count1=0;
       #w{i}=count1;
       docontinue=true;
       pointer=1;
       while docontinue
            if (((i+pointer)<=length(indsz)))
                if (indsz[i+pointer]==indsz[i+pointer-1]+1)
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
                if (indsz[i-pointer+1]==indsz[i-pointer]+1)
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
       w[i,1]=count2;
       w[i,2]=count1;
       end
     for i=1:emptyLen
          j=indsz[i];
          L[j,j]=0.0;
          ind1=trunc(Int,j-1- w[i,1])
          ind2=trunc(Int,j+1+ w[i,2])
          L[j,ind1]=-(w[i,2]+1)/(w[i,1]+w[i,2]+2);
          L[j,ind2]=-(w[i,1]+1)/(w[i,1]+w[i,2]+2);
     end
return L
end


 x=[1 0 4 0 8 0 3 0 2 0 0 3 0 0 0 9]
 x=vec(x)
 interpL(x)*x


function dimensionLaplacian(Data,direction,sigma)
    for i=1:size(Data,1)
        x=(1-i):(size(Data,1)-i)
        y=(1/((1/sqrt(2*pi*sigma.^2))*exp(-(1.^2)/2*sigma.^2)))*(1/sqrt(2*pi*sigma.^2))*exp(-(x.^2)/2*sigma.^2)
        A[i,:]=y
        A[i,i]=0.0
        # prev=i-1;
        # if (prev>=1)&&(prev<=size(Data,1))
        # A[i,prev]=1.0;
        # end
        # post=i+1;
        # if (post>=1)&&(post<=size(Data,1))
        # A[i,post]=1.0;
        # end
        # prev=i-2;
        # if (prev>=1)&&(prev<=size(Data,1))
        # A[i,prev]=0.0;
        # end
        # post=i+2;
        # if (post>=1)&&(post<=size(Data,1))
        # A[i,post]=0.0;
        # end
    end
    L=speye(size(Data,1),size(Data,1))-inv(diagm(vec(sum(A,2))))*A
    L=sparse(L);
    return L
end
