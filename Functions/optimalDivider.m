function divider = optimalDivider( HF )
num=size(HF,1)*size(HF,2);
i=floor(sqrt(num));
divider=num/i;
while((floor(divider)~=divider)&&(floor(divider)~=divider))
divider=num/i;
i=i-1;
end

end

