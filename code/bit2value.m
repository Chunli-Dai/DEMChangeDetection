function oflag=bit2value(oflagbit);
%   oflag (right to left): bit 1, bit2, bit3, bit4, bit5, bit6
%  given a vector, convert it to a number

oflag=2^5*oflagbit(6)+2^4*oflagbit(5)+2^3*oflagbit(4)+2^2*oflagbit(3)+2^1*oflagbit(2)+2^0*oflagbit(1);

return
end
