function rd=rela_diff(a,b)
%calculate the relative difference between two variables.

rd=(b-a)./((abs(a)+abs(b))/2);

return
end