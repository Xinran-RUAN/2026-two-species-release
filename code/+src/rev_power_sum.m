function s=rev_power_sum(a,b,g)
% (a^g-b^g)/(a-b), including the continuous definition at equal a,b.
s=ones(size(a));
for q=1:g-1,s=s.*a+b.^q;end
end
