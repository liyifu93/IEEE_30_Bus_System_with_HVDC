function a = insertBack(a,temp_ind,temp_a )
%------------------------------------------------------------------
% insert value "temp_a" with its original index number "temp_ind" back to
% array "a"
%------------------------------------------------------------------
% check if a is  an array
[m,n]= size(a);
if  (m-1)*(n-1)~=0
    error('insert back only valid for array')
else 
   a = reshape(a,[1,m*n]);
end


% check if there are same elements
if length(temp_ind) ~= length(unique(temp_ind) )
    error('there are same indexes during elements inserting back ')
end
%----sort insert indexes---------
[temp_ind,NewIndex] = sort(temp_ind);
temp_a = temp_a(NewIndex);
%-------------------



for j=1:length(temp_ind)
    
    if temp_ind(j) <= length(a) && temp_ind(j) > 1
        a=[a(1:temp_ind(j)-1),temp_a(j),a(temp_ind(j):end)];
    elseif temp_ind(j)==1
        a=[temp_a(j),a];
    else
        a=[a,temp_a(j)];
    end

end
a=a';
end