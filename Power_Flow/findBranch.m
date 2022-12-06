function [ sorted_branch, B_half ] = findBranch( Y )
%FINDBRANCH Summary of this function goes here
%   return branch pair in rows
%    [1 7;
%     4 8; 
%     6 8] 

                          
      [row,col]=size(Y);
     

     % get branch index in the Ybus matrix 
     branch=[];
     for n=2:col
         for k=1:(n-1)
            if nonzeros(abs(Y(k,n)))   
                branch=[branch;[k,n]];
            end
         end
     end 

sorted_branch =branch;  
%%
%           [from   to]
% branch =  [from   to]
%           [from   to]
%%
       for j=1:(length(sorted_branch)-1)
          for k=1:(length(sorted_branch)-j)
              if sorted_branch(k,1) > sorted_branch(k+1,1)
                  med=sorted_branch(k+1,:);
                  sorted_branch(k+1,:)=sorted_branch(k,:);
                  sorted_branch(k,:)=med;     
              end
          end
       end    
branch_switched = [sorted_branch(:,2),sorted_branch(:,1)];
sorted_branch =[sorted_branch;branch_switched];
       
       
%%   paird branch below    
%        for j=1:(length(sorted_branch)-1)
%           for k=1:(length(sorted_branch)-j)
%               if sorted_branch(k,2) > sorted_branch(k+1,2) && sorted_branch(k,1) == sorted_branch(k+1,1)
%                   med=sorted_branch(k+1,:);
%                   sorted_branch(k+1,:)=sorted_branch(k,:);
%                   sorted_branch(k,:)=med;     
%               end
%           end
%        end     
     

end

