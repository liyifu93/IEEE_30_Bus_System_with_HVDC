% This program will omit Q in volt-controlled bus

% 1. Parameters
% Branch=Branch30;
    % 1         2           3         4         5          6
 %  fromBus    toBus        R         X         B        ratio      phaseshtif
Branch=[
    1.0000    2.0000    0.0192    0.0575    0.0528         0         0
    1.0000    3.0000    0.0452    0.1652    0.0408         0         0
    2.0000    4.0000    0.0570    0.1737    0.0368         0         0
    3.0000    4.0000    0.0132    0.0379    0.0084         0         0
    2.0000    5.0000    0.0472    0.1983    0.0418         0         0
    2.0000    6.0000    0.0581    0.1763    0.0374         0         0
    4.0000    6.0000    0.0119    0.0414    0.0090         0         0
    5.0000    7.0000    0.0460    0.1160    0.0204         0         0
    6.0000    7.0000    0.0267    0.0820    0.0170         0         0
    6.0000    8.0000    0.0120    0.0420    0.0090         0         0
    6.0000    9.0000         0    0.2080         0    0.9780         0
    6.0000   10.0000         0    0.5560         0    0.9690         0
    9.0000   11.0000         0    0.2080         0         0         0
    9.0000   10.0000         0    0.1100         0         0         0
    4.0000   12.0000         0    0.2560         0    0.9320         0
   12.0000   13.0000         0    0.1400         0         0         0
   12.0000   14.0000    0.1231    0.2559         0         0         0
   12.0000   15.0000    0.0662    0.1304         0         0         0
   12.0000   16.0000    0.0945    0.1987         0         0         0
   14.0000   15.0000    0.2210    0.1997         0         0         0
   16.0000   17.0000    0.0524    0.1923         0         0         0
   15.0000   18.0000    0.1073    0.2185         0         0         0
   18.0000   19.0000    0.0639    0.1292         0         0         0
   19.0000   20.0000    0.0340    0.0680         0         0         0
   10.0000   20.0000    0.0936    0.2090         0         0         0
   10.0000   17.0000    0.0324    0.0845         0         0         0
   10.0000   21.0000    0.0348    0.0749         0         0         0
   10.0000   22.0000    0.0727    0.1499         0         0         0
   21.0000   22.0000    0.0116    0.0236         0         0         0
   15.0000   23.0000    0.1000    0.2020         0         0         0
   22.0000   24.0000    0.1150    0.1790         0         0         0
   23.0000   24.0000    0.1320    0.2700         0         0         0
   24.0000   25.0000    0.1885    0.3292         0         0         0
   25.0000   26.0000    0.2544    0.3800         0         0         0
   25.0000   27.0000    0.1093    0.2087         0         0         0
   28.0000   27.0000         0    0.3960         0    0.9680         0
   27.0000   29.0000    0.2198    0.4153         0         0         0
   27.0000   30.0000    0.3202    0.6027         0         0         0
   29.0000   30.0000    0.2399    0.4533         0         0         0
    8.0000   28.0000    0.0636    0.2000    0.0428         0         0
    6.0000   28.0000    0.0169    0.0599    0.0130         0         0];

% Bus=Bus30;
 %Bus Type
 % 0 - Unregulated (load, PQ)
 % 1 - Hold MVAR generation within voltage limits, (PQ)
 % 2 - Hold voltage within VAR limits (gen, PV)
 % 3 - Hold voltage and angle (swing, V-Theta) (must always have one)
     %1        2           3        4         5        6         7         8           9          10        11         12
   %Bus   V_desired(kV)  Type      v(pu)    angle    P-load    Q-load     P-Gen     Q-Gen   switch_shunt  shunt_G   shunt_B       S_dc   busWithHarm
Bus=[
    1.0000  132.0000    3.0000    1.0600   98.4320         0         0  260.0000  -16.7900         0         0         0            0    0
    2.0000  132.0000    2.0000    1.0431   93.0800   21.7000   12.7000   40.0000   50.0000         0         0         0            0    0
    3.0000  132.0000         0    1.0207   90.9000    2.4000    1.2000         0         0         0         0         0            0    0
    4.0000  132.0000         0    1.0118   89.1480    7.6000    1.6000         0         0         0         0         0            0    0
    5.0000  132.0000    2.0000    1.0100   84.2660   94.2000   19.0000         0   36.8000         0         0         0            0    0   
    6.0000  132.0000         0    1.0103   87.3670         0         0         0         0         0         0         0            0    0
    7.0000  132.0000         0    1.0024   85.5660   22.8000   10.9000         0         0         0         0         0            0    0
    8.0000  132.0000    2.0000    1.0100   86.6180   30.0000   30.0000         0   37.1000         0         0         0            0    0
    9.0000   33.0000         0    1.0509   84.3230         0         0         0         0         0         0         0            0    0
   10.0000   33.0000         0    1.0451   82.7320    5.8000    2.0000         0         0         0         0         0.19         0    0
   11.0000   33.0000    2.0000    1.0820   84.3230         0         0         0   16.2000         0         0         0            0    0
   12.0000   33.0000         0    1.0571   83.4880   11.2000    7.5000         0         0         0         0         0            0    0
   13.0000   33.0000    2.0000    1.0710   83.4880         0         0         0   10.6000         0         0         0            0    0
   14.0000   33.0000         0    1.0423   82.5960    6.2000    1.6000         0         0         0         0         0            0    0
   15.0000   33.0000         0    1.0377   82.5040    8.2000    2.5000         0         0         0         0         0            0    0
   16.0000   33.0000         0    1.0444   82.9050    3.5000    1.8000         0         0         0         0         0            0    0
   17.0000   33.0000         0    1.0399   82.5700    9.0000    5.8000         0         0         0         0         0            0    0
   18.0000   33.0000         0    1.0282   81.8900    3.2000    0.9000         0         0         0         0         0            0    0
   19.0000   33.0000         0    1.0257   81.7160    9.5000    3.4000         0         0         0         0         0            0    0
   20.0000   33.0000         0    1.0297   81.9130    2.2000    0.7000         0         0         0         0         0            0    0
   21.0000   33.0000         0    1.0327   82.2890   17.5000   11.2000         0         0         0         0         0            0    0
   22.0000   33.0000         0    1.0333   82.3030         0         0         0         0         0         0         0            0    0
   23.0000   33.0000         0    1.0272   82.1140       3.2       1.6         0         0         0         0         0            0    0
   24.0000   33.0000         0    1.0216   81.9370    8.7000    6.7000         0         0         0         0         0.043        0    0
   25.0000   33.0000         0    1.0173   82.3650         0         0         0         0         0         0         0            0    0
   26.0000   33.0000         0    0.9997   81.9450       3.5       2.3         0         0         0         0         0            0    0
   27.0000   33.0000         0    1.0232   82.8890         0         0         0         0         0         0         0            0    0
   28.0000  132.0000         0    1.0068   86.7430         0         0         0         0         0         0         0            0    0
   29.0000   33.0000         0    1.0034   81.6590    2.4000    0.9000         0         0         0         0         0            0    0
   30.0000   33.0000         0    0.9919   80.7760   10.6000    1.9000         0         0         0         0         0            0    0];
% Bus=[
%     1.0000  132.0000    3.0000    1.35542   132.4393         0         0  260.0000  -16.7900         0         0         0            0    0
%     2.0000  132.0000    2.0000    1.16728   98.7264   21.7000   12.7000   40.0000   50.0000         0         0         0            0    0
%     3.0000  132.0000         0    1.03545   90.900    2.4000    1.2000         0         0         0         0         0            0    0
%     4.0000  132.0000         0    1.02971   89.1480    7.6000    1.6000         0         0         0         0         0            0    1
%     5.0000  132.0000    2.0000    1.11216   83.9903  94.2000   19.0000         0   36.8000         0         0         0            0    0   
%     6.0000  132.0000         0    1.0103   87.3670         0         0         0         0         0         0         0            0    0
%     7.0000  132.0000         0    1.0024   85.5660   22.8000   10.9000         0         0         0         0         0            0    0
%     8.0000  132.0000    2.0000    1.0100   86.6180   30.0000   30.0000         0   37.1000         0         0         0            0    0
%     9.0000   33.0000         0    1.0509   84.3230         0         0         0         0         0         0         0            0    0
%    10.0000   33.0000         0    1.0451   82.7320    5.8000    2.0000         0         0         0         0         0.19         0    0
%    11.0000   33.0000    2.0000    1.0820   84.3230         0         0         0   16.2000         0         0         0            0    0
%    12.0000   33.0000         0    1.0571   83.4880   11.2000    7.5000         0         0         0         0         0            0    1
%    13.0000   33.0000    2.0000    1.0710   83.4880         0         0         0   10.6000         0         0         0            0    0
%    14.0000   33.0000         0    1.0423   82.5960    6.2000    1.6000         0         0         0         0         0            0    0
%    15.0000   33.0000         0    1.0377   82.5040    8.2000    2.5000         0         0         0         0         0            0    0
%    16.0000   33.0000         0    1.0444   82.9050    3.5000    1.8000         0         0         0         0         0            0    0
%    17.0000   33.0000         0    1.0399   82.5700    9.0000    5.8000         0         0         0         0         0            0    0
%    18.0000   33.0000         0    1.0282   81.8900    3.2000    0.9000         0         0         0         0         0            0    0
%    19.0000   33.0000         0    1.0257   81.7160    9.5000    3.4000         0         0         0         0         0            0    0
%    20.0000   33.0000         0    1.0297   81.9130    2.2000    0.7000         0         0         0         0         0            0    0
%    21.0000   33.0000         0    1.0327   82.2890   17.5000   11.2000         0         0         0         0         0            0    0
%    22.0000   33.0000         0    1.0333   82.3030         0         0         0         0         0         0         0            0    0
%    23.0000   33.0000         0    1.0272   82.1140       3.2       1.6         0         0         0         0         0            0    0
%    24.0000   33.0000         0    1.0216   81.9370    8.7000    6.7000         0         0         0         0         0.043        0    0
%    25.0000   33.0000         0    1.0173   82.3650         0         0         0         0         0         0         0            0    0
%    26.0000   33.0000         0    0.9997   81.9450       3.5       2.3         0         0         0         0         0            0    0
%    27.0000   33.0000         0    1.0232   82.8890         0         0         0         0         0         0         0            0    0
%    28.0000  132.0000         0    1.0068   86.7430         0         0         0         0         0         0         0            0    0
%    29.0000   33.0000         0    1.0034   81.6590    2.4000    0.9000         0         0         0         0         0            0    0
%    30.0000   33.0000         0    0.9919   80.7760   10.6000    1.9000         0         0         0         0         0            0    0];
% set initial stack angle to be zero
Bus(:,5)=Bus(:,5)-Bus(1,5);


[numOfBus,~]=size(Bus);
[numOfBranch,~]=size(Branch);
Y=zeros(numOfBus,numOfBus);
B_half=Y;
tap = ones(numOfBranch,1) ;
tap_index = find(Branch(:, 6));   % find index of tap
tap(tap_index) = Branch(tap_index, 6).* exp(1i*pi/180 * Branch(tap_index, 7));



% 2. YBus30
% 1) calculate diagonal element in Ybus matrix
for j=1:numOfBus
    for row=1:numOfBranch
        if Branch(row,1)==j
            Y(j,j)=Y(j,j)+(  1/( Branch(row,3)+1i*Branch(row,4) )+1i/2*Branch(row,5)  )/(tap(row) * conj(tap(row)));
        end
        
        if Branch(row,2)==j
            Y(j,j)=Y(j,j)+(  1/( Branch(row,3)+1i*Branch(row,4) )+1i/2*Branch(row,5)  );
        end
    end
end
   
% apply shunt capactitor
for j=1:numOfBus
    Y(j,j) = Y(j,j) + 1i * Bus(j,12);  % Bus(:,12) is shunt capacitor in Pu;
end  
   
% 2) calculate non-diagonal element in Yadmittance matrix
% Y(row,col)=  - y(row,col)= - 1/(R+jX)  
for row=1:numOfBranch
    Y(Branch(row,1),Branch(row,2))= - ( 1/( Branch(row,3)+1i*Branch(row,4) )+0i/2*Branch(row,5)  )/conj(tap(row));
    Y(Branch(row,2),Branch(row,1))= - ( 1/( Branch(row,3)+1i*Branch(row,4) )+0i/2*Branch(row,5)  )/  tap(row) ;
    
    B_half(Branch(row,1),Branch(row,2)) = ( 1i/2*Branch(row,5)+1i*Bus(j,12) )/conj(tap(row))  ;
    B_half(Branch(row,2),Branch(row,1)) = ( 1i/2*Branch(row,5)+1i*Bus(j,12) ) /  tap(row);
end


Ymag=abs(Y);
Yangle = angle(Y);
[numOfBus,col]=size(Y);
% x=[delta1; delta2 ... deltaN; v1; v2;... vN]   % x = [angle ; voltage]
x=zeros(numOfBus*2,1);
% y=[P1;P2;...PN;Q1;Q2;...QN]   % y = [P ; Q]
y=x;
PV=[];
PQ=[];

% 3. PF
% 1) Find bus type and initial 'x' 'y'   x=[θ ; V], y=[P ; Q]
%  the lengths of both 'x' and 'y'are unchanged in the entire process
for j=1:numOfBus
    if Bus(j,3)==3   % swing
        x(j)=0;  % x(j,1)=Bus(j,5)
        x(j+numOfBus)=Bus(j,4);
        slack=j;
    end
    
     if Bus(j,3)==2  % PV generator
        y(j,1)=Bus(j,8)/100-Bus(j,6)/100;
        y(j+numOfBus,1)=Bus(j,9)/100-Bus(j,7)/100;
        x(j+numOfBus,1)=Bus(j,4);
        PV=[PV,j];
        x(j)=0; % initial arbitrary phase angle in PV bus:  x(j)=Bus(j,5);
     end
    
     if Bus(j,3)==0  % PQ load
        y(j,1)=  - Bus(j,6)/100;
        y(j+numOfBus,1)= - Bus(j,7)/100;
        PQ=[PQ,j];
        x(j)=0; % initial phase angle in PQ bus:  x(j)=Bus(j,5)
        x(j+numOfBus)=Bus(j,4); % initial arbitrary voltage in PQ bus  
     end
end

%%
% 3)  -----  START  AC Solver ----
% Initial Variables
Jac1=zeros(numOfBus);
Jac2=Jac1;
Jac3=Jac1;
Jac4=Jac1;
% initial calculated P Q
P=zeros(numOfBus,1); 
Q=P;
% Separate Vangle and Vmag from x
    Vangle=x(1:numOfBus);
    Vmag=x(1+numOfBus:end);
% Separate P Q from y
    Pref=y(1:numOfBus);
    Qref=y(1+numOfBus:end);
 converge=false;   
 interateion=0;
 %
while ~converge && interateion<100
    
    interateion=interateion+1;    
    
%    generate Jacobian
    for k=1:numOfBus
       for n=1:numOfBus
           if k~=n
               Jac1(k,n) =  Vmag(k) *Ymag(k,n) *Vmag(n) *sin( Vangle(k) -Vangle(n) -Yangle(k,n) );
               Jac2(k,n) =  Vmag(k) *Ymag(k,n)          *cos( Vangle(k) -Vangle(n) -Yangle(k,n) );
               Jac3(k,n) = -Vmag(k) *Ymag(k,n) *Vmag(n) *cos( Vangle(k) -Vangle(n) -Yangle(k,n) );
               Jac4(k,n) =  Vmag(k) *Ymag(k,n)          *sin( Vangle(k) -Vangle(n) -Yangle(k,n) ); 
               
           else% k = n
               
               Jac1(k,n) = -Vmag(k) *( sum(    Ymag(k,:)' .*Vmag .*sin( Vangle(k) -Vangle -Yangle(k,:)' ) )- Ymag(k,n) *Vmag(k) *sin( Vangle(k) -Vangle(n) -Yangle(k,n) )    );
               Jac2(k,n) =  Vmag(k) *Ymag(k,k) *cos( Yangle(k,k)) + sum(   Ymag(k,:)' .*Vmag .*cos( Vangle(k) -Vangle -Yangle(k,:)' )   );
               Jac3(k,n) =  Vmag(k) *( sum(    Ymag(k,:)' .*Vmag .*cos( Vangle(k) -Vangle -Yangle(k,:)' ) )- Ymag(k,n) *Vmag(k) *cos( Vangle(k) -Vangle(n) -Yangle(k,n) )    );
               Jac4(k,n) = -Vmag(k) *Ymag(k,k) *sin( Yangle(k,k)) + sum(   Ymag(k,:)' .*Vmag .*sin( Vangle(k) -Vangle -Yangle(k,:)' )   );
           end
       end
    end
    
% Remove slack bus in Jac1 matrix 
Jac1(slack,:)=[];
Jac1(:,slack)=[];
Jac2(slack,:)=[];
Jac2(:,slack)=[];
Jac3(slack,:)=[];
Jac3(:,slack)=[];
Jac4(slack,:)=[];
Jac4(:,slack)=[];
% Remove Qref row of V-controlled bus in Jacobian matrix
Jac3(PV-1,:)=[];
Jac4(PV-1,:)=[];
% Remove Vmag column of V-controlled bus in Jacobian matrix 
Jac2(:,PV-1)=[];
Jac4(:,PV-1)=[];
%              [ J1 | J2 ]
% Jacobian  =  ----------
%              [ J3 | J4 ]
Jacob=[Jac1,Jac2;Jac3,Jac4];


% Calculate  P[x(i)] and Q[x(i)]
    for k=1:(numOfBus)  
             P(k)= Vmag(k) *sum(   Ymag(k,:)' .*Vmag .*cos( Vangle(k)-Vangle-Yangle(k,:)' )      );
             Q(k)= Vmag(k) *sum(   Ymag(k,:)' .*Vmag .*sin( Vangle(k)-Vangle-Yangle(k,:)' )      ) ;
    end
    


% 4) ---During interation, slack bus are omitted---     
%      Remove slack in calculated  P[x(i)] and Q[x(i)] 
        P(slack) = [];
        Q([slack,PV]) = [];
%      Remove slack in Vangle and Vmag 
%      and V in Vmag for v-controlled bus
       temp_Vangle = Vangle(slack);  Vangle(slack) = [];
       temp_Vmag = Vmag([slack,PV]); Vmag([slack,PV]) = [];
%      Remove slack in Pref and Qref   
%      and Q in Qref for V-controlled bus 
        Pref(slack)=[];
        Qref([slack,PV])=[];
        
% 5) ---START interation ---
dy=[Pref;Qref]-[P;Q];
dx=(Jacob)\dy;
 % Update Vangle for PV and PQ, Update Vmag for PQ
         Median = [Vangle;Vmag];
         Median = Median+dx;     
         Vangle = Median(1: length(Vangle));
         Vmag = Median(length(Vangle)+1:end);

if norm(dx,inf)<0.01
    converge=1;
end

        % x length unchanged
        % Update PV PQ Bus angle to new angle matrix
        Vangle = insertBack(Vangle,slack,temp_Vangle);
        
        % Update PQ Bus angle to new angle matrix
        Vmag = insertBack(Vmag,[slack,PV],temp_Vmag);
        % update 'x'
        x=[Vangle; Vmag];
        
            % Separate Vangle and Vmag from x again
            Vangle=x(1:numOfBus);
            Vmag=x(1+numOfBus:end);
            % Separate P Q from y again
            Pref=y(1:numOfBus);
            Qref=y(1+numOfBus:end);
            % Wrap angle
           Vangle=angle(exp(1i*Vangle));        
            % Update Qref for voltage controlled bus
    for j=1:length(PV)
        k=PV(j);
        Qv= Vmag(k) *sum(    Ymag(k,:)' .*Vmag .*sin( Vangle(k,1) -Vangle -Yangle(k,:)' ));
        y(k+numOfBus)=Qv;
        %disp({'Qv' k ,'update to y',k+numOfBus} )
    end
    
    
end

% 6) After converge 
% insert calculated P and Q for slack bus 
    y(slack)= Vmag(slack) *sum(   Ymag(slack,:)' .*Vmag .*cos( Vangle(slack)-Vangle-Yangle(slack,:)' )      );
    y(slack+numOfBus)= Vmag(slack) *sum(   Ymag(slack,:)' .*Vmag .*sin( Vangle(slack)-Vangle-Yangle(slack,:)' )      );
    y(slack)=y(slack)+Bus(slack,6);  % y(slack)=y(slack)+Bus(slack,15);
    y(slack+numOfBus)=y(slack+numOfBus)+Bus(slack,7);  % y(slack+numOfBus)=y(slack+numOfBus)+Bus(slack,16);
 
% calculate brach power flow
branch = findBranch(Y);
[rowB,colB] = size(branch) ;
S_branch = zeros( rowB,1 ) ;
for k=1:rowB
    theta = Vangle( branch(k,2) ) - Vangle( branch(k,1) ) ;
    S_branch(k,1) =  Vmag( branch(k,1) ) *conj( ( Vmag( branch(k,1) ) - Vmag( branch(k,2) )*exp(theta*1i)) * -Y( branch(k,1), branch(k,2) ) + Vmag( branch(k,1) )*B_half(branch(k,1), branch(k,2))  );
end 

disp( '    FromBUS    ToBUS    P(WM)    Q(MWAr) ' )
disp( [branch(1:rowB/2,:), real(S_branch(1:rowB/2))*100, imag(S_branch(1:rowB/2))*100 ] )
disp( '    FromBUS    ToBUS    P(WM)    Q(MWAr) ' )
disp( [branch(rowB/2+1:end,:), real(S_branch(rowB/2+1:end))*100, imag(S_branch(rowB/2+1:end))*100 ] )


% 7) Show final  Bus(node) result
Pgen=zeros(numOfBus,1);
Qgen=Pgen;
P_Load=Pgen;
Q_Load=Pgen;


Bus(4,6)=Bus(4,6)-real(Bus(4,13));
Bus(4,7)=Bus(4,7)-imag(Bus(4,13));
Bus(12,6)=Bus(12,6)-real(Bus(12,13));
Bus(12,7)=Bus(12,7)-imag(Bus(12,13));
for j=1:numOfBus
    if j==slack % Slack
        Pgen(j)=y(j)*100;
        Qgen(j)=y(j+numOfBus)*100;
        P_Load(j)=Bus(j,6);
        Q_Load(j)=Bus(j,7);
    elseif find(j==PV) % PV
        Pgen(j)=y(j)*100+Bus(j,6) ;
        Qgen(j)=y(j+numOfBus)*100+Bus(j,7); 
        P_Load(j)=Bus(j,6);
        Q_Load(j)=Bus(j,7);
    else % PQ
        Pgen(j)=Bus(j,8);
        Qgen(j)=Bus(j,9);
        P_Load(j)=Bus(j,6);
        Q_Load(j)=Bus(j,7);
    end

end

disp( '     BUS       Tpye     V(pu)   Angle(deg)  P-load    Q-load     P-gen   Q-gen' )
disp(  ( [(1:numOfBus)' ,(Bus(:,3)),Vmag ,Vangle*180/pi,P_Load,Q_Load,Pgen,Qgen ] )  )