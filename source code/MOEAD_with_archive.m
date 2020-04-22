%MOEA/D-TCH for 3-objective and 5 objective problems (M=3 or M=5)
%T=15 for small population (i.e., N=15)
%T=20 for standard population (i.e., N=91 for 3-obj and N=210 for 5-obj)
%T=200 for large population (i.e., N=990 for 3-obj and N=1001 for 5-obj)
%T=1000 for very large population (i.e., N=5050 for 3-obj and N=5985 for 5-obj)
%-evaluation =50000 for 3-obj and 200000 for 5-obj
%The code is adapted from the PlatEMO 

clear 
clc

parfor iter=1:2

z=GLOBAL('-problem',@WFG9,'-M',3,'-N',91,'-evaluation',50000,'-D',24);
 tic

    %% Generate the weight vectors
    [W,z.N] = UniformPoint(z.N,z.M);
    T=20;
    
     [Archive0, ~]=UniformPoint(15,3);
     [Archive1, ~]=UniformPoint(91,3);
     [Archive2, ~]=UniformPoint(990,3);
     [Archive3, ~]=UniformPoint(5050,3);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);
       
    
    %% Generate random population
    Population = z.Initialization();
   
    Z = min(Population.objs,[],1);
    
    Population0=z.Initialization(15);
    Population1=z.Initialization(91); 
    Population2=z.Initialization(990);
    Population3=z.Initialization(5050);
   
    %small archive for 3 obj (j=1:15 for 5-obj)
    for j=1:15
        Population0(j)=Population(1);
    end
    
     for i=2:z.N
        g_old_archive0 = max(abs(Population0.objs-repmat(Z,15,1)).*Archive0,[],2);
        g_new_archive0 = max(repmat(abs(Population(i).obj-Z),15,1).*Archive0,[],2); 
        Population0(g_old_archive0>=g_new_archive0) = Population(i);
     end
    
    %standard archive for 3 obj (j=1:1001 for 5-obj)
    for j=1:91  
        Population1(j)=Population(1); 
    end
    
    for i=2:z.N 
        g_old_archive1 = max(abs(Population1.objs-repmat(Z,91,1)).*Archive1,[],2);
        g_new_archive1 = max(repmat(abs(Population(i).obj-Z),91,1).*Archive1,[],2);
        Population1(g_old_archive1>=g_new_archive1) = Population(i);
    end
   
    %large archive for 3 obj (j=1:1001 for 5-obj)
    for j=1:990
        Population2(j)=Population(1); 
    end
    
    for i=2:z.N
        g_old_archive2 = max(abs(Population2.objs-repmat(Z,990,1)).*Archive2,[],2);
        g_new_archive2 = max(repmat(abs(Population(i).obj-Z),990,1).*Archive2,[],2);
        Population2(g_old_archive2>=g_new_archive2) = Population(i);
    end
        
    %very large archive for 3 obj (j=1:5985 for 5-obj)
    for j=1:5050
        Population3(j)=Population(1); 
    end
    
    for i=2:z.N
        g_old_archive3 = max(abs(Population3.objs-repmat(Z,5050,1)).*Archive3,[],2);
        g_new_archive3 = max(repmat(abs(Population(i).obj-Z),5050,1).*Archive3,[],2);
        Population3(g_old_archive3>=g_new_archive3) = Population(i);
    end
        
    %% Optimization
while z.NotTermination(Population)
        % For each solution
        
        for i = 1 : z.N      
            % Choose the parents
            P = B(i,randperm(size(B,2)));

            % Generate an offspring
            Offspring = GAhalf(Population(P(1:2)));

            % Update the ideal point
            Z = min(Z,Offspring.obj);

            % Tchebycheff approach
             g_old = max(abs(Population(P).objs-repmat(Z,T,1)).*W(P,:),[],2);
             g_new = max(repmat(abs(Offspring.obj-Z),T,1).*W(P,:),[],2);    
               
            Population(P(g_old>=g_new)) = Offspring;
            
           %Archive0 (small archive)
                g_old_archive0 = max(abs(Population0.objs-repmat(Z,15,1)).*Archive0,[],2);
                g_new_archive0 = max(repmat(abs(Offspring.obj-Z),15,1).*Archive0,[],2);
                Population0(g_old_archive0>=g_new_archive0) = Offspring;  
            
            
            %Archive1 (standard archive) (for 5-obj, 91 is changed to 210)
                g_old_archive1 = max(abs(Population1.objs-repmat(Z,91,1)).*Archive1,[],2);
                g_new_archive1 = max(repmat(abs(Offspring.obj-Z),91,1).*Archive1,[],2);
                Population1(g_old_archive1>=g_new_archive1) = Offspring;
            
            %Archive2 (large archive) (for 5-obj, 990 is changed to 1001)
                g_old_archive2 = max(abs(Population2.objs-repmat(Z,990,1)).*Archive2,[],2);
                g_new_archive2 = max(repmat(abs(Offspring.obj-Z),990,1).*Archive2,[],2);
                Population2(g_old_archive2>=g_new_archive2) = Offspring;
                             
            %Archive3 (very large archive) (for 5-obj, 5050 is changed to 5985)
                g_old_archive3 = max(abs(Population3.objs-repmat(Z,5050,1)).*Archive3,[],2);
                g_new_archive3 = max(repmat(abs(Offspring.obj-Z),5050,1).*Archive3,[],2);
                Population3(g_old_archive3>=g_new_archive3) = Offspring;                  
        end
end
parsave(sprintf('MOEADSmall_%s_M%d_D%d_%d.mat',class(z.problem),z.M, z.D,iter),Population,Population1,Population2,Population3, Population0);

 end