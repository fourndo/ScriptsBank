function [ruleID,rule_lift,rules,Lbound,Ubound] = get_rules_v1(rule_file)
    

load(rule_file);
count = 1 ;
ruleID = [];
rule_lift = [];
rules = [];
Lbound = [];
Ubound{1} = [];
cc = 1;
for ii = 1 : size(Rules,1)
    
    if isnan(Rules(ii,1))==1
        
        continue
        
    end
    
    if cc ==1
        
        ruleID(count) = Rules(ii,1);
        cc = cc + 1;
        
    else
        
        rule_lift(count) = Rules(ii,1);
        

        indx = [];
        Lb = [];
        Ub = [];
        ii = ii + 1;
        while isnan(Rules(ii,2))==0 
            
            
            indx = [indx Rules(ii,2)];
            Lb =  [Lb Rules(ii,3)];
            Ub =  [Ub Rules(ii,4)];
    
            ii = ii + 1;
            
            if ii > size(Rules,1)
                
                break
                
            end
        end
        
        rules{count} = indx;
        Ubound{count} = Ub;
        Lbound{count} = Lb;
        count = count + 1;
        cc = 1;

    end
            
end

