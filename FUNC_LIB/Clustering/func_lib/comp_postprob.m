function [Pprob,top_rule,OUTid,OUThit] = comp_postprob(d,rules,Ub,Lb,Wp,Wm,out_col,var1,var2)

top_rule = [];
ndata = size( d ,1 );
nrules = size(Ub,2);
flag = 0;

% FOG
OUTvar = d(:,out_col); %OUTvar(OUTvar==-99999) = flag;
% Targeting
% OUTvar = data(:,end); OUTvar(OUTvar==1) = 0;

OUTid = unique(OUTvar); OUTid = OUTid(2:end);
OUTpop = zeros(length(OUTid),1);

% Compute the number of training sites per FOG
for ii = 1 : length(OUTid);

    OUTpop(ii) = sum(OUTvar == OUTid(ii));
    
end

% Training sites
outcome = OUTvar;

nts = sum(outcome~=0);
nnts = ndata - nts;



pre_prob = nts / ndata;
pre_logit = ones(ndata,1) * log( pre_prob / (1 - pre_prob));


%% Compute post-probability score

% Add weights for all stations and rank
post_logit = pre_logit ;

OUThit = zeros( nrules , length(OUTid) );
coverage = spalloc( ndata , nrules, nrules*nts);

switch var1
    
    case 'WofE'
     
        for ii = 1 : size(rules,2)
       
            if rules{ii}(1,2) == 1 

               % Assign score using computed weights     
                logic =  d(:,rules{ii}(1,1)) >= Lb{ii}(1) ;

            elseif rules{ii}(1,2) == -1
                 
                logic =  d(:,rules{ii}(1,1)) <= Ub{ii}(1) ;

            elseif  rules{ii}(1,2) == 0 
                
                logic = zeros(size(d,1),1);
                for jj = 1 : length(Lb{ii})

                logic =  logic | d(:,rules{ii}(1,1)) == Lb{ii}(jj) ;
                
                end
                
                logic = logic>=1;

            end
            
            post_logit(logic) = post_logit(logic) + Wp(ii);
            post_logit(logic==0) = post_logit(logic==0) + Wm(ii);
            
        end
        
        Pprob = exp(post_logit) ./ (1+exp(post_logit));
        

        %% Find only the rules that matter
    otherwise
         
        
        
        purity = zeros(size(rules,2),1);

        fprintf('\nComputing coverage of given rules\n');
        % Print progress
        progress=-1;
        tic

        % First compute how many Training cells are detected by each rule        
        for jj = 1 : size(rules,2)

            logic = ones(ndata,1);

            for ii = 1 : size(rules{jj},2)

                % Check if continuous variable
                if rules{jj}(2,ii) == 1

                logic =  logic &...
                    ( d(:,rules{jj}(1,ii)) <= Ub{jj}(ii) ) &...
                    ( d(:,rules{jj}(1,ii)) >= Lb{jj}(ii) );

                elseif rules{jj}(2,ii) == 0

                logic =  logic &...
                    ( d(:,rules{jj}(1,ii)) == Lb{jj}(ii) );  

                end

            end

            % Save purity of each rule
            purity(jj) = sum(outcome(logic)~=0) / sum(logic);

            temp = unique(OUTvar(logic));
            temp = temp(2:end);

            coverage(logic,jj) = 1;

            for kk = 1 : length(temp)

                OUThit( jj , OUTid == temp(kk) ) = 1;

            end

             % Print progress
            d_iter = floor(jj/size(rules,2)*100);
            if  d_iter > progress

                fprintf('Processed %i pct of data in %8.5f sec\n',d_iter,toc)
                progress = d_iter;

            end


        end

        switch var2
            case 'optimize'
                % Rank the rules

                % Ranked by number of unique FOG
                hitters = sum(OUThit,2);                
                [~,cc] = sort(hitters,1,'descend');

                % Ranked by number of FOG training points
                %         hitters = sum(lifters,1);                
                %         [~,cc] = sort(hitters,1,'descend');


                % Take the top rule with max FOG detected, then find the next rule 
                % that will increase knowledge
                knowledge = zeros(1,size(OUThit,2));
                top_rule = [];
                count = 0;

        %         Loop until you find at least one point per FOG#
                while count < sum(sum(OUThit,1)~=0)

                    leader = -99;

                    for jj = 1 : size(OUThit,1)

                        temp = knowledge.*0;
                        temp(OUThit(jj,:)==1 & knowledge==1) = 0;
                        temp(OUThit(jj,:)==1 & knowledge==0) = 1;
        %                 temp( temp~=1) = 0;

                        if sum(temp) > leader

                            leader = sum(temp);
                            candidate = jj;
                            index = jj;

                        % If more than one rule with the same ranking...
                        elseif sum(temp) == leader

                            candidate = [candidate;jj];

                        end


                    end

                    % ...keep the rule with highest purity?
                    criteria = purity;
                    % ...keep the rule with highest coverage?
%                     criteria = sum(coverage,1);
                    
                    
                    if length(candidate)>1

                        filter = criteria(candidate(1));
                        index = candidate(1);
                        for ii = 1 : length(candidate)

                            
                            if criteria(candidate(ii)) > filter

                                filter = criteria(candidate(ii));
                                index = candidate(ii);

                            end

                        end

                        

                    end

                    top_rule = [top_rule index];
                    knowledge = knowledge + OUThit(index,:);
                    knowledge(knowledge~=0) = 1;
                    count = sum(knowledge);

                end
               
        
            otherwise
                
                top_rule = 1 : size(rules,2);
                             
        end
        
        r_idx =  top_rule;
        %% COMPUTE POST-PROB FOR ALL SELECT RULES
        
        % Print progress
        progress=-1;
        tic
        fprintf('\n!!!Computing Post-Probility!!!\n') 
        for rr = 1 : length(r_idx)
         %%    
            jj = r_idx(rr);
            logic = ones(ndata,1);

            for ii = 1 : size(rules{jj},2)    
                % Check if continuous variable
                if rules{jj}(2,ii) == 1
       
                    logic =  logic &...
                        ( d(:,rules{jj}(1,ii)) <= Ub{jj}(ii) ) &...
                        ( d(:,rules{jj}(1,ii)) >= Lb{jj}(ii) );
                
                elseif rules{jj}(2,ii) == 0
                    
                    logic =  logic &...
                        ( d(:,rules{jj}(1,ii)) == Lb{jj}(ii) );  
                
                end
                
            end

            post_logit(logic) = post_logit(logic) + Wp(jj);
            post_logit(logic==0) = post_logit(logic==0) + Wm(jj);
            

             % Print progress
            d_iter = floor(rr/length(r_idx)*100);
            if  d_iter > progress

                fprintf('Processed %i pct of data in %8.5f sec\n',d_iter,toc)
                progress = d_iter;

            end
        
        end
        
        Pprob = post_logit;

end

