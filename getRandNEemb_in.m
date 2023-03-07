function [emb]=getRandNEemb_in(A,paras, savefilename,TableNode, worktype)
    % A sample run on the BlogCatalog Dataset  
    % load data
    % load('BlogCatalog');
    if ~exist('A','var') 
        warning('test');
        nn =1000; 
        A= rand(nn); A = A + A';  savefilename = 't11111111111.emb.txt'; 
        paras.d       =128;
        paras.Ortho   = 1;
        paras.q       = 2 ;
        paras.weights = [1,0.1,0.001];  
        paras.seed    = 0;
        TableNode = table('ID'+string(1:nn)');
        TableNode =[]
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    d       = 128;
    Ortho   = 1;
    seed    = 0;
    % % q       = paras.q;
    % % weights = paras.weights; 
    %
    N = length(A);
    % Common parameters
    % d = 128;
    % Ortho = 1;
    % seed = 0;
    switch worktype
        % worktype = 'classification';  'reconstruction';  
        case 'reconstruction'
            % embedding for adjacency matrix for reconstruction 
            q = 2;
            weights = [1,0.1,0.001];
            U_list = RandNE_Projection(A,q,d,Ortho,seed);
            U = RandNE_Combine(U_list,weights);
            % prec = Precision_Np(A,sparse(N,N),U,U,1e6);
            % semilogx(1:1e6,prec);
        case 'classification'
            % % embedding for transition matrix for classification     
            q = 3;
            weights = [1,1e2,1e4,1e5];
            A_tran = spdiags(1 ./ sum(A,2),0,N,N) * A;
            U_list = RandNE_Projection(A_tran,q,d,Ortho,seed);
            U = RandNE_Combine(U_list,weights);
            % % normalizing
            U = spdiags(1 ./ sqrt(sum(U .* U,2)),0,N,N) * U;
            % % Some Classification method, such as SVM in http://leitang.net/social_dimension.html
        otherwise
            error('No defintion')
    end
    [n_node, n_feature] = size( U ); 
    % Delimiter is space
    if ~isempty(TableNode)
        IDstr = TableNode{:,1}; 
    % % %     tbl_emb = table(IDstr, U);
    % % %     writetable(tbl_emb, savefilename, 'Delimiter', ' ' ,'WriteVariableNames',false );
        if ~isempty( savefilename )
            fileID = fopen(savefilename,'w+');
            fprintf(fileID,'%d %d\n', n_node, n_feature ); 
            fmt =strjoin(['%s',repmat({'%f'},1, n_feature),'\n'], ' ' ); 
            for ii=1:n_node
                fprintf(fileID,fmt, IDstr{ii}, U(ii,:)); 
            end
            fclose(fileID);
        end

    else
        U1=(0:N-1)';
        emb=[U1 U];
        if ~isempty( savefilename )
            fileID = fopen(savefilename,'w+');
            fprintf(fileID,'%d %d\n', n_node, n_feature ); 
            fmt =strjoin(['%d',repmat({'%f'},1, n_feature),'\n'], ' ' );   
            fprintf(fileID,fmt, emb'); 
            fclose(fileID);
        end
    end
end
