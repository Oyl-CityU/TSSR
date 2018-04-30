function [A_final, B_final, Y_pre] = TSSR(Y, W, Sd, St, lambda_d, lambda_t, beta, n_iter, repeat_times) 
%
% Inputs:
%   Y: drug-targer association matrix.
%   W: weight matrix.
%   Sd: similarity matrix of drug.
%   St: similarity matrix of target.
%   lambda_d: parameter to control the distance between A and Sd. The default value
%   is 1.
%   lambda_t: parameter to control the distance between B and St. The default value
%   is 1.
%   beta: rate parameter of the sparse regularization. The default
%   value is 1.
%   n_iter: the number of iterations limited in RSRGM. The default value is
%   150.
%   repeat_times: repeat the algorithm with different inilization.

% Outputs:
%   A_final: sparse representation of drugs.
%   B_final: sparse representation of targets.
%   Y_pre: the predicted drug-target interaction matrix.

    if nargin < 9
        repeat_times = 5;
    end

    if nargin < 8
        n_iter = 150;
    end

    if nargin < 7
        beta = 1;
    end
    
    if nargin < 6
        lambda_t = 1;
    end
    
    if nargin < 5
        lambda_d = 1;
    end
    
    if nargin < 4
        error('Yor need input training sets and similarity matrices.');
    end
    lowest_score = inf;
    
 for i = 1: repeat_times
    fprintf(['This is the ',num2str(i), '-th repeat...'])
    fprintf('\n')

    [n,m] = size(Y); %n-drug and m-target

    A_old = Sd;
    B_old = St;

    A_old = A_old./repmat(sum(A_old,2),1,n);
    B_old = B_old./repmat(sum(B_old),m,1);

    for j  = 1: n_iter
        
        Dan = 2*(W.*Y)*(B_old'*Y') + 2*lambda_d*Sd;
        Dap = 2*(W.*(A_old*Y*B_old))*B_old'*Y' + 2*lambda_d*A_old + beta;
        aa = diag(sum(A_old./(Dap + eps),2));
        ba = sum(A_old.*(Dan./(Dap + eps)),2);
        A = A_old.*((aa*Dan + 1)./(aa*Dap + repmat(ba,1,n) + eps));
        
        Dbn = 2*(Y'*A')*(W.*Y) + 2*lambda_t*St;
        Dbp = 2*Y'*A'*(W.*(A*Y*B_old)) + 2*lambda_t*B_old + beta;
        ab = diag(sum(B_old./(Dbp + eps)));
        bb = sum(B_old.*(Dbn./(Dbp + eps)));
        B = B_old.*((Dbn*ab + 1)./(Dbp*ab + repmat(bb,m,1) + eps));

        if sum(sum(abs(A - A_old))) + sum(sum(abs(B - B_old))) < 1e-6
            break;
        else
            A_old = A;
            B_old = B;

        end
    end
    
    score = (norm(W.*(Y - A_old*Y*B_old),'fro'))^2 + lambda_d*(norm(Sd - A_old,'fro'))^2 + lambda_t*(norm(St - B_old,'fro'))^2 + beta*sum(sum(A_old)) + beta*sum(sum(B_old)); 

        if score < lowest_score

            A_final = (A + A')/2;
            B_final = (B + B')/2;
            lowest_score = score;
        end   
 end
 Y_pre = A_final*Y*B_final;
end