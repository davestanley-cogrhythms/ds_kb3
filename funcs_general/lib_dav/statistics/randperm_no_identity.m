
function ind = randperm_no_identity(N)
    % function ind = randperm_no_identity(N)
    % Works just like matlab's randperm(N) command, returning a random
    % permutation of the values 1:Ntrials. However, it will not return any
    % identity values.
    % For example, randperm(4) could return [2 1 4 3] or [2 1 3 4]
    % But randperm_no_identity(4) could only return [2 1 4 3] because the
    % array [2 1 3 4] contains an identity (the 3,4 part).
    % Uses a montecarlo type approach on the original matlab randperm
    % command
    % To do : modify to work with randperm(N,k)
    % Author: David Stanley, Boston University, 2015
%     tic
    ind = randperm(N);
    ind_id = 1:N;
    bads = find((ind - ind_id) == 0);       % Find the indicies of identities
    
    while ~isempty(bads)    % If there are any identities in the random permutation
        shuff_ind = floor(unifrnd(1,N+1-0.1,1,1));          % Find some other index at random
        ind([bads(1) shuff_ind]) = ind([shuff_ind bads(1)]);    % Swap the bad index with the other index, and hope this removes the identity; if not, we will repeat
        bads = find((ind - ind_id) == 0);           % Recalculate the number of identities remaining.
    end
%     ind;
%     toc
    

end