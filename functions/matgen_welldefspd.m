function S = matgen_welldefspd(N,den)
% generate random well-definiteness symmetric positive definite
% matrix of size N and of density den.
% den = 0.05/2; % Percentage of non-zero elements
den = den/2;
n_non_zero = round(N^2*den); % Total number of non-zero elements
non_diag_rows = randi([1 N], 1,n_non_zero); % Random string numbers for non-zeroes
non_diag_columns = randi([1 N], 1,n_non_zero); % Random column numbers for non-zeroes
non_diag_values = rand(1,n_non_zero); % Random values for non-zeroes

diag_rows = 1:1:N; % Diag rows and columns
diag_values(1, 1:N) = n_non_zero; % Diag values - any large number


S = sparse(cat(2,non_diag_rows,diag_rows),cat(2,non_diag_columns,diag_rows),cat(2,non_diag_values, diag_values));
S=(S+S');
S=S/max(max(S));
%checks for not so large matrices
% spy(S)
% condest(S)
% det(S)
end