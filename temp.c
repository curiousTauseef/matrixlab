function C = matmul(A, B, blocksize)

   %------------------------------------------------------------------
   % Compute C = A*B, where A is a n x n array. Do this by cutting
   % the arrays into subblocks of size blocksize x blocksize. This
   % allows enhancing cache re-use. Warning: this function is only
   % for when mod(n, blocksize) = 0.
   %------------------------------------------------------------------

   [m, n] = size(A);
   C = zeros(n, n);

   for ii = 1:blocksize:n
      % --------------------------------------------------
      % Indices in overall A array of current block row:
      % --------------------------------------------------
      I = ii:ii+blocksize-1;

      for jj = 1:blocksize:n
        % ----------------------------------------------------
        % Indices in overall A array of current block column:
        % ----------------------------------------------------
        J = jj:jj+blocksize-1;

        %--------------------------------------------------------------
        % Copy over A's submatrix of size (ii/blocksize, jj/blocksize)
        % to the array T, transposing it as it is copied.
        %--------------------------------------------------------------
        T = A(I, J)';

        %----------------------------------------
        % For each block column of C and B ...
        %----------------------------------------
        for kk = 1:blocksize:n

            %-----------------------------------------------------------
            % Define K as the array of indices in C and B corresponding
            % to the current block column being updated.
            %-----------------------------------------------------------
            K = kk:kk+blocksize-1;

            %-----------------------------------------------------------------
            % Compute product of block (I, J) of A with the block col K of B
            % Try replacing the following line with loops that access the
            % matrices T and B with stride one in a 1-d array layout of
            % the two arrays containing T and B, and the expressive power
            % of block indexing will become obvious.
            %------------------------------------------------------------------
            C(I, K) = C(I, K) + T'*B(J, K);

        end % kk block

      end % jj loop

   end % ii loop

