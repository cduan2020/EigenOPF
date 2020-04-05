function  [v,d] = speig(A,B,k,sigma,options)

% SPEIG	Generalized sparse eigenvalue problem
%
%    SPEIG(A) is a vector containing the 5 eigenvalues of largest
%    magnitude of the matrix A.
%
%    [V,D] = SPEIG(A) produces a diagonal matrix D of the 5
%    eigenvalues of largest magnitude of A and a full matrix V
%    whose columns are the corresponding eigenvectors so that
%    A*V = V*D.
%
%    SPEIG(A,B) is a vector containing the 5 generalized
%    eigenvalues of largest magnitude of the square matrices
%    A and B.  B must be a square matrix the same size as A,
%    and must be positive semidefinite if A is sparse.  
%
%    [V,D] = SPEIG(A,B) produces a diagonal matrix D of the 5 general-
%    ized eigenvalues of largest magnitude of A and B and a full
%    matrix V whose columns are the corresponding eigenvectors so
%    that A*V = B*V*D.
%
%    SPEIG(A,k) and [V,D] = SPEIG(A,B,k) find the k eigenvalues of
%    largest magnitude.
%
%    SPEIG(A,k,sigma) and [V,D] = SPEIG(A,B,k,sigma) find the k
%    eigenvalues closest to sigma in the following sense:
%
%    sigma                    "find the k eigenvalues ..."
%    -------------------------------------------------------------    
%    Real or complex scalar   closest to sigma in the sense of
%                               absolute value
%    'LM'                     of Largest Magnitude
%                               [Default]
%    'SM'                     of Smallest Magnitude <=> sigma = 0
%    'LR'                     of Largest Real part
%    'SR'                     of Smallest Real part
%    'BE'                     on Both Ends.  Computes k/2
%                               eigenvalues from each end of the
%                               spectrum (one more from the high
%                               and if k is odd.)
%                             
%
%    Additionally, A need not be a matrix but may be an m-file
%    or function specifying a matrix-vector product.
%    SPEIG('mvprod', ...) expects an m-file or function mvprod such that
%    mvprod(X) returns A*X for the desired operator A and a matrix X of
%    column vectors.  This option is useful when the matrix A is not
%    explicitly defined, or the action of A is predictable to compute
%    (i.e. SPEIG('fft', ...) is much faster than SPEIG(F, ...) where F is
%    the explicitly defined FFT matrix.)  In the case that the argument A
%    is a string, the user must also specify the dimension of the problem,
%    n, either in an OPTIONS structure or as the last argument in the call.
%
%    SPEIG(A, ... ,OPTIONS) operates as above with default algorithm
%    parameters replaced by values in OPTIONS, an argument created
%    with the SPEIGSET function.  See SPEIGSET and SPEIGGET for details.

%	Richard J. Radke and Dan Sorensen, 11/1/95 

        global op lbd ubd sig   % Global variables used in cheby.m
        global filtpoly iter
        tt = clock;             % Starts clock
        defk = 5;               % Default value of k
        defsig = 'DS';          % Default value of sigma
        defopt = zeros(1,8);    % Default value of options
        maxshown = 25;          % How many values are displayed in the
                                %     Command Window

%       ======  Check input values and set defaults where needed
%                   Numeric code begins at line xxx

        if nargin == 0
                error('Missing argument A.')
        end

        [m,n] = size(A);
        if m ~= n & ~isstr(A)
                error('Matrix A must be square.')
        end

        if nargin == 1
                if isstr(A)
                        if strcmp(A,'bain')
                                error(['Can you hear the shape of',...
                                        ' a Drummond?'])
                        else
                                error(['Options argument must be ',...
                                        'specified for non-matrix A.'])
                        end
                else
                        B = speye(size(A));
                        k = defk;
                        sigma = defsig;
                        options = defopt;
                end
        elseif nargin == 2
                sigma = defsig;
                if isstr(A)
                        if isscalar(B)
                                n = B;
                                options = defopt;
                                B = speye(n);
                                k = defk;
                        elseif validopt(B,1)
                                options = B;
                                n = options(1);
                                B = speye(n);
                                k = defk;
                        else
                                error(['Invalid or missing options',...
                                        ' argument (must be specified',...
                                        ' for non-matrix A.)'])
                        end
                elseif validopt(B,0)
                        options = B;
                        B = speye(n);
                        k = defk;
                elseif isscalar(B)
                        k = B;
                        B = speye(n);
                        options = defopt;
                else				
                        k = defk;
                        options = defopt;
                end
        elseif nargin == 3
                if isstr(A)
                        if isscalar(k)
                                n = k;
                                options = defopt;
                                if isscalar(B)
                                        k = B;
                                        B = speye(n);
                                        sigma = defsig;
                                else
                                        k = defk;
                                        sigma = defsig;
                                end
                        elseif validopt(k,1)
                                options = k;
                                n = options(1);
                                if isscalar(B)
                                        k = B;
                                        B = speye(n);
                                        sigma = defsig;
                                else
                                        k = defk;
                                        sigma = defsig;
                                end
                        else
                                error(['Invalid or missing options ',...
                                        'argument (must be specified',...
                                        ' for non-matrix A.)'])
                        end
                elseif validopt(k,0)
                        options = k;
                        if isscalar(B)
                                k = B;
                                B = speye(n);
                                sigma = defsig;
                        else
                                k = defk;
                                sigma = defsig;
                        end
                elseif isscalar(B)
                        sigma = k;
                        k = B;
                        B = speye(n);
                        options = defopt;
                else
                        sigma = defsig;
                        options = defopt;
                end
        elseif nargin == 4
                if isstr(A)
                        if isscalar(sigma)
                                n = sigma;
                                options = defopt;
                                if isscalar(B)
                                        sigma = k;
                                        k = B;
                                        B = speye(n);
                                else
                                        sigma = defsig;
                                end
                        elseif validopt(sigma,1)
                                options = sigma;
                                n = options(1);
                                if isscalar(B)
                                        sigma = k;
                                        k = B;
                                        B = speye(n);
                                else
                                        sigma = defsig;
                                end
                        else
                                error(['Invalid or missing options',...
                                        ' argument (must be specified',...
                                        ' for non-matrix A.)'])
                        end
                elseif validopt(sigma,0)
                        options = sigma;
                        if isscalar(B)
                                sigma = k;
                                k = B;
                                B = speye(n);
                        else
                                sigma = defsig;
                        end
                else
                        options = defopt;
                end
        elseif nargin == 5
                if isstr(A) & isscalar(options)
                        n = options(1);
                        options = defopt;
                elseif isstr(A) & validopt(options)
                        n = options(1);
                end
        end

        if ~validopt(options,1)
                error('Invalid options.')
        end

        ksave = k;
        k = k+3;
        sig = sigma;

%       The above is a trick to get faster convergence.  The
%       tail end (closest to cut off of the sort) will not
%       converge as fast as the leading end of the "wanted"
%       Ritz values.  Upon convergence of the leading ksave
%       Ritz values, k <- ksave and return.

        if any(B-speye(n))
                mess = 'eig(full(A),full(B))';
        else
                mess = 'eig(full(A))';
        end

        if n <= 100 & ~isstr(A)
                disp(' ')
                disp(['WARNING: The problem size n is small.  It will be',...
                      ' faster to use '])
                disp([mess, ' instead.  Use sparse eig when n > 100.'])
                disp(' ')
        elseif k >= .2*n & ~isstr(A)
                disp(' ')
                disp(['WARNING: The number of eigenvalues requested, k,',...
                        ' is 20% or more of the'])
                disp(['problem size n.  Either use ',mess,...
                        ' or try asking for a'])
                disp(['smaller number of eigenvalues around several',...
                        ' different shifts.'])
                disp('Use sigma = ''LM'' to get an idea of the spectrum.')
                disp(' ')
        end

        if options(2) == 0
                p = min(max(2*k,20),n);         % Default p
        else
                p = options(2);
        end
	
        if options(4) == 0
                maxit = max(300,2*n/p);         % Default maxit
        else
                maxit = options(4);
        end

        if options(5) > 0 & isstr(A)
                issym = 1;
        elseif ~isstr(A)
                if ~nnz(A-A')
                        issym = 1;
                else
                        issym = 0;
                end
        else
                issym = 0;
        end
	
        if options(3) == 0
                if issym
                        tol = 1e-10;             % Default tol
                else
                        tol = 1e-6;
                end		
        else
                tol = options(3);
        end

        %    Check to see if polynomial acceleration is requested

        if options(6) > 0 & isstr(A) & strcmp(sigma,'LR')
                sigma = 'LO';
                dopoly = 1;
        elseif options(6) > 0 & isstr(A) & strcmp(sigma,'SR')
                sigma = 'SO';
                dopoly = 1;
        elseif options(6) > 0 & ~isstr(sig)
                sigma = 'LO';
                dopoly = 1;
        else
                dopoly = 0;
        end

        if strcmp(sigma,defsig)                  %  Default sigma
                sigma = 'LM';
        end

        if isstr(sigma)
                sigma = upper(sigma);
        end

%       if strcmp(sigma,'SM') & ~isstr(A)
%               sigma = 0;
%       end

        if options(7) > 0                        % gui = 1: Display the large
                gui = 1;                         %   "Progress Report" window
                smallgui = 0;
        elseif options(7) < 0
                gui = 0;                         % smallgui = 1: Display the
		smallgui = 1;                    %   "Stop Button" window
        else
                gui = 0;
                smallgui = 0;
        end
	
        if options(8) > 0
                dispn = 1;
                format long e
        else
                dispn = 0;
        end

        if any(size(B) ~= [n n])
                error('size(B) must equal size(A) = [n n].')
        elseif isstr(A) & ~isstr(sigma) & ~dopoly
		disp(' ')
		disp(['If ',A,' is symmetric, try using an options'])
		disp('structure with issym = 1, dopoly = 1.')
		disp(' ')
		error('Non-matrix A and numeric sigma not compatible.   ')
    	elseif k > n | k < 0
                error('k must satisfy 0 <= k <= n-3.')
        elseif fix(k) ~= k
                error('k must be an integer.')
        elseif p > n
                error('p must satisfy p <= n.')
        elseif issym & imag(sigma)
                error(['A symmetric problem has real eigenvalues.  ',...
                       'Choose sigma again.'])
        end

        if ~isstr(A)
                perm = symamd(spones(A)+spones(B));
                A = A(perm,perm);
                B = B(perm,perm);
        end	
	
        B2 = B;

        [BL,b2] = chol(B);

        if b2
                disp('WARNING: B is not positive semidefinite.')
        elseif isstr(sigma)
                B = speye(n);
        end

        if p < 2*k
                disp('WARNING: p < 2*k.  Eigenvalues may not converge.')
        end

%       ======  Graphical interface
%               Numeric code resumes at line xxx

                                          if smallgui
                                              f = figure;
                                              set(f,...
                                               'units','normalized',...
                                               'position',[.5 .5 .4 .2],...
                                               'numbertitle','off',...
                                               'name','Stop Button',...
                                               'userdata',0);
                                              amain = axes(...
                                               'position',[0 0 1 1],...
                                               'visible','off');
                                              b1 = uicontrol(...
                                               'style','pushbutton',...
                                               'units','normalized',...
                                               'String','Stop!',...
                                               'callback',...
                                               'set(gcf,''userdata'',1);',...
                                               'position',[.1 .1 .8 .5]);
                                              t10 = text(...
                                               'horiz','left',...
                                               'position',[.1 .9]);
                                              t8 = text(...
                                               'horiz','left',...
                                               'position',[.1 .75]);
                                          end

                                          if gui
                                              f = figure;
                                              hs = hsv(48);
                                              hs = hs(1:32,:);

                                              set(gcf,...
                                               'units','normalized',...
                                               'name',...
                                               'Progress report',...
                                               'numbertitle','off',...
                                               'userdata',0,...
                                               'colormap',hs)
                                              amain = axes(...
                                               'position',[0 0 1 1],...
                                               'visible','off');

                                              t1 = text(...
                                               'string',...
                                               ['k: ',num2str(k-3)],...
                                               'horiz','left',...
                                               'position',[.1 .9]); 
                                              t2 = text(...
                                               'string',...
                                               ['n: ',num2str(n)],...
                                               'horiz','left',...
                                               'position',[.1 .8]);
                                              t3 = text(...
                                               'string',...
                                               ['p: ',num2str(p)],...
                                               'horiz','left',...
                                               'position',[.1 .7]); 
                                              t4 = text(...
                                               'string',['maxit: ',...
                                               num2str(maxit)],...
                                               'horiz','left',...
                                               'position',[.1 .6]);
                                              t5 = text(...
                                               'horiz','left',...
                                               'position',[.1 .3],...
                                               'string','Initializing...');
                                              t6 = text(...
                                               'horiz','left',...
                                               'position',[.1 .4], ...
                                               'string',...
                                               ['Elapsed time: ',...
                                               num2str(etime(clock,tt))]);
                                              t7 = text(...
                                               'horiz','left',...
                                               'position',[.5 .35], ...
                                               'string',...
                                               ['Tolerance: ',...
                                               num2str(tol)]);
                                              t8 = text(...
                                               'horiz','left',...
                                               'position',[.5 .3]);
                                              t10 = text(...
                                               'horiz','left',...
                                               'position',[.1 .5]);
	
                                              a1 = axes(...
                                               'units','normalized',...
                                               'position',[.4 .5 .5 .4],...
                                               'nextplot','add');
	
                                              for i=1:16;
                                               eval(['p',num2str(i),...
                                                     '=plot(0,0,''.'');'])
                                               eval(['set(p',num2str(i),...
                                                     ',''markersize'',',...
                                                     '20,',...
                                                     '''color'',',...
                                                     'hs(',num2str(i),',:))'])
                                              end
				
                                              if isstr(sigma)
                                               n1 = ['The ',num2str(ksave),...
                                                     ' eigenvalues'];
                                              if strcmp(sigma,'SM')
                                               n1 = [n1,' of smallest',...
                                                     ' magnitude.'];
                                              elseif strcmp(sigma,'SR') |...
                                                     strcmp(sigma,'SO')
                                               n1 = [n1,' of smallest',...
                                                     ' real part.'];
                                              elseif strcmp(sigma,'LM')
                                               n1 = [n1,' of largest',...
                                                     ' magnitude.'];
                                              elseif strcmp(sigma,'BE')
                                               n1 = [n1,' on both ends.'];
                                              else
			                       n1 = [n1,' of largest',...
                                                     ' real part.'];
                                              end

                                              else
                                               n1 = ['The ',num2str(ksave),...
                                                     ' eigenvalues',...
                                                     ' closest to ',...
                                                     num2str(sigma)];
                                              end

                                              title(n1)
                                              z = get(gca,'title');
                                              ex = get(z,'extent');
                                              z1=text(...
                                               'string','Unconverged',...
                                               'units','normalized',...
                                               'color',hs(16,:),...
                                               'position',[ex(1) 1],...
                                               'horiz','left');
                                              z2=text(...
                                               'string','Converged',...
                                               'units','normalized',...
                                               'color',hs(1,:),...
                                               'pos',[ex(1)+ex(3) 1],...
                                               'horiz','right');
                                              ex1 = get(z1,'extent');
                                              ex2 = get(z2,'extent');
                                              lend = ex1(1)+ex1(3);
                                              uend = ex2(1);
                                              sp = linspace(lend,uend,6);
                                              text(sp(2),1,'-',...
                                               'units','normalized',...
                                               'color',hs(13,:))
                                              text(sp(3),1,'-',...
                                               'units','normalized',...
                                               'color',hs(10,:))
                                              text(sp(4),1,'-',...
                                               'units','normalized',...
                                               'color',hs(7,:))
                                              text(sp(5),1,'>',...
                                               'units','normalized',...
                                               'color',hs(4,:))
			
                                              xlabel('Real part')
                                              ylabel('Imaginary part')
                                              set(a1,'nextplot','replace')

                                              b1 = uicontrol(...
                                               'style','pushbutton',...
                                               'units','normalized',...
                                               'String','Stop!',...
                                               'callback',...
                                               ['set(gcf,',...
                                                '''userdata'',1);'],...
                                               'position',[.1 .1 .2 .1]);
                                              drawnow

                                           end

                                           if gui
                                              set(t5,'string',...
                                                'Selecting initial vector')
                                              drawnow
                                           end

%       ======  Resume numerical code

        info = 0;
        v = rand(n,1)- .5*ones(n,1) ;
        p = p - k;
        psave = p;
        op = A;

        beta = 1.0/sqrt(v(:,1)'*B*v(:,1));
        v(:,1) = v(:,1)*beta;

%                                          ======  Graphical interface
                                           if gui
                                              set(t5,'string',...
                                               'Backsolving second vector')
                                              drawnow
                                           end

%       Compute w = Av;

        if ~isstr(A) & ~isstr(sigma),
                if gui
                        set(t5,'string','Factoring and backsolving')
                        drawnow
                end
                [L,U] = lu(A - sigma*B);
                v(:,1) = U \ (L \ (B*v(:,1)));
                beta = 1.0/sqrt(v(:,1)'*B*v(:,1));
                v(:,1) = v(:,1)*beta;
                v(:,2) = U \ (L \ (B*v(:,1)));
        elseif isstr(A) & b2
                v(:,2) = B \ feval(A,v(:,1));
        elseif isstr(A) & ~b2
                v(:,2) = BL \ feval(A,(BL'\v(:,1)));
        elseif ~b2
                v(:,2) = BL \ (A*(BL'\v(:,1)));
        else
                v(:,2) = B \ (A*v(:,1));
        end

%                                          ======  Graphical interface

                                           if gui
                                              set(t5,'string',['Performing',...
						' iterative refinement'])
                                              set(t6,'string',...
                                                ['Elapsed time: ',...
                                              num2str(etime(clock,tt))]);
                                                drawnow
                                           end

                                           if gui | smallgui
                                              pause(0);
                                              if get(gcf,'userdata')
                                                 close
                                                 error(['Stop button was',...
                                                 'clicked... no ',...
                                                 'eigenvalues computed.'])
                                              end
                                           end

        beta = 1.0/sqrt(v(:,1)'*B*v(:,1));
        v(:,1) = v(:,1)*beta;

        alpha = v(:,1)'*B*v(:,2);
        h(1,1) = alpha;

%       Compute the residual in the second column of V

        v(:,2) = v(:,2) - alpha*v(:,1);

%       Perform one step of iterative refinement to correct any
%       orthogonality problems
	
        alpha = v(:,1)'*B*v(:,2);
        v(:,2) = v(:,2) - alpha*v(:,1);
	
        h(1,1) = h(1,1) + alpha;

%                                          ======  Graphical interface
		
                                           if gui
                                              set(t6,'string',...
                                               ['Elapsed time: ', ...
                                               num2str(etime(clock,tt))])
                                              set(t5,'string', ...
                                               ['Computing the initial',...
                                               ' Arnoldi sequence'])
                                              drawnow
                                           end

                                           if gui | smallgui
                                              pause(0);
                                              if get(gcf,'userdata')
                                                 close
                                                 error(['Stop button was',...
                                                  'clicked... no ',...
                                                  'eigenvalues computed.'])
                                              end
                                           end

%       Compute k steps of the Arnoldi sequence

        kstart = 1;
        ritz = 1.0;
        kp1 = k + 1;
        kend = k + p;
        k1 = 1;

        if ~isstr(A) & ~isstr(sigma)
                [v,h] =  arnold2(k1,k,L,U,B,v,h,tol);
        elseif ~b2
                [v,h] =  arnold(k1,k,A,BL,v,h,1);
        else
                [v,h] =  arnold(k1,k,A,B,v,h,0);
        end

%                                          ======  Graphical interface
			
                                           if gui
                                              set(t6,'string', ...
                                               ['Elapsed time: ',...
                                               num2str(etime(clock,tt))])
                                              drawnow
                                           end
	
                                           if gui | smallgui
                                              pause(0);
                                              if get(gcf,'userdata')
                                                 close
                                                 error(['Stop button was',...
                                                  'clicked... no ',...
                                                  'eigenvalues computed.'])
                                              end
                                           end
	
%       Now update the Arnoldi sequence in place

        iter = 0;
        knew = k;
        ritzes = zeros(ksave,1);
        ritzests = ones(ksave,1);
        stopcrit = 1;
        beta = 1;
        betanew = 1;
        residest = 1;
        if ksave == 1
             second = 1;
        else
             second = 2;
        end
                                           if gui | smallgui
                                              stopclick =...
                                               get(gcf,'userdata');
                                           else
                                              stopclick = 0;
                                           end

%       ======  MAIN LOOP

        while ((stopcrit > tol & iter < maxit) | iter < 2) & ~stopclick

                iter = iter + 1;
		
%                                          ======  Graphical interface
		
                                           if gui
                                              set(t5,'string',...
                                               ['Extending the ',...
                                               'Arnoldi basis'])
                                              set(t10,'string',...
                                               ['Current Iteration: ',...
                                               num2str(iter)]);
                                              set(t6,'string',...
                                               ['Elapsed time: ',...
                                               num2str(etime(clock,tt))]);
                                              drawnow
                                           elseif smallgui
                                              set(t10,'string',...
                                               ['Current Iteration: ',...
                                               num2str(iter)]);
                                              drawnow
                                           end

%       Compute p additional steps of the Arnoldi sequence

                kold = k;
                k = knew;
		
                if ~isstr(A) & ~isstr(sigma)
                        [v,h,info] = arnold2(k,kend,L,U,B,v,h,tol);
                elseif ~b2
                        [v,h,info] = arnold(k,kend,A,BL,v,h,1);
                else
                        [v,h,info] =  arnold(k,kend,A,B,v,h,0);
                end

%               If we're doing polynomial acceleration, when the Ritz estimates
%               on the extreme values of the spectrum are good enough, then
%               apply the accelerant polynomial instead.

                if dopoly == 1 & (ritzests(1) < .1) & (ritzests(second) < .1)

                        disp('Starting polynomial acceleration')
                	
                        A = 'accpoly';
                        if gui
                                set(t5,'string','Restarting the basis')
                                drawnow
                        end
                        eigh = eig(h);			
                        iter = 1;			

                        %  Restart the basis with the largest/smallest
                        %  ritz vector, as appropriate.

                        v = v(:,2);
                        beta = 1.0/sqrt(v(:,1)'*B*v(:,1));
                        v(:,1) = v(:,1)*beta;

                        if isstr(A) & b2
                                v(:,2) = B \ feval(A,v(:,1));
                        elseif isstr(A) & ~b2
                                v(:,2) = BL \ feval(A,(BL'\v(:,1)));
                        elseif ~b2
                                v(:,2) = BL \ (A*(BL'\v(:,1)));
                        else
                                v(:,2) = B \ (A*v(:,1));
                        end

                        beta = 1.0/sqrt(v(:,1)'*B*v(:,1));
                        v(:,1) = v(:,1)*beta;
                        alpha = v(:,1)'*B*v(:,2);
                        h = alpha;
                        v(:,2) = v(:,2) - alpha*v(:,1);
                        alpha = v(:,1)'*B*v(:,2);
                        v(:,2) = v(:,2) - alpha*v(:,1);
                        h = h + alpha;
			
                        if ~b2
                                [v,h,info] =  arnold(1,kend,A,BL,v,h,1);
                        else
                                [v,h,info] =  arnold(1,kend,A,B,v,h,0);
                        end

                        k1 = 1;
                        knew = k;
                        ritzes = zeros(ksave,1);
                        dopoly = -1;
                        residest = 1;
                        sigma = 'LR';
			
                                           if gui
                                              set(t5,'string',...
                                               ['Calculating accelerant ',...
                                                'polynomial'])
                                              axes(a1)
                                              cla
                                              l1 = linspace(lbd,ubd);
                                              l2 = accpoly1(l1);
                                              plot(l1,l2)
                                              hold on
                                              plot(sort(eigh), ...
                                                 zeros(size(eigh)),'*')
                                              title(['Accelerant ',...
                                                 'polynomial applied',...
                                                 ' to original ritz values']);
                                              xlabel('Ritz value')
                                              ylabel('Value of polynomial')
                                              set(t10,'string',...
                                               'Polynomial acceleration');
                                              drawnow
                                           elseif smallgui
                                              set(t10,'string',...
                                               ['Starting polynomial',...
                                                ' acceleration']);
                                              drawnow
                                           end
                end

                k = kold;
                                           if gui
                                              set(t6,'string',...
                                               ['Elapsed time: ',...
                                               num2str(etime(clock,tt))]);
                                              drawnow
                                           end

%        Compute p shifts based on sigma

                t4b = clock;

                                           if gui
                                              set(t5,'string',...
                                               'Computing shifts')
                                              drawnow
                                           end

%               If A is symmetric, keep H tridiagonal (to avoid spurious
%               imaginary parts).

                if issym
                        for i=1:kend-2
                                h(i,i+2:kend) = zeros(1,kend-i-1);
                                hav = mean([h(i,i+1),h(i+1,i)]);
                                h(i,i+1) = hav;
                                h(i+1,i) = hav;
                        end
                        hav = mean([h(kend,kend-1),h(kend-1,kend)]);
                        h(kend,kend-1) = hav;
                        h(kend-1,kend) = hav;
                end	

                if dispn, clc, end

                [w,q1] = shftit(h,kstart,kend,sigma);

%       Update the command window with current eigenvalue estimates

                ritzold = ritzes;

                if ~isstr(sigma)
                        ritzes = sigma + ...
                                1./w((kend-kstart+1):-1:(kend-ksave+1));
                        ritzes = [(sigma+1./eig(h(1:kstart-1,1:kstart-1)));...
                                ritzes];
                else
                        ritzes = w((kend-kstart+1):-1:(kend-ksave+1));
                        ritzes = [eig(h(1:kstart-1,1:kstart-1));ritzes];
                end

                if dispn
                   disp(ritzes(1:min(maxshown,length(ritzes))))
                end

                if iter == 1 & dopoly > 0
                        lbd = min(ritzes);
                        ubd = max(ritzes);
                elseif dopoly > 0
                        lbd = min(lbd, min(ritzes));
                        ubd = max(ubd, max(ritzes));
                end

                [m1,m2]=size(q1);
                ritz = norm(q1(m1,p+2:m1));
                betanew = sqrt(v(:,kend+1)'*B*v(:,kend+1));
                ritznew = betanew*q1(m1,1:m1);
                jj = m1;
                kconv = 0;
                while(jj > 0),
                        if (abs(ritznew(jj)) < tol),
                                jj = jj - 1;
                                kconv = kconv+1;
                        else
                                jj = -1;
                        end
                end
                kkconv = kconv;

%               The while loop counts the number of converged ritz values.

                ritzests = [w abs(ritznew)'];
                ritzests = ritzests(size(ritzests,1) : -1 : ...
                size(ritzests,1)-ksave+1,2);

                stopcrit = max(ritzests);
                residest = norm(ritzests);		

%       At first, we use the Ritz estimates to estimate convergence.
%       However, as the algorithm converges, the Ritz estimates become
%       poor estimates of the actual error in each Ritz pair.  So when the
%       Ritz estimates become too small, we actually form the matrix of
%       errors || AV - VD || where V are the estimates for the eigenvectors
%       and eigenvalues.  This is expensive computationally but ensures
%       that the user gets the desired eigenpairs to the requested
%       tolerance.

                if max(ritzests) < tol*max(norm(h),1)

                   if ~b2 & isstr(sigma)
                      vee = BL \ v(:,1:kend);
                   else
                      vee = v(:,1:kend);
                   end

                   if ~isstr(A)
                      errmat = A*vee*q1(1:kend, kend:-1:kend-ksave+1) - ...
                          B2*vee*q1(1:kend,kend:-1:kend-ksave+1)*diag(ritzes);
                   else
                      errmat=feval(A,vee)*q1(1:kend,kend:-1:kend-ksave+1)-...
                          B2*vee*q1(1:kend,kend:-1:kend-ksave+1)*diag(ritzes);
                   end

                   if ~isstr(A)
                      residest = norm(errmat,1)/norm(A,1);
                   else
                      residest = norm(errmat,1)/norm(feval(A,vee)* ...
                                 q1(1:kend,kend:-1:kend-ksave+1),1);
                   end

                   for ii = 1:length(ritzes);
                      ritzests(ii) = norm(errmat(:,ii));
                   end
	
                   stopcrit = residest;

                end
%                                          ======  Graphical interface
%                                          Numerical code resumes at line xxx

                                           if gui
                                              set(t5,'string', ...
                                               'Computing residual')
                                              drawnow
                                              set(t6,'string',...
                                               ['Elapsed time: ',...
                                               num2str(etime(clock,tt))]);
                                              set(t8,'string',...
                                               ['Residual norm: ',...
                                               num2str(residest)]);
                                              drawnow
                                           elseif smallgui
                                              set(t8,'string',...
                                               ['Residual norm: ',...
                                               num2str(residest)]);
                                              drawnow
                                           end

%                                          If gui activated, display
%                                          the converging eigenvalues

                                           if gui & dopoly >= 0
                                              set(t5,'string',...
                                               'Generating graph')
                                              drawnow
	
                                              for ii=1:length(ritzes)
                                                 ritzmags(ii) = floor(...
                                                   log10(ritzests(ii)... 
                                                   +eps)) + 17;
                                              end

                                              ritzmags(ritzmags<1)=...
                                               zeros(size(find(ritzmags<1)));
                                              ritzmags(ritzmags>32)=...
                                               32*...
                                               ones(size(find(ritzmags>32)));

                                              for ii=1:16
                                                 idata = ritzes(...
                                                  find(ritzmags == ii));
                                                 eval(['set(p',...
                                                  num2str(max(1,ii)),...
                                                  ',''xdata''',...
                                                  ',real(idata)',...
                                                  ',''ydata''',...
                                                  ',imag(idata))'])
                                              end
			
                                              r1 = min(real(ritzes));
                                              r2 = max(real(ritzes));
                                              d1 = (r2-r1)/18;
                                              i1 = min(imag(ritzes));
                                              i2 = max(imag(ritzes));
                                              d2 = (i2-i1)/18;
                                              if i1 == i2 & i1 == 0
                                                i1 = -1;
                                                i2 = 1;
                                              end

                                              set(a1,...
                                               'xlim',[r1-d1,r2+d1],...
                                               'ylim',[i1-d2,i2+d2])
                                              set(t6,...
                                               'string',['Elapsed time: ',...
                                               num2str(etime(clock,tt))]);

                                              drawnow
                                           end

         if (stopcrit > tol | iter < 2)

%        Apply the p implicit shifts if convergence has not yet 
%        happened.  Otherwise don't apply them and get out of the
%        loop on next loop test.  We need to keep same test here as
%        in the main loop test to avoid applying shifts and then
%        quitting, which would lead to a wrong size factorization
%        on return.

                                           if gui
                                              set(t5,'string',...
                                               'Applying shifts')
                                              drawnow
                                           end
            
%        If some ritz values have converged then
%        adjust k and p to move the "boundary"
%        of the filter cutoff.

                if kconv > 0
                        kk = ksave + 3 + kconv;
                        p = max(ceil(psave/3),kend-kk);
                        k = kend - p;
                end

                if any(any(imag(v))) | any(any(imag(h)))
                        [v,h,knew] = apshft1(v,h,w,k,p);
                else
                        [v,h,knew] = apshft2(v,h,real(w),imag(w),k,p);
                end

                betanew = sqrt(v(:,kp1)'*B*v(:,kp1));
			
         end

%                                          ======  Graphical interface
                                           if gui
                                              set(t6,'string',...
                                               ['Elapsed time: ',...
                                               num2str(etime(clock,tt))]);
                                           end

                                           if gui | smallgui
                                              pause(0);
                                              stopclick = get(gcf,'userdata');
                                           end

         end             %  End of Arnoldi iteration  (main loop)

                                           if gui
                                              set(t5,'string',...
                                               ['Convergence... ',...
                                                'cleaning up'])
                                              drawnow
                                           end

%       Compute the eigenvalues and eigenvectors of h
	
        [w,q1] = shftit(h,1,kend,sigma);

        k = ksave;
        p = psave;
     	
                                           if gui
                                              set(t6,'string',...
                                               ['Elapsed time: ',...
                                               num2str(etime(clock,tt))])
                                           end 

%       Transform the converged eigenvalues back to 
%       the original problem and return them in the diagonal
%       k by k matrix d.

        if iter >= maxit
                disp('WARNING: Maximum number of iterations exceeded!')
        end

%       Set v equal to the wanted eigenvectors

                                           if gui
                                              set(t6,'string',...
                                               ['Elapsed time: ',...
                                               num2str(etime(clock,tt))])
                                           end
	
        v = v(:,1:kend) * q1(1:kend, kend:-1:kend-k+1);

        if ~b2 & isstr(sigma)
                v = BL \ v;
        end

        for i=1:k
                v(:,i)=v(:,i)/norm(v(:,i));
        end	
	
%       In polynomial acceleration, we recover the eigenvalues by Rayleigh
%       quotients.  Otherwise, the eigenvalues are recovered from the
%       last set of shifts w.

        if dopoly == -1
                d = [];
                for i=1:ksave
                    d = [d; ( (v(:,i)' * feval(op,v(:,i)) ) ./ ...
                         (v(:,i)' * B2 * v(:,i)))];
                end
                d = diag(d);
        elseif ~isstr(sigma)
                t = sigma + 1./w;
                d = diag(t(kend:-1:kend-k+1));
        else
                d = diag(w(kend:-1:kend-k+1));
        end
	
        if ~isstr(A)
                residnorm = norm(A*v - B2*v*d)
        else
                residnorm = norm(feval(op,v) - B2*v*d)
        end
	
                                          if gui
                                             set(t6,'string',...
                                              ['Elapsed time: ',...
                                              num2str(etime(clock,tt))])
                                          end

        if ~isstr(A)
                v(perm,:) = v;
        end

        if nargout <= 1
                v = diag(d);
        end

        if stopclick
                disp(['WARNING:  The stop button was pressed',...
                       ' before the routine was completed.'])
                disp('')
                disp(['This may result in inaccurate eigenvalues',...
                       ' and eigenvectors, and/or a loss'])
                disp('in B-orthogonality of V.')
                disp('')
        end
        if gui | smallgui,close(gcf),end