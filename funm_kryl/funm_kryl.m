function [f,out,param] = funm_kryl(A,b,param)
%FUNM_KRYL Computes Krylov approximation to f(A)*b.
%   f = FUNM_KRYL(A,b,param) is an approximation to f(A)*b
%   computed by a restarted Krylov subspace method.
%
%   see index.htm for details
%
%   09.05.07 - 24.11.08, M. Eiermann, S. Guettel


if nargin < 3
    param = struct;
end;

[param,modified] = param_init(param);

if (modified)
    disp('Some of the parameters were incorrect. The algorithm will now run with a corrected setting.');
end;

if (param.waitbar) 
    hand = waitbar(0,'Please wait...'); 
end;

beta = norm(b);
v = (1/beta)*b;

% allocate memory [ >= m+max(ell) ]
global V_big                       
alloc = param.restart_length + 20;
V_big = zeros(length(b),alloc);     

if (param.bound)
    alpha_1 = +inf; 
    alpha_2 = -inf;
else
    alpha_1 = 0; 
    alpha_2 = 0; 
end;

m = param.restart_length;
ell = 0;                             % thick restart param
ell_prev = 0;

if (param.H_full)
    f = zeros(size(b));
    H_full = [];
else
    nr_single = length(param.function.single_poles);   
    poles = [ param.function.single_poles; param.function.conj_poles ];
    coeff = [ param.function.single_coeff; param.function.conj_coeff ];
    f = param.function.absterm * b;
    B_bar = zeros(m+2,length(poles)); % with this initialization:
    B_bar(end-2,:) = 1;               % sweep_1 = sweep_n for all n>1
    s = 1;                            % -- " --
end;

if (param.V_full)
    out.V_full = [];
end;

out.stop_condition = '0 - maximal # of restarts exceeded';

% restart loop starts here
for k = 1:param.max_restarts,

    % check whether a stop condition is satisfied
    if str2double(out.stop_condition(1))
        break;
    end;
    
    % waitbar and cputime
    if (param.waitbar)
        waitbar(k/param.max_restarts,hand); 
    end;
    out.time(k) = cputime;

    % compute A-invariant subspace of prev. cycle (thick restart)
    if (~isempty(param.thick) && k > 1),
        ell_prev = ell;
        [ ell,U,T,D ] = param.thick( H,U,T,D,k,param.function );
        out.thick_replaced{k-1} = D(1:ell);
        out.thick_interpol{k-1} = D(ell+1:end);
        if ell,
            U = U(:,1:ell);
            V_big(:,1:ell) = V_big(:,1:m+ell_prev)*U;
            H_hat = U'*H*U;
            H = [H_hat; eta*U(end,:)];
        else
            H = [];
        end;
    else
        H = [];
    end
    
    V_big(:,ell+1) = v;
    
    % compute/extend Krylov decomposition    
    if param.hermitian,
        [ v,H,eta,breakdown ] = lanczos( A,m+ell,H,ell+1,param );
    else
        [ v,H,eta,breakdown ] = arnoldi( A,m+ell,H,ell+1,param );
    end;
    
    if isfinite(param.harmonic_target),
        z = abs(eta)^2*(((H - param.harmonic_target*eye(size(H)))')\...
            unit(m+ell,m+ell));
        H = H + z*unit(m+ell,m+ell)';
        v = eta*v - V_big(:,1:m+ell)*z;
        eta = sqrt(param.inner_product(v,v));
        v = (1/eta)*v;
    end;
    
    if (breakdown)
        error(['Arnoldi breakdown in sweep ' num2str(k)]);
    end;
    
    % store full Krylov basis? compute bounds?
    if param.V_full,
        out.V_full = [ out.V_full , V_big(:,1:m+ell) ];
    end;

    % Schur form of H (for thick restarts or error estimator)
    if (~isempty(param.thick) || param.bound),
        if isreal(H),
            [U,T] = schur(H,'real');
        else
            [U,T] = schur(H);
        end;
        D = ordeig(T);
    end;
    
    if (param.bound)
        alpha_1 = min([real(D) ; alpha_1]);
        alpha_2 = max([real(D) ; alpha_2]);
    end;
 
    % evaluate matrix function of H times unit coordinate vector
    if (param.H_full),
        H_full = blkdiag(H_full,H);
        if k>1,
            H_full(end-m+1,end-m-ell) = s;
        end;
        H_bar = blkdiag(H_full,[alpha_1 0 ; 1 alpha_2 ]);
        H_bar(end-1,end-2) = eta;
        if isstruct(param.function)
            h = eval_parfrac(param.function,H_bar,eye(size(H_bar,1),1));
        else
            h = param.function(H_bar);
            h = h(:,1);
        end;
        h = h(end-m-ell-1:end);
    else
        H_bar = blkdiag(H, [alpha_1 0 ; 1 alpha_2 ]);
        H_bar(m+ell+1,m+ell) = eta;
        for p = 1:length(poles),
            q_bar = s*B_bar(m+ell_prev,p)*unit(ell+1,m+ell+2);
            B_bar(1:m+ell+2,p) = ( poles(p)*eye(size(H_bar)) - H_bar ) \ q_bar;
        end;
        h = B_bar(1:m+ell+2,1:nr_single)*coeff(1:nr_single);
        if(nr_single<length(poles))
            h = h + 2*real( B_bar(1:m+ell+2,nr_single+1:end)*coeff(nr_single+1:end) );
        end;
    end;

    s = eta;
    
    % workaround due to matlab 'bug' (copies large submatrices)
    h_big = h(1:m+ell,1); 
    if size(V_big,2) > length(h_big),
        h_big(size(V_big,2),1) = 0;
    end;
    
    % update Krylov approximation
    f = beta*(V_big*h_big) + f;

    out.appr(:,k) = f;
    out.update(k) = beta*norm(h_big);  % norm of update

    % get cpu-time
    out.time(k) = cputime - out.time(k);
    
    % check stopping conditions
    if(param.bound),                            % stop by error bound?
        out.bound_1(k) = beta*abs( h(end-1) );
	    if isnumeric(A),
        	w = A*v;
    	else
        	w = A(v);
     	end;
        out.bound_2(k) = beta*norm( (h(end-1) - alpha_1*h(end))*v + h(end)*w ); % TODO: reuse w in next cycle
        if isempty(param.exact) && stopcondition(out.bound_2(2:end)./out.bound_2(1:end-1) >  param.min_decay)
            out.stop_condition = '4 - linear convergence rate of upper bound > min_decay';
        end;
        if out.bound_2(k) < param.stopping_accuracy && isempty(param.exact)
            out.stop_condition = '2 - upper bound below stopping accuracy';
        end;
    end;

    if  stopcondition(out.update < param.stopping_accuracy)  % stop by norm of update?
        out.stop_condition = '5 - norm of updates decayed below stopping accuracy';
    end;

    if ~isempty(param.exact),                   % stop by absolute error?
        out.err(k) = norm(f - param.exact);
        out.err(k) = out.err(k);
        if stopcondition(out.err(2:end)./out.err(1:end-1) >  param.min_decay)
            out.stop_condition = '3 - linear convergence rate of absolute error > min_decay';
        end;
        if out.err(k) < param.stopping_accuracy
            out.stop_condition = '1 - absolute error below stopping accuracy';
        end;
    end
        
end 
% restart loop ends here

if (param.H_full)
    out.H_full = H_full;
end;

if (param.waitbar) 
    close(hand);
end;

clear V_big




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function r = stopcondition(v)
% stop condition
% checks whether the entries of v are of the form 1,...,1,0,...0,1,1

r = 0;

if length(v) < 2
    return;
end;

if all(v(end-1:end)) && ~all(v)
    r = 1;
    return;
end;

function u = unit(pos,dim)

u = zeros(dim,1);
u(pos) = 1;