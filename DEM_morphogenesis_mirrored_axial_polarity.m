function DEM = DEM_morphogenesis_mirrored_axial_polarity
% This routine illustrates self-assembly or more for genesis under active
% inference (free energy minimisation).  It exploits the fact that one can
% express a systems (marginal) Lyapunov function in terms of a variational
% free energy.  This means that one can prescribe an attracting set in
% terms of the generative model that defines variational free energy.  In
% this example, the attracting set is a point attractor in the phase space
% of a multi-celled organism: where the states correspond to the location
% and (chemotactic) signal expression of 16 cells.  The generative model
% and process are remarkably simple; however, the ensuing migration and
% differentiation of the 8 cells illustrates self-assembly – in the sense
% that each cell starts of in the same location and releasing the same
% signals.  In essence, the systems dynamics rest upon each cell inferring
% its unique identity (in relation to all others) and behaving in accord
% with those inferences; in other words, inferring its place in the
% assembly and behaving accordingly.  Note that in this example there are
% no hidden states and everything is expressed in terms of hidden causes
% (because the attracting set is a point attractor)  Graphics are produced
% illustrating the morphogenesis using colour codes to indicate the cell
% type – that is interpreted in terms of genetic and epigenetic
% processing.
%
% This particular version shows mirrored anterior/posterior polarity 
% (head and tail positioning).
% We introduced a positive squared gradient in the generative process for 
% the signaling inputs that each cell receives from its environment, 
% resulting in double head formation. This option can be turned on by 
% seting the global boolian "double_head" to 1 in line 66.
% We also introduced a negative squared gradient in the generative process  
% for the sensory inputs that each cell receives from its environment, 
% resulting in double tail formation. This option can be turned on by 
% seting the global boolian "double_tail" to 1 in line 67.
% The postions of cells and all subsequent characteristics were flipped
% along the vertical axis with respect to its original form for illustration
% purposes (have the head of the cell ensemble up, not down) in all
% graphics.
% Published in:
% F. Kuchling, K. Friston, G. Georgiev, M. Levin, Physics of Life Review
% 2019, “Morphogenesis as bayesian inference: A variational approach to  
% pattern formation and control in complex biological systems” 
% _________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% By Franz Kuchling
% Modified from Karl Friston, as published in 
% K. Friston, M. Levin, B. Sengupta, and G. Pezzulo, Royal Society 2015, 
% “Knowing one ' s place : a free-energy approach to pattern regulation,” 
% $Id: DEM_morphogenesis.m 7145 2017-07-31 13:57:39Z karl $
 
% preliminaries
%--------------------------------------------------------------------------
clear global
rng('default')
SPLIT    = 0;                              % split: 1 = upper, 2 = lower
N        = 32;                             % length of process (bins)
 
% generative process and model
%==========================================================================
M(1).E.d  = 1;                             % approximation order
M(1).E.n  = 2;                             % embedding order
M(1).E.s  = 1;                             % smoothness
 
global tau;  tau = 0.25/log(2) ;            % sensitivity s=1-exp(-t/tau)

% only one of the following can be set to 1, the other must be zero
global double_head; double_head=1; % set this to 1 for the abberrant signaling phenotype
global double_tail; double_tail=0; % set this to 1 for the resuce of abberrant signaling phenotype

% priors (prototype)
%--------------------------------------------------------------------------
T =[0 0 2 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0 0 0;
    2 0 0 0 0 0 4 0 4 0 3;
    0 0 0 0 1 0 0 0 0 0 0;
    0 0 2 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0];
     
p(:,:,1) = T > 0;
p(:,:,2) = T == 2 | T == 1;
p(:,:,3) = T == 3 | T == 1;
p(:,:,4) = T == 4;
 
[y,x] = find(p(:,:,1));
P.x   = spm_detrend([x(:) y(:)])'/2; 
 
% signalling of each cell type
%--------------------------------------------------------------------------
n     = size(P.x,2);                      % number of cells
m     = size(p,3);                        % number of signals
j     = find(p(:,:,1));
for i = 1:m
    s        = p(:,:,i);
    P.s(i,:) = s(j);
end
P.s   = double(P.s);
P.c   = morphogenesis(P.x,P.s);           % signal sensed at each position
 
% initialise action and expectations
%--------------------------------------------------------------------------
v     = randn(n,n)/8;                     % states (identity)
g     = Mg([],v,P);
a.x   = g.x;                              % action (chemotaxis)
a.s   = g.s;                              % action (signal release)
 
 
% generative process 
%==========================================================================
R     = spm_cat({kron(eye(n,n),ones(2,2)) []; [] kron(eye(n,n),ones(4,4));
                 kron(eye(n,n),ones(4,2)) kron(eye(n,n),ones(4,4))});

% level 1 of generative process
%--------------------------------------------------------------------------
G(1).g  = @(x,v,a,P) Gg(x,v,a,P);
G(1).v  = Gg([],[],a,a);
G(1).V  = exp(16);                         % precision (noise)
G(1).U  = exp(2);                          % precision (action)
G(1).R  = R;                               % restriction matrix
G(1).pE = a;                               % form (action)
 
 
% level 2; causes (action)
%--------------------------------------------------------------------------
G(2).a  = spm_vec(a);                      % endogenous cause (action)
G(2).v  = 0;                               % exogenous  cause
G(2).V  = exp(16);
 
 
% generative model
%==========================================================================
 
% level 1 of the generative model: 
%--------------------------------------------------------------------------
M(1).g  = @(x,v,P) Mg([],v,P);
M(1).v  = g;
M(1).V  = exp(3);
M(1).pE = P;
 
% level 2: 
%--------------------------------------------------------------------------
M(2).v  = v;
M(2).V  = exp(-2);
 
 
% hidden cause and prior identity expectations (and time)
%--------------------------------------------------------------------------
U     = zeros(n*n,N);
C     = zeros(1,N);
 
% assemble model structure
%--------------------------------------------------------------------------
DEM.M = M;
DEM.G = G;
DEM.C = C;
DEM.U = U;
 
% solve
%==========================================================================
DEM   = spm_ADEM(DEM);
spm_DEM_qU(DEM.qU,DEM.pU);


% split half simulations
%==========================================================================
if SPLIT
    
    % select (partially diferentiated cells to duplicate
    %----------------------------------------------------------------------
    t    = 8;
    v    = spm_unvec(DEM.pU.v{1}(:,t),DEM.M(1).v);
    if SPLIT > 1
        [i j] = sort(v.x(1,:), 'ascend');
    else
        [i j] = sort(v.x(1,:),'descend');
    end
    j    = [j(1:n/2) j(1:n/2)];
    
    % reset hidden causes and expectations
    %----------------------------------------------------------------------
    v    = spm_unvec(DEM.qU.v{2}(:,t),DEM.M(2).v);
    g    = spm_unvec(DEM.qU.v{1}(:,t),DEM.M(1).v);
    a    = spm_unvec(DEM.qU.a{2}(:,t),DEM.G(1).pE);
    
    v    = v(:,j);
    g.x  = g.x(:,j);
    g.s  = g.s(:,j);
    g.c  = g.c(:,j);
    a.x  = a.x(:,j) + randn(size(a.x))/512;
    a.s  = a.s(:,j) + randn(size(a.s))/512;
    
    DEM.M(1).v = g;
    DEM.M(2).v = v;
    DEM.G(2).a = spm_vec(a);
    
    % solve
    %----------------------------------------------------------------------
    DEM   = spm_ADEM(DEM);
    spm_DEM_qU(DEM.qU,DEM.pU);
    
end



 
% Graphics
%==========================================================================
 
% expected signal concentrations
%--------------------------------------------------------------------------
subplot(2,2,2); cla
A     = max(abs(P.x(:)))*3/2;
h     = 2/3;
 
x     = linspace(-A,A,32);
[x,y] = ndgrid(x,x);
x     = spm_detrend([x(:) y(:)])';
c     = morphogenesis(P.x,P.s,x);
c     = c - min(c(:));
c     = c/max(c(:));

x1=rot90(rot90(x(1,:)));  % rotate display from original paper Friston 2015
x2=rot90(rot90(x(2,:)));  % rotate display from original paper Friston 2015

for i = 1:size(c,2)
    col = c(end - 2:end,i);
%     plot(x(2,i),x(1,i),'.','markersize',32,'color',col); hold on
    plot(x2(i),x1(i),'.','markersize',32,'color',col); hold on  % rotate display from original paper Friston 2015
end
 
title('target signal','Fontsize',16)
xlabel('location')
ylabel('location')
set(gca,'Color','k');
axis([-1 1 -1 1]*A*(1+1/16))
axis square, box off
 
 
% free energy and expectations
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
colormap pink
subplot(2,2,1); cla
 
plot(-DEM.J)
title('Free energy','Fontsize',16)
xlabel('time')
ylabel('Free energy')
axis square tight
grid on
 
subplot(2,2,2); cla
v      = spm_unvec(DEM.qU.v{2}(:,end),DEM.M(2).v);
[i j]  = max(v);
v(:,j) = v;
imagesc(spm_softmax(v))
title('softmax expectations','Fontsize',16)
xlabel('cell')
ylabel('cell')
axis square tight
 
 
% target morphology
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf

subplot(2,2,1); cla
Px=-P.x(2,:); % rotate display from original paper Friston 2015
Py=-P.x(1,:); % rotate display from original paper Friston 2015
Ps=P.s;
for i = 1:m
    for j = 1:n
        x = Px(j);
        y = Py(j) + i/6;

        if Ps(i,j)
            plot(x,y,'.','markersize',24,'color','c'); hold on
        else
            plot(x,y,'.','markersize',24,'color','k'); hold on
        end
    end
end
xlabel('cell')
title('Encoding','Fontsize',16)
axis image off
hold off
 
subplot(2,2,2); cla

for i = 1:n
%     x = -P.x(:,i);
    x = -P.x(:,i);    % rotate display from original paper Friston 2015
    c = P.s(end - 2:end,i);
    c = full(max(min(c,1),0));
    plot(x(2),x(1),'.','markersize',16,'color',c);   hold on
    plot(x(2),x(1),'h','markersize',12,'color',h*c); hold on
end
 
title('morphogenesis','Fontsize',16)
xlabel('location')
ylabel('location')
set(gca,'Color','k');
axis([-1 1 -1 1]*A)
axis square, box off
hold off
 
 
% graphics
%--------------------------------------------------------------------------
subplot(2,2,3); cla;
for t = 1:N
    v = spm_unvec(DEM.qU.a{2}(:,t),a);
    for i = 1:n
%         x = -v.x(1,i);
        x = -v.x(1,i); % rotate display from original paper Friston 2015
        c = v.s(end - 2:end,i);
        c = full(max(min(c,1),0));
        plot(t,x,'.','markersize',16,'color',c); hold on
    end
end
 
title('morphogenesis','Fontsize',16)
xlabel('time')
ylabel('location')
set(gca,'Color','k');
set(gca,'YLim',[-1 1]*A)
axis square, box off
hold off
 
% movies
%--------------------------------------------------------------------------
subplot(2,2,4);hold off, cla;
for t = 1:N
    v = spm_unvec(DEM.qU.a{2}(:,t),a);
    
    for i = 1:n
        x = v.x(:,i);
        c = v.s(end - 2:end,i);
        c = max(min(c,1),0);
        plot(x(2),-x(1),'.','markersize',8,'color',full(c)); hold on
        
        % destination
        %------------------------------------------------------------------
        if t == N
%             plot(x(2),x(1),'.','markersize',16,'color',full(c));   hold on
%             plot(x(2),x(1),'h','markersize',12,'color',full(h*c)); hold on
             % rotate display from original paper Friston 2015
            plot(x(2),-x(1),'.','markersize',16,'color',full(c));   hold on
            plot(x(2),-x(1),'h','markersize',12,'color',full(h*c)); hold on
        end
    end
    set(gca,'Color','k');
    axis square, box off
    axis([-1 1 -1 1]*A)
    drawnow
    
    % save
    %----------------------------------------------------------------------
    Mov(t) = getframe(gca);
    
end
 
set(gca,'Userdata',{Mov,8})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('Extrinsic (left click for movie)','FontSize',16)
xlabel('location') 
 
return
 
 
% Equations of motion and observer functions
%==========================================================================
 
% sensed signal
%--------------------------------------------------------------------------
function c = morphogenesis(x,s,y)
% x - location of cells
% s - signals released
% y - location of sampling [default: x]
%__________________________________________________________________________
% preliminaries
%--------------------------------------------------------------------------
if nargin < 3; y = x; end                  % sample locations
n     = size(y,2);                         % number of locations
m     = size(s,1);                         % number of signals
k     = 1;                                 % signal decay over space 
c     = zeros(m,n);                        % signal sensed at each location
for i = 1:n
    for j = 1:size(x,2)
        
        % distance
        %------------------------------------------------------------------
        d      = y(:,i) - x(:,j);
        d      = sqrt(d'*d);
        
        % signal concentration
        %------------------------------------------------------------------

        c(:,i) = c(:,i) + exp(-k*d).*s(:,j);

 
    end
end
 
 
% first level process: generating input
%--------------------------------------------------------------------------
function g = Gg(x,v,a,P)
global t
global double_head
global double_tail
if isempty(t);
    s = 0;
else
    s = (1 - exp(-t*2));
end
a        = spm_unvec(a,P);

if double_head
    g.x(1,:)   = - a.x(1,:).^2 ; 
    g.x(2,:)   = a.x(2,:)  ; 
elseif double_tail    
    g.x(1,:)   = a.x(1,:).^2 ; 
    g.x(2,:)   = a.x(2,:)  ; 
else
    g.x(1,:) = a.x(1,:);                     % position  signal
    g.x(2,:) = a.x(2,:);                     % position  signal
end

g.s      = a.s;                          % intrinsic signal
g.c      = s*morphogenesis(a.x,a.s);     % extrinsic signal
 
% first level model: mapping hidden causes to sensations
%--------------------------------------------------------------------------
function g = Mg(x,v,P)
global t
global tau

if isempty(t);
    s = 0;
else
    s=1-exp(-t/tau);
end

p    = spm_softmax(v);                   % expected identity

g.x  = P.x*p;                            % position
g.s  = P.s*p;                            % intrinsic signal
g.c  = s*P.c*p;                          % extrinsic signal

