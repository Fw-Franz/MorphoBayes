function DEM = DEM_morphogenesis
% This routine illustrates self-assembly or more for genesis under active
% inference (free energy minimisation).  It exploits the fact that one can
% express a systems (marginal) Lyapunov function in terms of a variational
% free energy.  This means that one can prescribe an attracting set in
% terms of the generative model that defines variational free energy.  In
% this example, the attracting set is a point attractor in the phase space
% of a multi-celled organism: where the states correspond to the location
% and (chemotactic) signal expression of 16 cells.  The generative model
% and process are remarkably simple; however, the ensuing migration and
% differentiation of the 16 cells illustrates self-assembly – in the sense
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
% _________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_morphogenesis.m 6290 2014-12-20 22:11:50Z karl $

tic
 
%% preliminaries
%--------------------------------------------------------------------------
clear global
rng('default')
% rng(0,'v5uniform')
% rng(1,'twister');
SPLIT    = 0;                              % split: 1 = upper, 2 = lower, 3=all in random, 4=mid
N        = 32;                             % length of process (bins)
% global tau;  tau = 1/log(2) ;                     % sensitivity s=1-exp(-t/tau)
global tau;  tau = 0.4/log(2) ;                     % sensitivity s=1-exp(-t/tau)
fprintf('final sensitivity coefficient s=1-exp(-1/tau): %f \n',1-exp(-1/tau));
global doublehead; doublehead=0;
global cancer; cancer=1;
global rescue; rescue=0;
global c_n; c_n=1;
flipped=0;
disp_mov_sig=0;
disp_mov_P=0;
belief_1_cell=0;
legend=0;
ni=4; % number of cell to plot beliefs for
movie_path='D:\Box Sync\_PhD\Movies\';

addpath(genpath('D:\Box Sync\_PhD\My_Documents\Matlab_Scripts\'));

% which model fname_subs={'generic'=1; 'planaria_3cells'=2; 'planaria_3cells_packed=3'}
L=1;

% preperations for saving files and other set global parameters
% fname_subs={'generic_'; 'planaria_3cells_'; 'planaria_3cells_packed_'; 'planaria_realistic_';'large_cluster_'};
% fname_sub=fname_subs{L};

fname_sub='linear_pulses';

global split; split=SPLIT;



aa=[ones(1,8) zeros(1,8)];
bb=randperm(16);
cc=aa(bb(1:16))
dd=ones(1,4)
global pulses; pulses=kron(cc,dd)

global pulsed_s; pulsed_s=0; % 0 for none, 1 for cos, 2 for random
global linear; linear = 0;

%% generative process and model
%==========================================================================
M(1).E.d  = 2;                             % approximation order
M(1).E.n  = 2;                             % embedding order
M(1).E.s  = 1;                             % smoothness
 
% priors (prototype)
%--------------------------------------------------------------------------

% template:
% if L == 
%     T =[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% end

 
% if L == 1
%     T =[0 0 0 0 0 0 0 0 0 0 0 0 0;
%         0 0 0 2 0 0 0 0 0 0 0 0 0;
%         0 0 0 0 0 1 0 0 0 0 0 0 0;
%         0 2 0 0 0 0 0 4 0 4 0 3 0;
%         0 0 0 0 0 1 0 0 0 0 0 0 0;
%         0 0 0 2 0 0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 0 0 0 0];
% end


% original
if L == 1
    T =[0 0 2 0 0 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0 0 0;
        2 0 0 0 0 0 4 0 4 0 3;
        0 0 0 0 1 0 0 0 0 0 0;
        0 0 2 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0];
end



% if L == 1
%     T =[0 0 0 0 0 0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 2 0 0 0;
%         0 0 0 0 0 0 0 1 0 0 0 0 0;
%         0 3 0 4 0 4 0 0 0 0 0 2 0;
%         0 0 0 0 0 0 0 1 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 2 0 0 0;
%         0 0 0 0 0 0 0 0 0 0 0 0 0];
% end

% if L == 1
%     T =[0 0 0 0 0 0 0 0 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 2 0 0 0;
%         0 0 0 0 0 0 0 1 0 0 0 0 0;
%         0 4 0 3 0 3 0 0 0 0 0 2 0;
%         0 0 0 0 0 0 0 1 0 0 0 0 0;
%         0 0 0 0 0 0 0 0 0 2 0 0 0;
%         0 0 0 0 0 0 0 0 0 0 0 0 0];
% end

% eliptic morphology, 3 cell types
if L == 2
    T =[0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 2 0 0;
        0 0 3 0 1 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 2 0;
        0 0 3 0 1 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 0 2 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0];
end



% eliptic morphology, 3 cell types, packed
if L == 3
    T =[0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 2 0 4 0 4 0 3 0 0 0;
        0 2 0 0 0 0 0 0 0 0 0 3 0;
        0 0 2 0 4 0 4 0 4 0 3 0 0;
        0 2 0 0 0 0 0 0 0 0 0 3 0;
        0 0 0 2 0 4 0 4 0 3 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0];
    T=rot90(rot90(T));
end

% planaraia realistic morphology
if L == 4
    T =[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 2 0 0 0;
        0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 4 0 0 0 0 0 0 3 0 0 3 0 0 0 0 2 0 0 0 0 0 2;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 2 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
end



if L == 5
    T =[0 0 0 0 0 0 0 0 3 0 0 3 0 0 3 0 0 0 0 0 0 0;
        0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 2 0 0 0 3 0 0 0 3 0 0 0 3 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 2 0 0 0 4 0 0 0 4 0 0 0 4 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 4 0 0 4 0 0 4 0 0 0 0 0 0 0];
end
 
% define logical array corresponding to expression of each of the 4 signals at each given position.     
p(:,:,1) = T > 0;
p(:,:,2) = T == 2 | T == 1;
p(:,:,3) = T == 3 | T == 1;
p(:,:,4) = T == 4;
 
[y x] = find(p(:,:,1));
P.x   = spm_detrend([x(:) y(:)])'/2;


global split_t; split_t=0;
 
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
v     = randn(n,n)/8;                       % states (identity)
g     = Mg([],v,P);
a.x   = g.x;                              % action (chemotaxis)
% a.x   = g.x/10;                              % action (chemotaxis)
a.s   = g.s;                              % action (signal release)
 
 
%% generative process 
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
 
 
%% generative model
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
 
%% solve
%==========================================================================
   
DEM   = spm_ADEM(DEM);
spm_DEM_qU(DEM.qU,DEM.pU);

%% split half simulations
%==========================================================================
if SPLIT
    
    % select (partially diferentiated cells to duplicate
    %----------------------------------------------------------------------
%     t    = N/2;
    t    = N;
    split_t=1;
    v    = spm_unvec(DEM.pU.v{1}(:,t),DEM.M(1).v);
    if SPLIT==2
        [i j] = sort(v.x(1,:), 'ascend');
    elseif SPLIT==1
        [i j] = sort(v.x(1,:),'descend');
    elseif SPLIT==3
        j = randperm(length(v.x(1,:)))
    elseif SPLIT==4
        [i j] = sort(v.x(1,:), 'ascend');
        j = [j(round(n/4)+1:3*round(n/4)) j(8) j(round(n/4)+1:3*round(n/4))];
        j
    end
%     j    = [j(1:n/2) j(1:n/2)];
%     j    = [j(n/4:3*n/4-1) j(n/4:3*n/4-1)];
    
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
    a.x  = a.x(:,j);
    a.s  = randn(size(a.s))/512;
    
    DEM.M(1).v = g;
    DEM.M(2).v = v;
    DEM.G(2).a = spm_vec(a);
    
    % solve
    %----------------------------------------------------------------------
    
    DEM   = spm_ADEM(DEM);
    spm_DEM_qU(DEM.qU,DEM.pU);
end

%% Graphics
%==========================================================================
 
% expected signal concentrations
%--------------------------------------------------------------------------
subplot(2,2,2); cla
A     = max(abs(P.x(:)))*3/2;
h     = 2/3;
 
x     = linspace(-A,A,32);
[x y] = ndgrid(x,x);
x     = spm_detrend([x(:) y(:)])';
c     = morphogenesis(P.x,P.s,x);
c     = c - min(c(:));
c     = c/max(c(:));
for i = 1:size(c,2)
    col = c(end - 2:end,i);
    if flipped        
        plot(x(1,i),x(2,i),'.','markersize',32,'color',col); hold on
    else        
        plot(x(2,i),x(1,i),'.','markersize',32,'color',col); hold on
    end
end

if legend
    title('target signal','Fontsize',16)
    xlabel('location')
    ylabel('location')
end
set(gca,'Color','k');
axis([-1 1 -1 1]*A*(1+1/16))
axis square, box off

for i = 1:n
    x = P.x(:,i);
    c = P.s(end - 2:end,i);
    c = full(max(min(c,1),0));
    
    if flipped
        plot(x(1),x(2),'.','markersize',16,'color',c);   hold on
        plot(x(1),x(2),'h','markersize',12,'color',h*c); hold on
    else
        plot(x(2),x(1),'.','markersize',16,'color',c);   hold on
        plot(x(2),x(1),'h','markersize',12,'color',h*c); hold on
    end
end


fname='target_signal_';
movie_path_full=[movie_path fname fname_sub '.png'];
export_fig(movie_path_full, gcf);
 
 
% free energy and expectations
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
colormap pink
subplot(2,2,1); cla
 
plot(-DEM.J)

if legend
    title('Free energy','Fontsize',16)
    xlabel('time')
    ylabel('Free energy')
end
axis square tight
grid on

fname='Free_E_';
movie_path_full=[movie_path fname fname_sub '.png'];
export_fig(movie_path_full, gcf);

subplot(2,2,2); cla
v      = spm_unvec(DEM.qU.v{2}(:,end),DEM.M(2).v);
[i j]  = max(v);
v(:,j) = v;
imagesc(spm_softmax(v))
if legend
    title('softmax expectations','Fontsize',16)
    xlabel('cell')
    ylabel('cell')
end
axis square tight
 
fname='softmax_';
movie_path_full=[movie_path fname fname_sub '.png'];
export_fig(movie_path_full, gcf);
 
% target morphology
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf

subplot(2,2,1); cla
for i = 1:m
    for j = 1:n
        x = P.x(2,j);
        y = P.x(1,j) + i/6;
        
        if flipped  
            if P.s(i,j)
                plot(y,x,'.','markersize',24,'color','k'); hold on
            else
                plot(y,x,'.','markersize',24,'color','c'); hold on
            end      
        else   
            if P.s(i,j)
                plot(x,y,'.','markersize',24,'color','k'); hold on
            else
                plot(x,y,'.','markersize',24,'color','c'); hold on
            end     
        end
    end
end
xlabel('cell')
% title('Encoding','Fontsize',16)
axis image off
hold off
 
subplot(2,2,2); cla
for i = 1:n
    x = P.x(:,i);
    c = P.s(end - 2:end,i);
    c = full(max(min(c,1),0));
    
    if flipped
        plot(x(1),x(2),'.','markersize',16,'color',c);   hold on
        plot(x(1),x(2),'h','markersize',12,'color',h*c); hold on
    else
        plot(x(2),x(1),'.','markersize',16,'color',c);   hold on
        plot(x(2),x(1),'h','markersize',12,'color',h*c); hold on
    end
end
 
if legend
    title('morphogenesis','Fontsize',16)
    xlabel('location')
    ylabel('location')
end
% set(gca,'Color','k');
axis([-1 1 -1 1]*A)
axis image off, box off

hold off
 
fname='morpho_';
movie_path_full=[movie_path fname fname_sub '.png'];
export_fig(movie_path_full, gcf, '-transparent');

 
% graphics
%--------------------------------------------------------------------------
subplot(2,2,3); cla;
global tt;
for tt = 1:N
    v = spm_unvec(DEM.qU.a{2}(:,tt),a);
    for i = 1:n
        x = v.x(1,i);
        c = v.s(end - 2:end,i);
        c = full(max(min(c,1),0));
        plot(tt,x,'.','markersize',16,'color',c); hold on
    end
end
 
title('morphogenesis','Fontsize',16)
xlabel('time')
ylabel('location')
set(gca,'Color','k');
set(gca,'YLim',[-1 1]*A)
axis square, box off
hold off
 
%% movies
%--------------------------------------------------------------------------

v = spm_unvec(DEM.qU.a{2}(:,N),a);  %   qU.a    = Action
%     v = spm_unvec(DEM.qU.v{2}(:,N),a);   %   qU.v    = Conditional expectation of causal states
A     = max(abs(v.x(:)))*3/2;

%----- signal c plot

if disp_mov_sig
    subplot(2,2,4); hold off, cla;

        h     = 2/3;

        x     = linspace(-A,A,32);
        [x y] = ndgrid(x,x);
        x     = spm_detrend([x(:) y(:)])';
    for t = 1:N
        v = spm_unvec(DEM.qU.a{2}(:,t),a);

        c     = morphogenesis(v.x,v.s,x);
        c     = c - min(c(:));
        c     = c/max(c(:));
        for i = 1:size(c,2)
            col = c(end - 2:end,i);
            if flipped
                plot(x(1,i),x(2,i),'.','markersize',32,'color',col); hold on
            else
                plot(x(2,i),x(1,i),'.','markersize',32,'color',col); hold on                
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
 
    % implay(Mov);

    fname='signal_c_mov_';
    movie_path_full=[movie_path fname fname_sub];

    v = VideoWriter(movie_path_full);
    vid.Quality=100
    % v.CompressionRatio = 3;
    open(v);
    writeVideo(v,Mov);
    close(v);

    % set(gca,'Userdata',{Mov,8})
    % set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
    % title('expected signal (left click for movie)','FontSize',16)
    % xlabel('location') 
end

if disp_mov_P
    subplot(2,2,4); hold off, cla;
    v      = spm_unvec(DEM.qU.v{2}(:,N),DEM.M(2).v);
    [i j]  = max(v);
    v(:,j) = v;
    
    for t = 1:N
        v      = spm_unvec(DEM.qU.v{2}(:,t),DEM.M(2).v);
        [i j]  = max(v);
        v(:,j) = v;
        sm=spm_softmax(v);
        imagesc(sm);
        title('softmax expectations','Fontsize',16)
        xlabel('cell')
        ylabel('cell')
        axis square tight
     
        % save
        %--------------------------------------------------------------
        Mov(t) = getframe(gca);

    end
 
    fname='softmax_mov_';
    if doublehead
        movie_path_full=[movie_path fname 'doublehead_' fname_sub];
    else
        movie_path_full=[movie_path fname fname_sub];
    end
    
    vid = VideoWriter(movie_path_full);
    vid.Quality=100
    % v.CompressionRatio = 3;
    open(vid);
    writeVideo(vid,Mov);
    close(vid);    
end

if belief_1_cell
    subplot(2,2,4); hold off, cla;
    v      = spm_unvec(DEM.qU.v{2}(:,N),DEM.M(2).v);
    [i j]  = max(v);
    v(:,j) = v;
    
    Col=zeros(length(T(:,1)),length(T(1,:)),N);
    
    [yi xi] = find(p(:,:,1));
    
    xmax=length(T(:,1));
    ymax=length(T(1,:));

    x=1:1:xmax;
    y=1:1:ymax;

    xq=1:0.25:xmax;
    yq=1:0.25:ymax;
    
    Vq=zeros(length(xq),length(yq),N);

    [X,Y] = meshgrid(y,x);

    [Xq,Yq] = meshgrid(yq,xq);
    
    
    for t = 1:N
        v      = spm_unvec(DEM.qU.v{2}(:,t),DEM.M(2).v);
        [i j]  = max(v);
        v(:,j) = v;
        sm=spm_softmax(v);     
        
        
        % 3d gaussian smoothed landscape
        %--------------------------------------------------------------
        
        va      = spm_unvec(DEM.qU.a{2}(:,t),a);
        xr = va.x(:,ni);  % real positions of cells:
        
        for i = 1:n
            x = P.x(:,i);
%             col = sm(ni,i);
            col = sm(i,ni);
%             c = full(max(min(c,1),0));
            c = [1-col 1-col 1-col];
%             if flipped
%                 plot(x(1),x(2),'.','markersize',50,'color',c); hold on %full(c)); hold on
%             else
%                 plot(x(2),x(1),'.','markersize',50,'color',c); hold on %full(c)); hold on
%             end
%             axis square, box off
%             axis([-1 1 -1 1]*A)
%             drawnow
            

            Col(yi(i),xi(i),t)=col;
        end

        Vq(:,:,t) = interp2(X,Y,Col(:,:,t),Xq,Yq);

        surfl(Vq(:,:,t)); %hold on

        shading flat
        colormap winter
        zlim([0 0.8])
        axis off

        % save
        %--------------------------------------------------------------
        Mov(t) = getframe(gca);
    end
    

    fname=sprintf('Cell_%i_belief_mov_',ni);
    if doublehead
        movie_path_full=[movie_path fname 'doublehead_' fname_sub];
    else
        movie_path_full=[movie_path fname fname_sub];
    end
    
    vid = VideoWriter(movie_path_full);
    vid.Quality=100
    open(vid);
    writeVideo(vid,Mov);
    close(vid);    
end

%%%%%%%%
subplot(2,2,4);hold off, cla;
for tt = 1:N
    v = spm_unvec(DEM.qU.a{2}(:,tt),a);
    for i = 1:n
        x = v.x(:,i);
        c = v.s(end - 2:end,i);
        c = max(min(c,1),0);
        if flipped
            plot(x(1),x(2),'.','markersize',8,'color',full(c)); hold on
        else
            plot(x(2),x(1),'.','markersize',8,'color',full(c)); hold on
        end
        
        % destination
        %------------------------------------------------------------------
        if tt == N
            if flipped
                plot(x(1),x(2),'.','markersize',16,'color',full(c));   hold on
                plot(x(1),x(2),'h','markersize',12,'color',full(h*c)); hold on
            else
                plot(x(2),x(1),'.','markersize',16,'color',full(c));   hold on
                plot(x(2),x(1),'h','markersize',12,'color',full(h*c)); hold on
            end
        end
    end
    set(gca,'Color','k');
    axis square, box off
    axis([-1 1 -1 1]*A)
    drawnow
    
    % save
    %----------------------------------------------------------------------
    Mov(tt) = getframe(gca);
    
end
 

set(gca,'Userdata',{Mov,8})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('Extrinsic (left click for movie)','FontSize',16)
xlabel('location') 

fname='cell_migration_mov_';
movie_path_full=[movie_path fname fname_sub];

vid = VideoWriter(movie_path_full);
vid.Quality=100
open(vid);
writeVideo(vid,Mov);
close(vid);



toc

return

end


%% Equations of motion and observer functions
%==========================================================================
 
% sensed signal
%--------------------------------------------------------------------------


function c = morphogenesis(x,s,y)

    global rescue
% x - location of cells
% s - signals released
% y - location of sampling [default: x]
%__________________________________________________________________________
 
% preliminaries
%--------------------------------------------------------------------------
if nargin < 3; y = x; end                  % sample locations
n     = size(y,2);                         % number of locations
m     = size(s,1);                         % number of signals
% k     = 0.4;                         
k     = 1;                                 % signal decay over space 
c     = zeros(m,n); 

%disp(y(:,1))

% signal sensed at each location
for i = 1:n
    
    
    for j = 1:size(x,2)
        
        % distance
        %------------------------------------------------------------------
        
        d      = y(:,i) - x(:,j);
       
        d      = sqrt(d'*d);
        
        % signal concentration
        %------------------------------------------------------------------
        if rescue
    %------- rescue of 1-cell misbehavior      
            if x(1,j)<0 && x(2,j)<0;
                c(:,i) = c(:,i) + exp(-k*d^0.5).*s(:,j);
            else
                c(:,i) = c(:,i) + exp(-k*d).*s(:,j);
            end


         %------- alterantive to rescue of 1-cell misbehavior  
    %         if y(1,j)<0 && y(2,j)<0;
    %             c(:,i) = c(:,i) + exp(-k*d^0.5).*s(:,j);
    %         else
    %             c(:,i) = c(:,i) + exp(-k*d).*s(:,j);
    %         end

%         elseif j=1
%             c(:,i) = c(:,i) + exp(-k*d).*s(:,j)+;
        %------- Original                
        else
            c(:,i) = c(:,i) + exp(-k*d).*s(:,j);
        end

    end
end


% % graphic expected signal concentrations
% %--------------------------------------------------------------------------
% subplot(2,2,2); cla
% A     = max(abs(x(:)))*3/2
% h     = 2/3;
%  
% xx     = linspace(-A,A,32);
% [xx yy] = ndgrid(xx,xx);
% xx     = spm_detrend([xx(:) yy(:)])';
% ci     = c - min(c(:));
% ci     = ci/max(ci(:));
% for i = 1:size(ci,2)
%     col = ci(end - 2:end,i);
%     plot(xx(2,i),xx(1,i),'.','markersize',32,'color',col); %hold on
% end
%  
% title('expected signal','Fontsize',16)
% xlabel('location')
% ylabel('location')
% set(gca,'Color','k');
% axis([-1 1 -1 1]*A*(1+1/16))
% axis square, box off

end
    
% first level process: generating input
%--------------------------------------------------------------------------
function g = Gg(x,v,a,P)
    global t
    global tau
    global doublehead
    global cancer
    global split
    global split_t
    global c_n
    global linear
    
    % k     = diag([2 1]);                   % perturbations
    % k     = diag([1 1 1/4 1]);             % perturbations

    a     = spm_unvec(a,P);
    
    if linear 
        s=t/tau; 
    else
        s=1-exp(-t/tau); 
    end

    

    
%------- begin of alternatives for action-response

%       %------- Single Head with reversed polarity
% 
%     g.x(1,:)   = -a.x(1,:) ; 
%     g.x(2,:)   = a.x(2,:)  ;
    
      %------- Double Tail - increase sensitvity to 1.8

%     g.x(1,:)   = a.x(1,:).^2 ; 
%     g.x(2,:)   = a.x(2,:)  ; 
      
      %------- Elongated Double Tail

%     g.x(1,1:4)   = -a.x(1,1:4) ; 
%     g.x(2,1:4)   = a.x(2,1:4) ;  
%     
%     g.x(1,5:8)   = a.x(1,5:8) ; 
%     g.x(2,5:8)   = a.x(2,5:8) ;  

      %------- phenotype

%     g.x(1,:)   = 2*a.x(1,:)  ; 
%     g.x(2,:)   = a.x(2,:)  ;  

      %------- Elongated phenotype

%     g.x(1,:)   = 0.5*a.x(1,:)  ; 
%     g.x(2,:)   = a.x(2,:)  ;  
        
      %------- "cancerous" cell - can bes rescued with alternatives to
      %signal concentration
if cancer 
    n=length(a.x(1,:));
    g.x(1,1:c_n-1)   = a.x(1,1:c_n-1) ; 
    g.x(2,1:c_n-1)   = a.x(2,1:c_n-1)  ;     

    
    g.x(1,c_n)   = a.x(1,c_n).^2 ; 
    g.x(2,c_n)   = a.x(2,c_n)  ; 

    g.x(1,c_n+1:n)   = a.x(1,c_n+1:n) ; 
    g.x(2,c_n+1:n)   = a.x(2,c_n+1:n)  ; 


  %------- Double Head
elseif doublehead
    if split
        if  split_t            
            g.x(1,:)   = -a.x(1,:).^2 ; 
            g.x(2,:)   = a.x(2,:)  ;  
        else            
            g.x(1,:)   = a.x(1,:)  ; 
            g.x(2,:)   = a.x(2,:)  ; 
        end
    else
        g.x(1,:)   = a.x(1,:).^2 ; 
        g.x(2,:)   = a.x(2,:)  ; 
    end
    %------- end of alternatives for action-response (comment out next 2 lines if any if the above statemens outside the if/else conditions were used)
else     
    a.x(1,:)
        g.x(1,:)   = -(a.x(1,:)).^2 ; 
        g.x(2,:)   = a.x(2,:)  ;   
end
    
    g.s   = a.s;                             % intrinsic signal
    g.c   = s*morphogenesis(a.x,a.s);        % extrinsic signal

end
% first level model: mapping hidden causes to sensations
%--------------------------------------------------------------------------
function g = Mg(x,v,P)
    global tt
    global t
    global tau
    global pulses
    global pulsed_s
    global linear

    if isempty(tt); tt = 0; end
    if isempty(t); t = 0; end


    
    if linear
        if pulsed_s==1
            s=cos(tt+pi/2)/4+t/tau;
        elseif pulsed_s==2
            if tt==0
                s=(pulses(tt+1)/2-0.25)+t/tau;
            else
                s=(pulses(tt)/2-0.25)+t/tau;
            end
        else
            s=t/tau;
        end
    end
        
    if pulsed_s==1
        s=cos(tt+pi/2)/4+(1.5-exp(-t/(10/log(2))));
    elseif pulsed_s==2
        if tt==0
            s=(pulses(tt+1)/2-0.25)+(1.5-exp(-t/(10/log(2))));
        else
            s=(pulses(tt)/2-0.25)+(1.5-exp(-t/(10/log(2))));
        end
    else
        s=1-exp(-t/tau);
    end

    
    
    p    = spm_softmax(v);                   % expected identity
    
    for i = 8:-1:1
        if pulsed_s==1
            s_mat(i)=cos(tt+4*i+pi/2)/4+(1.5-exp(-t/(10/log(2))));
        elseif pulsed_s==2
            if tt==0
                s_mat(i)=(pulses(tt+32-4*i+1)/2-0.25)+(1.5-exp(-t/(10/log(2))));
            else
                s_mat(i)=(pulses(tt+32-4*i)/2-0.25)+(1.5-exp(-t/(10/log(2))));
            end
        else
            s_mat(i)=1-exp(-t/tau);
        end
        
    end
    
    g.x  = P.x*p;                            % position
    g.s  = P.s*p;                            % intrinsic signal
    g.c  = s_mat.*P.c*p;                          % extrinsic signal
    
end
