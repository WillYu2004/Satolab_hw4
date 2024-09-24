%Echebarria B, Karma A.
%Mechanisms for initiation of cardiac discordant alternans.
%The European Physical Journal Special Topics 2007; 146: 217‐31.
%http://link.springer.com/article/10.1140/epjst/e2007-00182-y
%

clear all;
close all;
colormap(jet(1024));

set(gcf, 'Renderer', 'zbuffer')

tic

%tissue size
sizex=600;
sizey=600;

% initial values
v=zeros(sizex,sizey);
vnew=zeros(sizex,sizey);
h=0.8*ones(sizex,sizey);
f=0.5*ones(sizex,sizey);
stim=zeros(sizex,sizey);

% constants
tauso=15;
taufi=0.8;
tauh1=4.8;
tauh2=10.0;

tausi=4.0;
tauf1=100;
tauf2=30;

dfu=0.0005;%diffusion constant
dx=0.015;

dt=0.1;

tmax=1000;%total simulation time
tectopic=200;%timing of an ectopic beat
ectsize=200;%size of an ectopic beat

%convert to integer
tnmax=round(tmax/dt);
tnectopic=round(tectopic/dt);

%stimulate the edge
v(:,1:5)=1.0;

%movie frame
frm=1;


% main loop
for tn=0:tnmax

    minf=((v/0.2).^6)./(1+((v/0.2).^6));
    hinf=1./(1+((v/0.1).^6));
    dinf=((v/0.4).^4)./(1+((v/0.4).^4));
    finf=1./(1+((v/0.1).^4));

    tauh=tauh1+tauh2*exp(-20*((v-0.1).^2));
    tauf=tauf2+(tauf1-tauf2).*v.^3;

    jfi=h.*minf.*(v-1.3)/taufi;
    jsi=f.*dinf.*(v-1.4)/tausi;
    jso=(1-exp(-4*v))/tauso;
    ion=-(jfi+jsi+jso-stim);

    % update variables
    v=v+ion*dt;
    h=h+dt*(hinf-h)./tauh;
    f=f+dt*(finf-f)./tauf;


    %non-flux boundary condition
    v(1,:)=v(3,:);
    v(sizex,:)=v(sizex-2,:);
    v(:,1)=v(:,3);
    v(:,sizey)=v(:,sizey-2);
    for cx=2:sizex-1
        for cy=2:sizey-1
            vnew(cx,cy)=v(cx,cy)+(v(cx-1,cy)+v(cx+1,cy)+v(cx,cy+1)+v(cx,cy-1)-4*v(cx,cy))*dfu*dt/(dx*dx)/2;
        end
    end

    vnew(1,:)=vnew(3,:);
    vnew(sizex,:)=vnew(sizex-2,:);
    vnew(:,1)=vnew(:,3);
    vnew(:,sizey)=vnew(:,sizey-2);
    for cx=2:sizex-1
        for cy=2:sizey-1
            v(cx,cy)=vnew(cx,cy)+(vnew(cx-1,cy)+vnew(cx+1,cy)+vnew(cx,cy+1)+vnew(cx,cy-1)-4*vnew(cx,cy))*dfu*dt/(dx*dx)/2;
        end
    end

    if (mod(tn,50)==0)
        v(1,:)=v(3,:);
        v(sizex,:)=v(sizex-2,:);
        v(:,1)=v(:,3);
        v(:,sizey)=v(:,sizey-2);
        surf(v);
        shading interp;
        axis equal;
        axis([-inf inf -inf inf -0.1,2]);
        view(2);
        caxis([0,1.2]);
        M(frm)=getframe;
        savefname=sprintf('v%d.jpg',frm);
        saveas(gcf,savefname)
        frm=frm+1;
    end;

    %ectopic beat
    if (tn==tnectopic)
        for ix=0:sizex
            for iy=0:sizey
                if (ix-sizex/2)*(ix-sizex/2)+(iy-sizey/2)*(iy-sizey/2)<ectsize*ectsize
                    v(ix,iy)=1.0;
                end
            end
        end
    end

end

toc
