function advdiff(ka,R)
m = 0;
x = [0:0.00005:1]; 
t = [0:0.00005:2];
pars=[ka,R];
%    [t,y] = ode15s(@vdp1000,tspan,y0,opts);

interupt_time = 10;
%@(t,sol,flag)interuptFun(t,sol,flag,interupt_time);
opts = odeset('OutputFcn',@interuptFun);

try
%    [t,y] = ode15s(@vdp1000,tspan,y0,opts);
    tic
    sol = pdepe(m,@advdiffpde,@advdiffic,@advdiffbc,x,t,[]);
    save(strcat('logsensing_Ka_smalldeltax=',num2str(ka,'%.5f'),'_R=',num2str(R,'%.5f'),'.mat'),'sol','t','x','ka','R')
catch ME
    if strcmp(ME.identifier,'interuptFun:Interupt')
        disp(ME.message);
        % Do other things
    else
        rethrow(ME); % It's possible the error was due to something else
    end
end

% B = sol(:,:,1);
% C = sol(:,:,2);


function status = interuptFun(t,y,flag,interupt_time)   %#ok<INUSL>
persistent INIT_TIME;
fprintf('%3f', (INIT_TIME));
status = 0;
    switch(flag)
        case 'init'
            INIT_TIME = tic;
            tic
        case 'done'
            clear INIT_TIME;
            1
        otherwise
            elapsed_time = toc(INIT_TIME);
            fprintf('%3f', toc(INIT_TIME));
            if elapsed_time > interupt_time
                clear INIT_TIME;
                str = sprintf('%.6f',elapsed_time);
                error('interuptFun:Interupt',...
                     ['Interupted integration. Elapsed time is ' str ' seconds.']);
            end
    end
end

function [c,f,s] = advdiffpde(x,t,u,DuDx)
    %fprintf('%5d %2d %5d\n',pars(1),pars(2),t)
    %toc
    if toc<6000
        Ki=0.001;
        C=pars(2)*x; %change concentration
        Ka=pars(1); %change Ka

        Vc=DuDx(2)*(Ka-Ki)/((Ki+C)*(Ka+C));

        c = [1; 1];
        f = [1;1].*DuDx+[-Vc; 0].*u; % changed DC DB =1;
        s = [0;0];
    else
        str = sprintf('%.6f',toc);
        error('interuptFun:Interupt',...
                     ['Interupted integration. Elapsed time is ' str ' seconds.']);
    end       
end

function u0 = advdiffic(x)

    
  
a=0.00001;

f=@(x)(1./(sqrt(pi*a))).*exp(-((x-0.5).^2)./(1*a));

 

hh=f(x);
u0 = [hh; pars(2)*x];  %change concentration
end

function [pl,ql,pr,qr] = advdiffbc(xl,ul,xr,ur,t)
    pl = [0; ul(2)];  
    ql = [1; 0]; 
    pr = [0; ur(2)-pars(2)]; %change concentration
    qr = [1; 0]; 
end

end