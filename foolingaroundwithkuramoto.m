% foolingaroundwithkuramoto

fsample = 400;
dt      = 1/fsample;
f       = 40;
we_all      = 1:0.2:10;

% LOAD STRUCTURE
load ~/Downloads/SC_AAL_90.mat    
areas 	= 90;
C     	= SC/max(max(SC))*0.2;
all_delay = 5;

omega = 2*pi*f;
figure;
% INITIALIZATION
for delay = 1 : length(all_delay)
  for wii = 1 : length(we_all)

    we=we_all(wii);
    clear x y amp pha

    pha   = 2*pi*rand(areas,delay+1);
    x     = rand(delay+1,90);
    y     = rand(delay+1,90);
    amp   = rand(delay+1,90);

    tt = delay;

    for t=0:dt:2

      tt = tt+1;

      for i = 1 : 90

        pha(i,tt) = pha(i,tt-1) + dt*(omega + we*sum(C(i,:)'.*sin(pha(i,tt-1)-pha(1:90,tt-delay))));

      end

      x(tt,:) = x(tt-1)+sin(pha(:,tt));
      y(tt,:) = y(tt-1)+cos(pha(:,tt));
      
      plot(x(:,1),y(:,1),'.');
      drawnow
      
      amp(tt,:) = sqrt(x(tt,:).*x(tt,:)+y(tt,:).*y(tt,:));


    end
  % subplot(4,3,we)
  fc_sim = corrcoef(amp);
  fff(delay,wii)=mean(fc_sim(:));

  end
end
  %sqrt(dt)*sig*randn(areas,1)
  
  
  %%
  pha = pha + omega + we*sum(C(n,p)*sin(pha(p)-pha(n)))%sqrt(dt)*sig*randn(areas,1)

  
for t=0:dt:50
  
  sumax = C*x-diag(C*repmat(x',areas,1));
  sumay = C*y-diag(C*repmat(y',areas,1));
  
  x = x+dt*(a.*x-y.*omega-x.*(x.*x+y.*y)+we*sumax)+sqrt(dt)*sig*randn(areas,1);
  y = y+dt*(a.*y+x.*omega-y.*(x.*x+y.*y)+we*sumay)+sqrt(dt)*sig*randn(areas,1);
  
end

% RUN FULL MODEL

fprintf('Running model...\n');
ts_cnt = 0;

for t=0:dt:time_steps  
  
  ts_cnt = ts_cnt+1;
  
  sumax = C*x-diag(C*repmat(x',areas,1));
  sumay = C*y-diag(C*repmat(y',areas,1));

  x = x+dt*(a.*x-y.*omega-x.*(x.*x+y.*y)+(gmod(ts_cnt)+we)*sumax)+sqrt(dt)*sig*randn(areas,1);
  y = y+dt*(a.*y+x.*omega-y.*(x.*x+y.*y)+(gmod(ts_cnt)+we)*sumay)+sqrt(dt)*sig*randn(areas,1);
  
end