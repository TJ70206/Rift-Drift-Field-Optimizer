function [nDim, LB, UB, Vio, GloMin, Obj, EvalObj] = ProbInfo(n)
problemDims = [3, 4, 4, 4, 11];
EvalObj = [];
if n == 1
    % Tension/compression spring design
    LB = [0.05,0.25,2.00];
    UB = [2,1.3,15.0];
    Vio = [3 4 0.1 0.1 1];
    Obj = @f1;
    GloMin = 0.012665232788;
elseif n == 2
    % Pressure vessel design
    LB = [0.51,0.51,10,10];
    UB = [99.49,99.49,200,200];
    Vio = [17000 14000 1 1 1];
    Obj = @f2;
    GloMin = 6059.714335048436 ;
elseif n == 3
    % Design of gear train
    LB = 12*ones(1,problemDims(n));
    UB = 60*ones(1,problemDims(n));
    Vio = [1 1];
    Obj = @f3;
    GloMin = 2.70085714e-12-1e-4;
elseif n == 4
    % Corrugated bulkhead design
    LB = [0 0 0 0];
    UB = [100 100 100 5];
    Vio = [10000 10000 500 500 10000 100 1];
    Obj = @f4;
    GloMin = 6.8429580100808;
elseif n == 5
    % Car side impact design
    LB = [0.50 0.50 0.50 0.50 0.50 0.50 0.50 0 0 -30 -30];
    UB = [1.50 1.50 1.50 1.50 1.50 1.50 1.50 1 1 +30 +30];
    Vio = [1 1 1 1 1 1 10 150 2 20 1];
    Obj = @f5;
    GloMin = 22.84296954;
end
if n < 1 || n > numel(problemDims) || ~isfinite(problemDims(n))
    error('ProbInfo:UnsupportedProblem', 'Unsupported engineering problem number: %d', n);
end
if isempty(EvalObj)
    EvalObj = Obj;
end
nDim = problemDims(n);
end

%% Objective and constraint functions
function [z, g, h] = f1(x)
 % Tension/compression spring design
 z = x(:,1).^2.*x(:,2).*(x(:,3)+2);
 h = 0;
 g(:,1) = 1-(x(:,2).^3.*x(:,3))./(71785.*x(:,1).^4);
 g(:,2) = (4.*x(:,2).^2-x(:,1).*x(:,2))./(12566.*(x(:,2).*x(:,1).^3-x(:,1).^4))....
     + 1./(5108.*x(:,1).^2)-1;
 g(:,3) = 1-140.45.*x(:,1)./(x(:,2).^2.*x(:,3));
 g(:,4) = (x(:,1)+x(:,2))./1.5-1;
end

function [z, g, h] = f2(x)
 % Pressure vessel design
 x(:,1) = 0.0625.*round(x(:,1));
 x(:,2) = 0.0625.*round(x(:,2));
 z = 0.6224.*x(:,1).*x(:,3).*x(:,4)+1.7781.*x(:,2).*x(:,3).^2....
     +3.1661.*x(:,1).^2.*x(:,4)+19.84.*x(:,1).^2.*x(:,3);
 g(:,1) = -x(:,1)+0.0193.*x(:,3);
 g(:,2) = -x(:,2)+0.00954.*x(:,3);
 g(:,3) = -pi.*x(:,3).^2.*x(:,4)-4/3.*pi.*x(:,3).^3+1296000;
 g(:,4) = x(:,4)-240;
 h = 0;
end

function [z, g, h] = f3(x)
 % Design of gear train
 x=round(x); term1=1/6.931; term2=(x(3)*x(2))/(x(1)*x(4));
 z = (term1-term2)^2;
 g = 0;
 h = 0;
end

function [z, g, h] = f4(x)
 % Corrugated bulkhead design
 b=x(1); h=x(2); l=x(3); t=x(4); ABD = abs(l^2-h^2);
 z = (5.885*t*(b+l))/(b+(abs(l^2-h^2))^0.5);
 g(1)=-t*h*(0.4*b+l/6)+8.94*(b+(ABD)^0.5);
 g(2)=-t*h^2*(0.2*b+l/12)+2.2*(8.94*(b+(ABD)^0.5))^(4/3);
 g(3)=-t+0.0156*b+0.15;
 g(4)=-t+0.0156*l+0.15;
 g(5)=-t+1.05;
 g(6)=-l+h;
 h = 0;
end

function [z, g, h] = f5(x)
 % Car side impact design
 Sec8 = [0.192 0.345];
 Sec9 = [0.192 0.345];
 nSec8 = numel(Sec8);
 nSec9 = numel(Sec9);
 x(8) = Sec8(min(floor(x(8)*nSec8+1),nSec8));
 x(9) = Sec9(min(floor(x(9)*nSec9+1),nSec9));

 % Objective
 z=1.98+4.90*x(1)+6.67*x(2)+6.98*x(3)+4.01*x(4)+1.78*x(5)+2.73*x(7);

 % Constraints
 Fa =1.16-0.3717*x(2)*x(4)-0.00931*x(2)*x(10)-0.484*x(3)*x(9)+0.01343*x(6)*x(10);
 VCu =0.261-0.0159*x(1)*x(2)-0.188*x(1)*x(8)-0.019*x(2)*x(7)+0.0144*x(3)*x(5)+0.0008757*x(5)*x(10)+0.08045*x(6)*x(9)+0.00139*x(8)*x(11)+0.00001575*x(10)*x(11);
 VCm =0.214+0.00817*x(5)-0.131*x(1)*x(8)-0.0704*x(1)*x(9)+0.03099*x(2)*x(6)-0.018*x(2)*x(7)+0.0208*x(3)*x(8)+0.121*x(3)*x(9)-0.00364*x(5)*x(6)+0.0007715*x(5)*x(10)-0.0005354*x(6)*x(10)+0.00121*x(8)*x(11)+0.00184*x(9)*x(10)-0.02*x(2)^2;
 VCl=0.74-0.61*x(2)-0.163*x(3)*x(8)+0.001232*x(3)*x(10)-0.166*x(7)*x(9)+0.227*x(2)^(2);
 Dur=28.98+3.818*x(3)-4.2*x(1)*x(2)+0.0207*x(5)*x(10)+6.63*x(6)*x(9)-7.7*x(7)*x(8)+0.32*x(9)*x(10);
 Dmr=33.86+2.95*x(3)+0.1792*x(10)-5.057*x(1)*x(2)-11*x(2)*x(8)-0.0215*x(5)*x(10)-9.98*x(7)*x(8)+22*x(8)*x(9);
 Dlr=46.36-9.9*x(2)-12.9*x(1)*x(8)+0.1107*x(3)*x(10);
 Fp=4.72-0.5*x(4)-0.19*x(2)*x(3)-0.0122*x(4)*x(10)+0.009325*x(6)*x(10)+0.000191*x(11)^(2);
 VMBP=10.58-0.674*x(1)*x(2)-1.95*x(2)*x(8)+0.02054*x(3)*x(10)-0.0198*x(4)*x(10)+0.028*x(6)*x(10);
 VFD=16.45-0.489*x(3)*x(7)-0.843*x(5)*x(6)+0.0432*x(9)*x(10)-0.0556*x(9)*x(11)-0.000786*x(11)^(2);
 g = [Fa-1, VCu-0.32, VCm-0.32, VCl-0.32, Dur-32, Dmr-32, Dlr-32, Fp-4, VMBP-9.9, VFD-15.7];
 h = 0;
end