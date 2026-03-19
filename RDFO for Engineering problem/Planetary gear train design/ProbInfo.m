function [nDim, LB, UB, Vio, GloMin, Obj, EvalObj] = ProbInfo(n)
if nargin < 1
    n = 1;
end
if n ~= 1
    error('ProbInfo:UnsupportedProblem', 'Only Planetary gear train design is available in this submission folder.');
end

nDim = 9;
LB = [17, 14, 14, 17, 14, 48, 1, 1, 1];
UB = [96, 54, 51, 46, 51, 124, 6, 6, 3];
Vio = 1.0e101 * ones(1, 11);
GloMin = NaN;
Obj = @f1;
EvalObj = @f1_eval;
end

%% Planetary gear train objective and constraint functions
function [z, g, h] = f1(x)
 % Planetary gear train design
 [z, g, h] = f1_penalty(x);
end

function [z, g, h] = f1_eval(x)
 [z, g, ~, h] = f1_components(x);
end

function [z, g, h] = f1_penalty(x)
 [z, gRaw, hPenalty, ~] = f1_components(x);
 g = gRaw;
 violatedMask = g > 0;
 g(violatedMask) = g(violatedMask).^2;
 h = hPenalty;
end

function [z, g, hPenalty, hEval] = f1_components(x)
 x = round(reshape(x, 1, []));
 mSet = [1.75, 2.0, 2.25, 2.5, 2.75, 3.0];
 pSet = [3, 4, 5];

 N1 = x(1);
 N2 = x(2);
 N3 = x(3);
 N4 = x(4);
 N5 = x(5);
 N6 = x(6);
 m1 = mSet(min(max(x(7), 1), numel(mSet)));
 m3 = mSet(min(max(x(8), 1), numel(mSet)));
 p = pSet(min(max(x(9), 1), numel(pSet)));

 i1 = N6 / N4;
 i2 = N6 * (N1 * N3 + N2 * N4) / (N1 * N3 * (N6 - N4));
 iR = -(N2 * N6 / (N1 * N3));
 z = max([abs(i1 - 3.11), abs(i2 - 1.84), abs(iR + 3.11)]);

 Dmax = 220;
 delta22 = 0.5;
 delta33 = 0.5;
 delta55 = 0.5;
 delta35 = 0.5;
 delta34 = 0.5;
 delta56 = 0.5;

 betaCos = ((N4 + N5)^2 + (N6 - N3)^2 - (N3 + N5)^2) / (2 * (N6 - N3) * (N4 + N5));
 betaCos = min(max(betaCos, -1), 1);
 beta = acos(betaCos);

 g(1) = m3 * (N6 + 2.5) - Dmax;
 g(2) = m1 * (N1 + N2) + m1 * (N2 + 2) - Dmax;
 g(3) = m3 * (N4 + N5) + m3 * (N5 + 2) - Dmax;
 g(4) = abs(m1 * (N1 + N2) - m3 * (N6 - N3)) - m1 - m3;
 g(5) = -(N1 + N2) * sin(pi / p) + N2 + 2 + delta22;
 g(6) = -(N6 - N3) * sin(pi / p) + N3 + 2 + delta33;
 g(7) = -(N4 + N5) * sin(pi / p) + N5 + 2 + delta55;
 g(8) = (N3 + N5 + 2 + delta35)^2 - (N6 - N3)^2 - (N4 + N5)^2 + 2 * (N6 - N3) * (N4 + N5) * cos(2 * pi / p - beta);
 g(9) = N4 - N6 + 2 * N5 + 2 * delta56 + 4;
 g(10) = 2 * N3 - N6 + N4 + 2 * delta34 + 4;
 hEval = rem(N6 - N4, p);
 hPenalty = ((N6 - N4) / p)^2 * double(hEval ~= 0);
end