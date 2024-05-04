clc;
clear all;
close all;
%% Ue1 Aufgabe 1

% Anfangswert

a = [1.03
    1.54];
Q_a = [5.34 3.84
      3.84 2.80];
% Einfaches Runden
aeinfach = round(a);
%Bootstrapping vor(2.10)
a_Bvor(1) = round(a(1));
a_Bvor(2) = a(2)-Q_a(2,1)/Q_a(1,1)*(a(1)-a_Bvor(1));
a_Bvor(2) = round(a_Bvor(2));
%Bootstrapping back(2.11)
a_Bback(2) = round(a(2));
a_Bback(1) = a(1)-Q_a(1,2)/Q_a(2,2)*(a(2)-a_Bback(2));
a_Bback(1) = round(a_Bback(1));
%Z-Transformation
Q{1} = Q_a;
r{1} = Q{1}(2)/Q{1}(1)*Q{1}(3)
disp(r{1});
for i =1:10  
Z{i} = [-round(Q{i}(2)/Q{i}(1)) 1
       1 0];

Z_product{1} = Z{1};
if i > 1
        
        Z_product {i}= Z{i} * Z_product{i-1};
        
end


 if Z{i}(1) == 0
        break;

    end
Q{i+1} = Z{i}*Q{i}*Z{i}';
r{i+1} = Q{i+1}(2)/Q{i+1}(1)*Q{i+1}(3)
disp(r{i+1});
i =i+1;

end

z = Z_product {i}*a;

% Einfaches Runden
zeinfach = round(z);
%Bootstrapping vor(2.10)
z_Bvor(1) = round(z(1));
z_Bvor(2) = z(2)-Q{i}(2,1)/Q{i}(1,1)*(z(1)-z_Bvor(1));
z_Bvor(2) = round(z_Bvor(2));
%Bootstrapping back(2.11)
z_Bback(2) = round(z(2));
z_Bback(1) = z(1)-Q{i}(1,2)/Q{i}(2,2)*(z(2)-z_Bback(2));
z_Bback(1) = round(z_Bback(1));

% Ergebnisse zusammenstellen
results = [aeinfach a_Bvor' a_Bback' zeinfach z_Bvor' z_Bback'];
% Ergebnisse plotten
figure;
bar(results(1,:), 'grouped');
set(gca, 'XTickLabel', {'a', 'a\_Bvor', 'a\_Bback', 'z', 'z\_Bvor', 'z\_Bback'});
ylabel('Wert der Komponenten');
title('a1 Vergleich der Fixing Methoden');
grid on;
figure;
bar(results(2,:), 'grouped');
set(gca, 'XTickLabel', {'a', 'a\_Bvor', 'a\_Bback', 'z', 'z\_Bvor', 'z\_Bback'});
title('a2 Vergleich der Fixing Methoden');
ylabel('Wert der Komponenten');
grid on;

%% Aufgabe 2
clear all;
% Demo routine for LAMBDA
%
% An example data file will be used with:
%
% a       float ambiguity vector (n x 1)
% Q       variance matrix of float ambiguities (n x n)
%
% OUTPUT
%
% a_ILS   ILS solution (n x 2) (best and second-best candidate vector)
% a_B     Bootstrapping solution (n x 1)
% a_R     Rounding solution (n x 1)
% a_PAR   PAR solution (n x 1) with min.required success rate of 0.995
% a_RT    Solution of ILS with fixed failure rate Ratio Test
% sqnorm  Squared norms of ambiguity residuals (1 x 2) of best and
%         second-best ILS solution
% Ps      Bootstrapping success rate
% PsPAR   Bootstrapping success rate with PAR
% Qzhat   Variance matrix of decorrelated float ambiguities
% Z       Transformation matrix
% nPAR    Number of fixed ambiguities with PAR
% nRT     Number of fixed ambiguities with Ratio Test 
%         (0 if rejected, n if accepted)
% mu      Threshold value used for Ratio Test
%
% Below the 'return' command, more useful examples can be found on how to 
% use the main LAMBDA routine with the different options
%
%------------------------------------------------------------------
% DATE    : 04-MAY-2012                                          
% Author  : Sandra VERHAGEN                                             
%           GNSS Research Centre, Curtin University
%           Mathematical Geodesy and Positioning, Delft University of
%           Technology 
%------------------------------------------------------------------
% Modifiziert für Satellitennavigation Übung 1
% 25.4.2023 D.Becker INS

% Mordifiziert für Satellitennavigation Übung 1
% 28.4.23 Richard Ertmann
%------------------------------------------------------------------
addpath("Lambda_tools\data");
% Specify the file which contains the float solution: a and Q
load amb-ss24   %a-Vektor, Q-Matrix

%   method: 1: ILS method based on search-and-shrink [DEFAULT]
%           2: ILS method based enumeration in search
%           3: integer rounding method
%           4: integer bootstrapping method

%ILS Methode (Method 2)
[a_ILS,sqnorm]                    = LAMBDA(a,Q,1,'ncands',2);
[a_ILS2,sqnorm2]                  = LAMBDA(a,Q,2,'ncands',2);
%integer rounding method (Method 3)
[a_R]                             = LAMBDA(a,Q,3);
%integer bootstrapping method (Method 4)
[a_B]                             = LAMBDA(a,Q,4);

%b) ratio-Test
t0 = sqnorm(1)/sqnorm(2);

% Ergebnisse zusammenstellen
results_ILS = [a a_ILS(:,1) a_ILS2(:,1) a_R a_B];

% Ergebnisse plotten für die erste Komponente
figure;
bar(results_ILS(1,2:5), 'grouped');
set(gca, 'XTickLabel', {'a', 'a_ILS', 'a_ILS2', 'a_R', 'a_B'});
ylabel('Wert der Komponenten');
title('a1 Vergleich der ILS Fixing Methoden');
grid on;

% Ergebnisse plotten für die zweite Komponente
figure;
bar(results_ILS(2,2:5), 'grouped');
set(gca, 'XTickLabel', {'a', 'a_ILS', 'a_ILS2', 'a_R', 'a_B'});
title('a2 Vergleich der ILS Fixing Methoden');
ylabel('Wert der Komponenten');
grid on;

%% Aufgabe 3
clc;
clear all;
close all;

max_size = 45;
for n = 2:max_size   
    tic
    a{n} = rand(n,1)
    randsqmat = rand(n,n);  
    Q{n} = randsqmat*randsqmat';
    [afixed{n},sqnorm{n}] = LAMBDA(a{n},Q{n},2,'ncands',2);
    toc
    compute_times{n} = toc;  
end
compute_times = cell2mat(compute_times(2:max_size));

% 绘制计算时间与向量大小的关系
figure;
plot(2:max_size, compute_times, 'o-');
xlabel('向量大小 Vektorgröße');
ylabel('计算时间 (秒) Time(s)');
title('ILS方法2的计算时间与向量大小的关系 Abhängigkeit der Vektorgröße für die ILS Methode 2');

% 拟合一条曲线
fit_curve = fit((2:max_size)', compute_times', 'smoothingspline');
hold on;
plot(fit_curve);
legend('计算时间 Vektorgröße', '拟合曲线 Polyfit');
%13th Gen Intel(R) Core(TM) i9-13900H   2.60 GHz
%RAM 16.0 GB 
%Matlab R2023a