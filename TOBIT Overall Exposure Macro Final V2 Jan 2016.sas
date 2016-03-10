************************************************************************************************** ;
*	ESTIMATION OF OVERALL EXPOSURE EFFECTS FOR TOBIT REGRESSION MODELS                             ;
*   Version 2, Uploaded 01/2016 												               ;
*																					               ;
*	Code written by Wei Wang (2014)													               ;
* 																					               ;
*	Reference:																		               ;
* 	Wang W, and Griswold ME (2015). Natural interpretations in Tobit regression models using       ;
*   marginal estimation methods. Statistical Methods in Medical Research                           ;
*   DOI: 10.1177/0962280215602716                                                                  ; 
* 																					               ;
***************************************************************************************************;
*																					               ;
*	PROGRAM DESCRIPTION:															               ;
*   -----------------------------------------------------------------------------------------------;
*	THIS PROGRAM PRODUCES THE ESTIMATION OF OVERALL EXPSOURE EFFECTS FOR TOBIT                     ; 
*	REGRESSION MODELS. IT CAN ACCOMODATE SINGLE BOUNDARY AND DOUBLE BOUNDARY                       ; 
*   RESPONSE VARIABLES ANA ALLOW MULTIPLE COVARIATES IN THE DATA SET. THIS MACRO                   ;
*   WILL USE TWO MAJOR APPROACHES TO ESTIMATE OVERALL EXPOSURE EFFECTS FOR TOBIT                   ;
*   MODELS, DIRECT-MARGINALIZATION APPROACH AND AVERAGE PREDICTED VALUE APPROACH. THE FORMER       ;
*	PROVIDES TWO ALTERNATIVE METHODS (NEWTON-RAPHSON METHOD OR BRENT'S METHOD) TO                  ;
*   FIND THE SOLUTION FOR THE LATENT DEPENDENT VARIABLE MEAN USING THE NONLINEAR                   ;
*   EQUATION OF MARGINAL MEAN AND LATENT DEPENDENT MEAN IN MODEL ESTIMATION. THE FINAL ESTIMATED   ;
*   OVERALL EXPOSURE EFFECT ESTIMATES FROM NEWTON-RAPHSON BASED OR BRENT'S METHOD BASED            ;
*   DIRECT-MARGINALIZATION APPROACH ARE USUALLY VERY CLOSE. AS WHAT WE MENTIONED IN THE REFERENCE  ;
*   PAPER, THE DIRECT-MARGINALIZATION APPROACH AND AVERAGE PREDICTED VALUE APPROACH ASSUME DIFFERNT;
*   MODEL SETTINGS, IN WHICH THE FORMER ASSUMES HOMOGENUOUS EXPOSURE AND COVARIATE EFFECTS         ;
*   FOR ORIGINAL RESPONSE VARIABLE AND THE LATTER ASSUMES HOMOGENUOUS EXPOSURE AND COVARIATE       ;
*   EFFECTS FOR LATENT DEPENDENT VARIABLE. INFORMATION CRITERION FIT STATISTICS AIC WILL BE        ; 
*   PROVIDED FOR EACH METHOD, AND THE OVERALL EXPOSURE EFFECTS WITH LEAST AIC                      ;
*	IS SUGGESTED TO BE CHOSEN. FOR AVERAGE PREDICTED VALUE APPROACH, EMPIRICAL DISTRIBUTION        ;
*	OF THE BASELINE COVARIATES ARE USED FOR THE ESTIMATION OF OVERALL EXPOSURE                     ;  
*   EFFECTS FOR TARGET POPULATION OR REFERENCE POPULATION. THE EXPOSED GROUP                       ;
*	WAS USED AS TARGET POPULATION IN THIS MACRO.                                                   ;
*																					               ;
*	DATASET REQUIREMENTS:															               ;
*   -----------------------------------------------------------------------------------------------;
*	THE DATASET, MUST HAVE ONE LINE PER SUBJECT WHERE EACH SUBJECT MUST CONTAIN 	               ;
*	ONE BINARY EXPOSURE VARIABLE, ANY NUBMER OF COVARIATES (ZERO TO INFINITY)                      ;
*   AND ONE TRUNCATED NORMALLY DISTRIBUTED OUTCOME WITH INFLATION AT BOUNDARY VALUE.               ;
*   OF NOTE, BASELINE COVARIATES CAN ONLY BE BINARY OR CONTINUOUS, SO CATEGORICAL                  ;
*   BASELINE COVARIATES WITH MORE THAN TWO CATEGORIES SHOULD BE TRANSFORMED                        ;
*   TO DUMMY VARIABLES FIRST BEFORE APPLYINGCAN THIS MACRO. FOR EXAMPLE, CATEGORICAL               ;
*   COVARIATES WITH THREE DISTINCT VALUES SHOULD BE TRANSFORMED TO TWO DUMMY VARIABLES             ;
*   BEFORE THE CURRENT MACRO CAN BE USED.                                                          ;
*	DIRECT-MARGINALIZATION APPROACH USING BRENT'S METHOD NEEDS FROOT FUNCTION IN SAS/IML,          ;
*	WHICH MAY NOT BE AVAILABLE FOR SOME SAS SOFTWARES, SO PLEASE CHOOSE CORRECT SAS MACRO          ;
*	VARIABLE "ANALYSIS" TO EXCLUDE RESULTS FROM THAT METHOD.						               ;
*	THIS MACRO CAN ONLY WORK FOR DATASET WITHOUT MISSING VALUES. THE MISSING VALUES SHOULD BE      ;  
*   IMPUTED OR ANY SUBJECT WITH MISSING DATA SHOULD BE DELETED. ONE SAMPLE                         ;
*	DATA SET WITH TWO COVARIATES AND OUTCOME WITH LOWER MEASUREMENT LIMIT 3.0 IS                   ;
*	LISTED AS FOLLOWING,																	       ;
*																					               ;
*			SUBJECT  EXPOSURE COVARIATE1 COVARIATE2 OUTCOME	            			               ;
*				1		1		21		  1		      13.2        					               ;
*				2		0		28		  0		       3.0          				               ;
*				3		1		25		  1		      12.5	       		    		               ;
*				4		0		24		  1		       3.0	       			    	               ;
*				5		0		26		  1		      10.1	       				                   ;
*				6		0		32		  1		       7.8	       				                   ;
*				7		1		16		  1		       6.5	       				                   ;
*																					               ;
*																					               ;
*	MODEL SPECIFICATION																               ;
*   -----------------------------------------------------------------------------------------------;
*	T: BINARY EXPOSURE Y*: NORMALLY DISTRIBUTED LATENT DEPENDENT VARIALBE WITHOUT                  ;
*   BOUNDARY Y: TRUNCATED NORMALLY DISTRIBUTED OUTCOME WITH INFLATION AT BOUNDAY                   ;
*   VALUE W: BASELINE COVARIATE. Y* CAN ONLY BE OBSERVED WITHIN DETECTABLE RANGE                   ;
*																					               ;
*   REGULAR TOBIT REGRESSION MODEL													               ;
*	Y* = BETA0 + BETA1 * W + BETA2 * T + EPSILON   	WHERE EPSILON~N(0, SIGMASQUARE)                ;
*																					               ;
*   DIRECT-MARGINALIZATION TOBIT REGRESSION MODEL   								               ;
*	E(Y) = GAMMA0 + GAMMA1 * W + GAMMA2 * T                                                        ;
*																					               ;
*													                 				               ;
*	MACRO VARIABLES:																               ;
*   -----------------------------------------------------------------------------------------------;
*	ANALYSIS: 		NUMERIC VALUE 1 AND 2 INDICATING WHETHER RESULTS FROM DIRECT-MARGINALIZATION   ;
*                   APPROACH USING BRENT'S METHOD WILL BE PROVIDED OR NOT                          ;
*																					               ;
*					1 = THREE OVERALL EXPOSURE EFFECTS ESTIMATION                                  ;
*                       DIRECT-MARGINALIZATION APPROACH USING NEWTON-RAPHSON METHOD                ;
*                       DIRECT-MARGINALIZATION APPROACH USING BRENT'S METHOD                       ;
*                       AVERAGE-PREDICTED-VALUE APPROACH                                           ;
*																					               ;
*					2 = TWO OVERALL EXPOSURE EFFECTS ESTIMATION                                    ;
*                       DIRECT-MARGINALIZATION APPROACH USING NEWTON-RAPHSON METHOD                ;
*                       AVERAGE-PREDICTED-VALUE APPROACH                                           ;
*																					               ;
*	DATASET:        MAIN INPUT DATASET												               ;
*																					               ;
*	Y:		SHOULD BE CONTINUOUS VARIABLE WITH OR WITHOUT BOUNDAIES              	               ;
*																					               ;
*	X:		COVARIATE VECTOR INCLUDING BOTH BASELINE COVARIATES AS WELL AS BINARY                  ; 
*           EXPOSURE VARIABLE. OF NOTE, THE MACRO ALLOWS ANY NUMBER OF BASELINE                    ; 
*           COVARIATES, BUT CATEGORICAL COVARIATES WITH MORE THAN TWO DISTINCT                     ; 
*           VALUES SHOULD BE TRANSFORMED TO DUMMY VARIABLES FIRST. AND ALSO BASELINE               ;
*			COVARIATES SHOULD BE SPECIFIED FIRST AND LAST VARIABLE IN COVARIATE                    ;
*           VECTOR X SHOULD BE THE EXPOSURE VARIABLE. 								               ;
*																					               ;
*	CUTOFF1:LOWER DETECTION LIMIT IF APPLICABLE 									               ;
*																					               ;
*	CUTOFF2:UPPER DETECTION LIMIT IF APPLICABLE 									               ;
*																					               ;
*	OF NOTE, IF CUTOFF1 OR CUTOFF2 DOES NOT EXIST (SINGLE BOUNDARY RESPONSE VARIABLE               ;
*   ), ONLY THE AVAILABLE BOUNDARY NEEDS TO BE SPECIFIED, AND THE ONE THAT DOES NOT                ;
*   EXIST CAN BE OMITTED (SEE EXAMPLE BELOW).                                                      ; 
*																					               ;
*	OUT: NAME OF THE OUTPUT DATASET CONTAINING THE OVERALL EXPOSURE EFFECTS ESTIMATE               ;
*																					               ;
*	OF NOTE, WHEN THE MAGNITUDE OF SOME VARIABLES ARE LARGE, THE ESTIMATION METHOD IN THIS MACRO   ;
*	MAY FAIL BY REACHING ABSGCONV CONVERGENCE CRITERION WITH STARTING VALUES. WHEN THIS            ;
*	PROBLEM HAPPENS, ONE METHOD IS TO RESCALE THE INPUT VARIABLES BEFORE APPLYING THE MACRO.       ;
*	FOR EXAMPLE, ONE BASELINE COVARIATE WITH MEAN 100 AND STD 10 CAN BE DIVIDED BY 100 BEFORE      ;
*	APPLYING THE MACRO. WITH THE TRANSFORMED VARIABLE, THE FINAL COEFFICIENT SHOULD BE TRANSFORMED ;
*	BACK TO ORIGINAL SCALE BY DIVIDING 100.      									               ;
*																					               ;
*	EXAMPLE CODE:																	               ;
*																					               ;
*   %include 'TOBITOVERALL.sas' 														           ;
*																					               ;
*	DOUBLE BOUNDARY:    															               ;
*   %TOBITOVERALL(analysis = 1, dataset=DSN, y=RESPONSE, x=COVARIATE1 COVARIATE2 EXPOSURE,         ; 
*   cutoff1 = 0, cutoff2 =100, out =OUT1)                                                         ;
*																					               ;
*	SINGLE BOUNDARY:																               ;
*   %TOBITOVERALL(analysis = 1, dataset=DSN, y=RESPONSE, x=COVARIATE1 COVARIATE2 EXPOSURE,         ; 
*   cutoff1 = 0, out =OUT2)                                                                       ;
*																					               ;
***************************************************************************************************;

options ls=200 ps=64;
dm  'log;clear;out;clear;';

***************************************************************************************************;
*************************************Final Macro***************************************************;
***************************************************************************************************;

%macro tobitoverall(analysis = ., dataset=' ', y=' ', x=' ', cutoff1 =., cutoff2 =., out =' ');

/* General linear model without considering the boundary issue */

%macro tobit0(dsn=, y=,x=, theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
* define var indicating censored ohs;
limit = y;
start LL(theta) global (yx, limit);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do m = 1 to n0;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - (x0[m,]*beta))**2/(2*sigma2);
end;

f = sum1;
return(f);
finish LL ;

* define gradient;
start GRAD(theta) global (yx, limit);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do m = 1 to n0;
   mu = x0[m,]*beta;

   g[1:k1] = g[1:k1] + (y0[m] - mu)/(max(sigma2, exp(-150))) * t(x0[m,]); 
   g[k] = g[k] - 1/sigma + (y0[m] - mu) ** 2/sigma**3;
end;
return(g);
finish GRAD ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
k = ncol(yx) + 1;
n = nrow(y);
theta0 = j(1,k,0);
xx = j(n,1,1) || x;
beta = inv(t(xx)*xx)*t(xx) *y;
e = y - xx*beta;
s = sqrt(ssq(e)/(n- k));
theta0[1:(k-1)] = t(beta);
theta0[k] = s;
print 'OLS starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
       con1 = {  .    . , 
                .    . };
*				theta = {1 1 1 1};
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t3 = t2//t1;
con = t3;
*call nlpqn(rc, theta, 'LL',theta0,optn, con, tc);
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
print rc theta;

yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do m = 1 to n0;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - (x0[m,]*beta))**2/(2*sigma2);
end;

aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

sum1 = theta[k-1];
finalsd = sd[k-1];

postprobs=aic||bic||sum1||finalsd||rc;
 cname = {"NOBOUNDAIC" "NOBOUNDBIC" "NOBOUNDE" "NOBOUNDSD" "NOBOUNDST"};

create ff11 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* APV method: Two Boundary using T = 1 as reference*/

%macro tobit1(dsn=, y=,x=, cutoff1=, cutoff2 =, theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
a = &cutoff1;
b = &cutoff2;
* define var indicating censored ohs;
limit = y;
limit1 = y;
limit2 = y;
limit1 [loc(y=a)] = .;
limit2 [loc(y=b)] = .;
limit3 = limit1;
limit3 [loc(limit1=b)] = .;
n0 = nrow(y);
n1 = ncol(loc( limit1=.));
n2 = ncol(loc( limit2=.));
print n1 "observations are censored at &cutoff1 and" n2 "observations are censored at &cutoff2 out of " n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;
start LL(theta) global (yx, limit, limit1, limit2, limit3, a, b);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit3 ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit1 = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

yx2 = yx[loc(( limit2 = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   chk1 = probnorm((a - x1[i,]*beta)/sigma);
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do j = 1 to n2;
   chk2 = probnorm((x2[j,]*beta - b)/sigma);
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do m = 1 to n0;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - (x0[m,]*beta))**2/(2*sigma2);
end;

f = sum1;
return(f);
finish LL ;

* define gradient;
start GRAD(theta) global (yx, limit, limit1, limit2, limit3, a, b);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit3 ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit1 = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

yx2 = yx[loc(( limit2 = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu = x1[i,]*beta;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   g[1:k1] = g[1:k1] - t4/(max((t2 * sigma), exp(-150))) * t(x1[i,]); 
   g[k] = g[k] + t4 * (mu - a)/(max((sigma2 * t2), exp(-150)));
end;

do j = 1 to n2;
   mu = x2[j,]*beta;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   g[1:k1] = g[1:k1] - tt1/(max((1-t1), exp(-150))) * t(x2[j,]); 
   g[k] = g[k] + t3 * (b - mu)/(max((sigma2 * (1-t1)), exp(-150)));

end;

do m = 1 to n0;
   mu = x0[m,]*beta;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   g[1:k1] = g[1:k1] + (y0[m] - mu)/(max(sigma2, exp(-150))) * t(x0[m,]); 
   g[k] = g[k] - 1/sigma + (y0[m] - mu) ** 2/sigma**3;
end;
return(g);
finish GRAD ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
k = ncol(yx) + 1;
n = nrow(y);
theta0 = j(1,k,0);
xx = j(n,1,1) || x;
beta = inv(t(xx)*xx)*t(xx) *y;
e = y - xx*beta;
s = sqrt(ssq(e)/(n- k));
theta0[1:(k-1)] = t(beta);
theta0[k] = s;
print 'OLS starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
       con1 = {  .    . , 
                .    . };
*				theta = {1 1 1 1};
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t3 = t2//t1;
con = t3;
*call nlpqn(rc, theta, 'LL',theta0,optn, con, tc);
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
print rc theta;

yx0 = yx[loc((limit3 ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit1 = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

yx2 = yx[loc(( limit2 = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   chk1 = probnorm((a - x1[i,]*beta)/sigma);
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do j = 1 to n2;
   chk2 = probnorm((x2[j,]*beta - b)/sigma);
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do m = 1 to n0;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - (x0[m,]*beta))**2/(2*sigma2);
end;

aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
x1 = x[loc((x[,k-2] = 1)),];
n1 = nrow(x1);
if k >= 4 then x0 = j(n1,1,1) || x1[,1:(k-3)]||j(n1,1,0);
else if k = 3 then x0 = j(n1,1,1) ||j(n1,1,0);
x1 = j(n1,1,1) || x1;
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
sum1 = 0;
*print k n1 x1 x0 beta;

do i = 1 to n1;
   mu = x1[i,]*beta;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f11 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b; 
   f12 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   f13 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a - ttt1 * b;

   mu = x0[i,]*beta;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f01 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b; 
   f02 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   f03 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a - ttt1 * b;

   sum1 = sum1 + (f11 - f01);
   g[1:k1] = g[1:k1] + (f12  * t(x1[i,]) - f02 * t(x0[i,])); 
   g[k] = g[k] + (f13 - f03);
end;

sum1 = sum1/n1;
g = g/n1;
finalvar = g * var * t(g);
finalsd = sqrt(finalvar);
print sum1 finalvar finalsd;

postprobs=aic||bic||sum1||finalsd||rc;
 cname = {"APVAIC" "APVBIC" "APVE" "APVSD" "APVST"};

create ff11 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* Direct Marginalization method: two Boundary using Newton Raphson Method*/

%macro tobit2(dsn=, y=,x=, cutoff1=, cutoff2 =, theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
a = &cutoff1;
b = &cutoff2;
* define var indicating censored ohs;
limit = y;
limit1 = y;
limit2 = y;
limit1 [loc(y=a)] = .;
limit2 [loc(y=b)] = .;
limit3 = limit1;
limit3 [loc(limit1=b)] = .;
n0 = nrow(y);
n1 = ncol(loc( limit1=.));
n2 = ncol(loc( limit2=.));
print n1 "observations are censored at &cutoff1 and" n2 "observations are censored at &cutoff2 out of " n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;
start LL(theta) global (yx, limit, limit1, limit2, limit3, a, b);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit3 ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit1 = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

yx2 = yx[loc(( limit2 = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   mu = mu1;
   diff = 2;
   do r = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk1 = probnorm((a - mu)/sigma);
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   diff = 2;
   do s = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk2 = probnorm((mu - b)/sigma);
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   diff = 2;
   do t = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/(max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - mu)**2/(2*sigma2);
end;

f = sum1;
return(f);
finish LL ;
* define gradient;


start GRAD(theta) global (yx, limit, limit1, limit2, limit3, a, b);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit3 ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit1 = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

yx2 = yx[loc(( limit2 = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   mu = mu1;
   diff = 2;
   do r = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/(max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;

   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a - ttt1 * b;

   g[1:k1] = g[1:k1] + sign(f2) * tt2/(max(abs(t2 * f2), exp(-150))) * t(x1[i,]); 
   g[k] = g[k] + t4 * (mu - a)/(max(abs(sigma2 * t2), exp(-150))) - t4 * sign(f2) * ff2/(max(abs(sigma * t2 * f2), exp(-150)));
end;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   diff = 2;
   do s = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/(max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;

   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a - ttt1 * b;

   g[1:k1] = g[1:k1] - sign(f2) * tt1/(max(abs(f2 * (1-t1)), exp(-150))) * t(x2[j,]); 
   g[k] = g[k] + t3 * (b - mu)/(max(abs(sigma2 * (1-t1)), exp(-150))) 
               + t3 * sign(f2) * ff2/(max(abs(f2 * sigma * (1-t1)), exp(-150)));
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   diff = 2;
   do t = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/(max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a - ttt1 * b;

   g[1:k1] = g[1:k1] + sign(f2) * (y0[m] - mu)/(max(abs(f2 * sigma2), exp(-150))) * t(x0[m,]); 
   g[k] = g[k] - 1/sigma + (y0[m] - mu) ** 2/sigma**3 + (y0[m] - mu) * sign(f2)/(max(abs(f2 * sigma2), exp(-150))) * ff2;

end;

return(g);
finish GRAD ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
k = ncol(yx) + 1;
n = nrow(y);
theta0 = j(1,k,0);
xx = j(n,1,1) || x;
beta = inv(t(xx)*xx)*t(xx) *y;
e = y - xx*beta;
s = sqrt(ssq(e)/(n- k));
theta0[1:(k-1)] = t(beta);
theta0[k] = s;
print 'OLS starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
       con1 = {  .    . , 
                .    . };
*				theta = {1 1 1 1};
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t3 = t2//t1;
con = t3;
*call nlpnrr(rc, theta, 'LL',theta0,optn) grd='GRAD';
*call nlpqn(rc, theta, 'LL',theta0,optn) grd='GRAD';
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
*call nlpqn(rc, theta, 'LL',theta0,optn, , tc);
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
print rc;

yx0 = yx[loc((limit3 ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit1 = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

yx2 = yx[loc(( limit2 = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   mu = mu1;
   diff = 2;
   do r = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk1 = probnorm((a - mu)/sigma);
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   diff = 2;
   do s = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk2 = probnorm((mu - b)/sigma);
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   diff = 2;
   do t = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/(max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - mu)**2/(2*sigma2);
end;

aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

sum1 = theta[k-1];
finalsd = sd[k-1];
postprobs=aic||bic||sum1||finalsd||rc;
 cname = {"MARNRAIC" "MARNRBIC" "MARNRE" "MARNRSD" "MARNRST"};

create ff12 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* Direct Marginalization approach: two Boundary using Brent’s numerical root-finding method*/

%macro tobit3(dsn=, y=,x=, cutoff1=, cutoff2 =, theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
a = &cutoff1;
b = &cutoff2;
* define var indicating censored ohs;
limit = y;
limit1 = y;
limit2 = y;
limit1 [loc(y=a)] = .;
limit2 [loc(y=b)] = .;
limit3 = limit1;
limit3 [loc(limit1=b)] = .;
n0 = nrow(y);
n1 = ncol(loc( limit1=.));
n2 = ncol(loc( limit2=.));
print n1 "observations are censored at &cutoff1 and" n2 "observations are censored at &cutoff2 out of " n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;

   start Func(mu) global (mu1, sigma, a, b);
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
   return(f1); 
   /*   return(((probnorm((b - mu)/sigma) - probnorm((a - mu)/sigma)) * mu + (pdf ('normal', (a - mu)/sigma, 0, 1) - pdf ('normal', (b - mu)/sigma, 0, 1)) * sigma */
/*+ probnorm((a - mu)/sigma) * a + (1 - probnorm((b - mu)/sigma)) * b - mu1)); */
   finish;

start LL(theta) global (yx, limit, limit1, limit2, limit3, a, b, mu1, sigma);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit3 ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit1 = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

yx2 = yx[loc(( limit2 = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   /* Found Boundary of mu */
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);   
   chk1 = probnorm((a - mu)/sigma);
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  
   chk2 = probnorm((mu - b)/sigma);
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - mu)**2/(2*sigma2);
end;

f = sum1;
return(f);
finish LL ;
* define gradient;


start GRAD(theta) global (yx, limit, limit1, limit2, limit3, a, b, mu1, sigma);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit3 ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit1 = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

yx2 = yx[loc(( limit2 = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  

   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a - ttt1 * b;

   g[1:k1] = g[1:k1] + sign(f2) * tt2/(max(abs(t2 * f2), exp(-150))) * t(x1[i,]); 
   g[k] = g[k] + t4 * (mu - a)/(max(abs(sigma2 * t2), exp(-150))) - t4 * sign(f2) * ff2/(max(abs(sigma * t2 * f2), exp(-150)));
end;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  

   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a - ttt1 * b;

   g[1:k1] = g[1:k1] - sign(f2) * tt1/(max(abs(f2 * (1-t1)), exp(-150))) * t(x2[j,]); 
   g[k] = g[k] + t3 * (b - mu)/(max(abs(sigma2 * (1-t1)), exp(-150))) 
               + t3 * sign(f2) * ff2/(max(abs(f2 * sigma * (1-t1)), exp(-150)));
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a - ttt1 * b;

   g[1:k1] = g[1:k1] + sign(f2) * (y0[m] - mu)/(max(abs(f2 * sigma2), exp(-150))) * t(x0[m,]); 
   g[k] = g[k] - 1/sigma + (y0[m] - mu) ** 2/sigma**3 + (y0[m] - mu) * sign(f2)/(max(abs(f2 * sigma2), exp(-150))) * ff2;

end;

return(g);
finish GRAD ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
k = ncol(yx) + 1;
n = nrow(y);
theta0 = j(1,k,0);
xx = j(n,1,1) || x;
beta = inv(t(xx)*xx)*t(xx) *y;
e = y - xx*beta;
s = sqrt(ssq(e)/(n- k));
theta0[1:(k-1)] = t(beta);
theta0[k] = s;
print 'OLS starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
       con1 = {  .    . , 
                .    . };
*				theta = {1 1 1 1};
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t3 = t2//t1;
con = t3;
*call nlpnrr(rc, theta, 'LL',theta0,optn) grd='GRAD';
*call nlpqn(rc, theta, 'LL',theta0,optn) grd='GRAD';
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
*call nlpqn(rc, theta, 'LL',theta0,optn, , tc);
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
print rc;

yx0 = yx[loc((limit3 ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit1 = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

yx2 = yx[loc(( limit2 = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   /* Found Boundary of mu */
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);   
   chk1 = probnorm((a - mu)/sigma);
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  
   chk2 = probnorm((mu - b)/sigma);
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = probnorm((a - mu)/sigma);
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = (-1) * t3/sigma;
   tt2 = (-1) * t4/sigma;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = probnorm((a - mu)/sigma);
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - mu)**2/(2*sigma2);
end;

aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

sum1 = theta[k-1];
finalsd = sd[k-1];
postprobs=aic||bic||sum1||finalsd||rc;
 cname = {"MARBRAIC" "MARBRBIC" "MARBRE" "MARBRSD" "MARBRST"};

create ff13 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* APV method: Single Upper Boundary using T = 1 as reference*/

/*********** Tobit Model*********************/

%macro tobit4(dsn=, y=,x=, cutoff2 = , theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
b = &cutoff2;
* define var indicating censored ohs;
limit = y;
limit [loc(y=b)] = .;
n0 = nrow(y);
n2 = ncol(loc( limit=.));
print n2 "observations are censored at &cutoff2 out of " n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;
start LL(theta) global (yx, limit, b);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx2 = yx[loc(( limit = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do j = 1 to n2;
   chk2 = probnorm((x2[j,]*beta - b)/sigma);
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do m = 1 to n0;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - (x0[m,]*beta))**2/(2*sigma2);
end;

f = sum1;
return(f);
finish LL ;

* define gradient;
start GRAD(theta) global (yx, limit, b);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx2 = yx[loc(( limit = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do j = 1 to n2;
   mu = x2[j,]*beta;
   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;

   g[1:k1] = g[1:k1] - tt1/(max((1-t1), exp(-150))) * t(x2[j,]); 
   g[k] = g[k] + t3 * (b - mu)/(max((sigma2 * (1-t1)), exp(-150)));

end;

do m = 1 to n0;
   mu = x0[m,]*beta;
   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;

   g[1:k1] = g[1:k1] + (y0[m] - mu)/(max(sigma2, exp(-150))) * t(x0[m,]); 
   g[k] = g[k] - 1/sigma + (y0[m] - mu) ** 2/sigma**3;
end;
return(g);
finish GRAD ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
k = ncol(yx) + 1;
n = nrow(y);
theta0 = j(1,k,0);
xx = j(n,1,1) || x;
beta = inv(t(xx)*xx)*t(xx) *y;
e = y - xx*beta;
s = sqrt(ssq(e)/(n- k));
theta0[1:(k-1)] = t(beta);
theta0[k] = s;
print 'OLS starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
       con1 = {  .    . , 
                .    . };
*				theta = {1 1 1 1};
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t3 = t2//t1;
con = t3;
*call nlpqn(rc, theta, 'LL',theta0,optn, con, tc);
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
print rc theta;

yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx2 = yx[loc(( limit = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do j = 1 to n2;
   chk2 = probnorm((x2[j,]*beta - b)/sigma);
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do m = 1 to n0;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - (x0[m,]*beta))**2/(2*sigma2);
end;

aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
x1 = x[loc((x[,k-2] = 1)),];
n1 = nrow(x1);
if k >= 4 then x0 = j(n1,1,1) || x1[,1:(k-3)]||j(n1,1,0);
else if k = 3 then x0 = j(n1,1,1) ||j(n1,1,0);
x1 = j(n1,1,1) || x1;
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
sum1 = 0;
*print k n1 x1 x0 beta;

do i = 1 to n1;
   mu = x1[i,]*beta;
   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = 0;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = 0;

   f11 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b; 
   f12 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   f13 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma - ttt1 * b;

   mu = x0[i,]*beta;
   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = 0;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = 0;

   f01 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b; 
   f02 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   f03 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma - ttt1 * b;

   sum1 = sum1 + (f11 - f01);
   g[1:k1] = g[1:k1] + (f12  * t(x1[i,]) - f02 * t(x0[i,])); 
   g[k] = g[k] + (f13 - f03);
end;

sum1 = sum1/n1;
g = g/n1;
finalvar = g * var * t(g);
finalsd = sqrt(finalvar);
print sum1 finalvar finalsd;

postprobs=aic||bic||sum1||finalsd||rc;
 cname = {"APVAIC" "APVBIC" "APVE" "APVSD" "APVST"};

create ff11 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* Direct Marginalization method: Single Upper Boundary using Newton Raphson Method*/

%macro tobit5(dsn=, y=,x=, cutoff2 = , theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
b = &cutoff2;
* define var indicating censored ohs;
limit = y;
limit [loc(y=b)] = .;
n0 = nrow(y);
n2 = ncol(loc( limit=.));
print n2 "observations are censored at &cutoff2 out of " n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;
start LL(theta) global (yx, limit, b);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx2 = yx[loc(( limit = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   diff = 2;
   do s = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk2 = probnorm((mu - b)/sigma);
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   diff = 2;
   do t = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - mu)**2/(2*sigma2);
end;

f = sum1;
return(f);
finish LL ;
* define gradient;

start GRAD(theta) global (yx, limit, b);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx2 = yx[loc(( limit = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   diff = 2;
   do s = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;

   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = 0;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = 0;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma - ttt1 * b;

   g[1:k1] = g[1:k1] - sign(f2) * tt1/(max(abs(f2 * (1-t1)), exp(-150))) * t(x2[j,]); 
   g[k] = g[k] + t3 * (b - mu)/(max(abs(sigma2 * (1-t1)), exp(-150))) 
               + t3 * sign(f2) * ff2/(max(abs(f2 * sigma * (1-t1)), exp(-150)));
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   diff = 2;
   do t = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = 0;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = 0;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma - ttt1 * b;

   g[1:k1] = g[1:k1] + sign(f2) * (y0[m] - mu)/(max(abs(f2 * sigma2), exp(-150))) * t(x0[m,]); 
   g[k] = g[k] - 1/sigma + (y0[m] - mu) ** 2/sigma**3 + (y0[m] - mu) * sign(f2)/(max(abs(f2 * sigma2), exp(-150))) * ff2;

end;

return(g);
finish GRAD ;

optn={1 2};
tc={3000 15000 . 1E-6};
*tc={3000 15000 . 1E-12};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
k = ncol(yx) + 1;
n = nrow(y);
theta0 = j(1,k,0);
xx = j(n,1,1) || x;
beta = inv(t(xx)*xx)*t(xx) *y;
e = y - xx*beta;
s = sqrt(ssq(e)/(n- k));
theta0[1:(k-1)] = t(beta);
theta0[k] = s;
print 'OLS starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
       con1 = {  .    . , 
                .    . };
*				theta = {1 1 1 1};
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t3 = t2//t1;
con = t3;
*call nlpnrr(rc, theta, 'LL',theta0,optn) grd='GRAD';
*call nlpqn(rc, theta, 'LL',theta0,optn) grd='GRAD';
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
*call nlpqn(rc, theta, 'LL',theta0,optn, , tc);
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
print rc;

* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx2 = yx[loc(( limit = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   diff = 2;
   do s = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk2 = probnorm((mu - b)/sigma);
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   diff = 2;
   do t = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - mu)**2/(2*sigma2);
end;

aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

sum1 = theta[k-1];
finalsd = sd[k-1];
postprobs=aic||bic||sum1||finalsd||rc;
 cname = {"MARNRAIC" "MARNRBIC" "MARNRE" "MARNRSD" "MARNRST"};

create ff12 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* Direct Marginalization approach: Single Upper Boundary using Brent’s numerical root-finding method*/

%macro tobit6(dsn=, y=,x=, cutoff2 = , theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
b = &cutoff2;
* define var indicating censored ohs;
limit = y;
limit [loc(y=b)] = .;
n0 = nrow(y);
n2 = ncol(loc( limit=.));
print n2 "observations are censored at &cutoff2 out of " n0 " obs";

* define log likelihood function;
*print limit1, limit2, limit3;

   start Func(mu) global (mu1, sigma, b);

   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 

   return(f1); 
   /*   return(((probnorm((b - mu)/sigma) - probnorm((a - mu)/sigma)) * mu + (pdf ('normal', (a - mu)/sigma, 0, 1) - pdf ('normal', (b - mu)/sigma, 0, 1)) * sigma */
/*+ probnorm((a - mu)/sigma) * a + (1 - probnorm((b - mu)/sigma)) * b - mu1)); */
   finish;

start LL(theta) global (yx, limit, b, mu1, sigma);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx2 = yx[loc(( limit = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;
   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
   mu = froot( "Func", bound);  
   chk2 = probnorm((mu - b)/sigma);
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;
   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - mu)**2/(2*sigma2);
end;

f = sum1;
return(f);
finish LL ;
* define gradient;

start GRAD(theta) global (yx, limit, b, mu1, sigma);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx2 = yx[loc(( limit = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;
   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  

   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = 0;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = 0;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma - ttt1 * b;

   g[1:k1] = g[1:k1] - sign(f2) * tt1/(max(abs(f2 * (1-t1)), exp(-150))) * t(x2[j,]); 
   g[k] = g[k] + t3 * (b - mu)/(max(abs(sigma2 * (1-t1)), exp(-150))) 
               + t3 * sign(f2) * ff2/(max(abs(f2 * sigma * (1-t1)), exp(-150)));
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;
   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  
   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;
   ttt1 = (-1) * (b - mu) * t3/sigma2;
   ttt2 = 0;
   ttt3 = t3 * (b - mu)**2 / sigma**3;
   ttt4 = 0;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma - ttt1 * b;

   g[1:k1] = g[1:k1] + sign(f2) * (y0[m] - mu)/(max(abs(f2 * sigma2), exp(-150))) * t(x0[m,]); 
   g[k] = g[k] - 1/sigma + (y0[m] - mu) ** 2/sigma**3 + (y0[m] - mu) * sign(f2)/(max(abs(f2 * sigma2), exp(-150))) * ff2;

end;

return(g);
finish GRAD ;

optn={1 2};
tc={3000 15000 . 1E-6};
*tc={3000 15000 . 1E-12};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
k = ncol(yx) + 1;
n = nrow(y);
theta0 = j(1,k,0);
xx = j(n,1,1) || x;
beta = inv(t(xx)*xx)*t(xx) *y;
e = y - xx*beta;
s = sqrt(ssq(e)/(n- k));
theta0[1:(k-1)] = t(beta);
theta0[k] = s;
print 'OLS starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
       con1 = {  .    . , 
                .    . };
*				theta = {1 1 1 1};
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t3 = t2//t1;
con = t3;
*call nlpnrr(rc, theta, 'LL',theta0,optn) grd='GRAD';
*call nlpqn(rc, theta, 'LL',theta0,optn) grd='GRAD';
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
*call nlpqn(rc, theta, 'LL',theta0,optn, , tc);
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
print rc;

yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx2 = yx[loc(( limit = .)),];
n2 = nrow(yx2);
y2 = yx2[,1];
x2 = j(n2,1,1) || yx2[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do j = 1 to n2;
   mu1 = x2[j,]*beta;
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;
   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
   mu = froot( "Func", bound);  
   chk2 = probnorm((mu - b)/sigma);
   if chk2 <= 0 then chk2 = 0.1e8;
   sum1 = sum1 + log(chk2);
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   t1 = probnorm((b - mu)/sigma);
   t2 = 0;
   t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
   t4 = 0;
   tt1 = (-1) * t3/sigma;
   tt2 = 0;
   tt3 = t3 * (b - mu) / sigma2;
   tt4 = 0;
   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma - tt1 * b;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
           t1 = probnorm((b - mu)/sigma);
           t2 = 0;
           t3 = pdf ('normal', (b - mu)/sigma, 0, 1);
           t4 = 0;
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + (1 - t1) * b - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - mu)**2/(2*sigma2);
end;

aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

sum1 = theta[k-1];
finalsd = sd[k-1];
postprobs=aic||bic||sum1||finalsd||rc;
 cname = {"MARBRAIC" "MARBRBIC" "MARBRE" "MARBRSD" "MARBRST"};

create ff13 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* APV method: Single Lower Boundary using T = 1 as reference*/

%macro tobit7(dsn=, y=,x=, cutoff1=, theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
a = &cutoff1;
* define var indicating censored ohs;
limit = y;
limit [loc(y=a)] = .;
n0 = nrow(y);
n1 = ncol(loc( limit=.));
print n1 "observations are censored at &cutoff1 out of" n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;
start LL(theta) global (yx, limit, a);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   chk1 = probnorm((a - x1[i,]*beta)/sigma);
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do m = 1 to n0;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - (x0[m,]*beta))**2/(2*sigma2);
end;

f = sum1;
return(f);
finish LL ;

* define gradient;
start GRAD(theta) global (yx, limit, a);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu = x1[i,]*beta;
   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   g[1:k1] = g[1:k1] - t4/(max((t2 * sigma), exp(-150))) * t(x1[i,]); 
   g[k] = g[k] + t4 * (mu - a)/(max((sigma2 * t2), exp(-150)));
end;

do m = 1 to n0;
   mu = x0[m,]*beta;
   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   g[1:k1] = g[1:k1] + (y0[m] - mu)/(max(sigma2, exp(-150))) * t(x0[m,]); 
   g[k] = g[k] - 1/sigma + (y0[m] - mu) ** 2/sigma**3;
end;
return(g);
finish GRAD ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
k = ncol(yx) + 1;
n = nrow(y);
theta0 = j(1,k,0);
xx = j(n,1,1) || x;
beta = inv(t(xx)*xx)*t(xx) *y;
e = y - xx*beta;
s = sqrt(ssq(e)/(n- k));
theta0[1:(k-1)] = t(beta);
theta0[k] = s;
print 'OLS starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
       con1 = {  .    . , 
                .    . };
*				theta = {1 1 1 1};
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t3 = t2//t1;
con = t3;
*call nlpqn(rc, theta, 'LL',theta0,optn, con, tc);
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
print rc theta;

yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   chk1 = probnorm((a - x1[i,]*beta)/sigma);
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do m = 1 to n0;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - (x0[m,]*beta))**2/(2*sigma2);
end;

aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
x1 = x[loc((x[,k-2] = 1)),];
n1 = nrow(x1);
if k >= 4 then x0 = j(n1,1,1) || x1[,1:(k-3)]||j(n1,1,0);
else if k = 3 then x0 = j(n1,1,1) ||j(n1,1,0);
x1 = j(n1,1,1) || x1;
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
sum1 = 0;
*print k n1 x1 x0 beta;

do i = 1 to n1;
   mu = x1[i,]*beta;
   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = 0;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = 0;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f11 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a; 
   f12 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   f13 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a;

   mu = x0[i,]*beta;
   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = 0;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = 0;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f01 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a; 
   f02 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   f03 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a;

   sum1 = sum1 + (f11 - f01);
   g[1:k1] = g[1:k1] + (f12  * t(x1[i,]) - f02 * t(x0[i,])); 
   g[k] = g[k] + (f13 - f03);
end;

sum1 = sum1/n1;
g = g/n1;
finalvar = g * var * t(g);
finalsd = sqrt(finalvar);
print sum1 finalvar finalsd;

postprobs=aic||bic||sum1||finalsd||rc;
 cname = {"APVAIC" "APVBIC" "APVE" "APVSD" "APVST"};

create ff11 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* Direct Marginalization method: Single Lower Boundary using Newton Raphson Method*/

%macro tobit8(dsn=, y=,x=, cutoff1=, theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
a = &cutoff1;
* define var indicating censored ohs;
limit = y;
limit [loc(y=a)] = .;
n0 = nrow(y);
n1 = ncol(loc( limit=.));
print n1 "observations are censored at &cutoff1 out of" n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;
start LL(theta) global (yx, limit, a);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   mu = mu1;
   diff = 2;
   do r = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk1 = probnorm((a - mu)/sigma);
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   diff = 2;
   do t = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - mu)**2/(2*sigma2);
end;

f = sum1;
return(f);
finish LL ;
* define gradient;

* define gradient;
start GRAD(theta) global (yx, limit, a);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   mu = mu1;
   diff = 2;
   do r = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;

   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = 0;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = 0;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a;

   g[1:k1] = g[1:k1] + sign(f2) * tt2/(max(abs(t2 * f2), exp(-150))) * t(x1[i,]); 
   g[k] = g[k] + t4 * (mu - a)/(max(abs(sigma2 * t2), exp(-150))) - t4 * sign(f2) * ff2/(max(abs(sigma * t2 * f2), exp(-150)));
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   diff = 2;
   do t = 1 to 200 while (abs(diff) > 1E-8);
   temp = mu;

   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;

   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = 0;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = 0;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a;

   g[1:k1] = g[1:k1] + sign(f2) * (y0[m] - mu)/(max(abs(f2 * sigma2), exp(-150))) * t(x0[m,]); 
   g[k] = g[k] - 1/sigma + (y0[m] - mu) ** 2/sigma**3 + (y0[m] - mu) * sign(f2)/(max(abs(f2 * sigma2), exp(-150))) * ff2;

end;

return(g);
finish GRAD ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
k = ncol(yx) + 1;
n = nrow(y);
theta0 = j(1,k,0);
xx = j(n,1,1) || x;
beta = inv(t(xx)*xx)*t(xx) *y;
e = y - xx*beta;
s = sqrt(ssq(e)/(n- k));
theta0[1:(k-1)] = t(beta);
theta0[k] = s;
print 'OLS starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
       con1 = {  .    . , 
                .    . };
*				theta = {1 1 1 1};
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t3 = t2//t1;
con = t3;
*call nlpnrr(rc, theta, 'LL',theta0,optn) grd='GRAD';
*call nlpqn(rc, theta, 'LL',theta0,optn) grd='GRAD';
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
*call nlpqn(rc, theta, 'LL',theta0,optn, , tc);
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
print rc;

yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   mu = mu1;
   diff = 2;
   do r = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   chk1 = probnorm((a - mu)/sigma);
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   diff = 2;
   do t = 1 to 200 while (abs(diff) > 1E-8);
*   do while (abs(diff) > 1E-8);
   temp = mu;

   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   mu = mu - sign(f2) * (f1 - mu1)/( max(abs(f2), exp(-150)));

   diff = mu - temp;
   end;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - mu)**2/(2*sigma2);
end;

aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

sum1 = theta[k-1];
finalsd = sd[k-1];
postprobs=aic||bic||sum1||finalsd||rc;
 cname = {"MARNRAIC" "MARNRBIC" "MARNRE" "MARNRSD" "MARNRST"};

create ff12 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

/* Direct Marginalization approach: Single Lower Boundary using Brent’s numerical root-finding method*/

%macro tobit9(dsn=, y=,x=, cutoff1=, theta= 0);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
a = &cutoff1;
* define var indicating censored ohs;
limit = y;
limit [loc(y=a)] = .;
n0 = nrow(y);
n1 = ncol(loc( limit=.));
print n1 "observations are censored at &cutoff1 out of" n0 " obs";
* define log likelihood function;
*print limit1, limit2, limit3;

   start Func(mu) global (mu1, sigma, a);
   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
   return(f1); 
   /*   return(((probnorm((b - mu)/sigma) - probnorm((a - mu)/sigma)) * mu + (pdf ('normal', (a - mu)/sigma, 0, 1) - pdf ('normal', (b - mu)/sigma, 0, 1)) * sigma */
/*+ probnorm((a - mu)/sigma) * a + (1 - probnorm((b - mu)/sigma)) * b - mu1)); */
   finish;

start LL(theta) global (yx, limit, a, mu1, sigma);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   /* Found Boundary of mu */
   mu = mu1;
   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);   
   chk1 = probnorm((a - mu)/sigma);
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - mu)**2/(2*sigma2);
end;

f = sum1;
return(f);
finish LL ;
* define gradient;

* define gradient;
start GRAD(theta) global (yx, limit, a, mu1, sigma);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   mu = mu1;
   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  

   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = 0;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = 0;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a;

   g[1:k1] = g[1:k1] + sign(f2) * tt2/(max(abs(t2 * f2), exp(-150))) * t(x1[i,]); 
   g[k] = g[k] + t4 * (mu - a)/(max(abs(sigma2 * t2), exp(-150))) - t4 * sign(f2) * ff2/(max(abs(sigma * t2 * f2), exp(-150)));
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  
   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;
   ttt1 = 0;
   ttt2 = (-1) * (a - mu) * t4/sigma2;
   ttt3 = 0;
   ttt4 = t4 * (a - mu)**2 / sigma**3;

   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   ff2 = (ttt1 - ttt2) * mu + (t4 - t3) + (ttt4 - ttt3) * sigma + ttt2 * a;

   g[1:k1] = g[1:k1] + sign(f2) * (y0[m] - mu)/(max(abs(f2 * sigma2), exp(-150))) * t(x0[m,]); 
   g[k] = g[k] - 1/sigma + (y0[m] - mu) ** 2/sigma**3 + (y0[m] - mu) * sign(f2)/(max(abs(f2 * sigma2), exp(-150))) * ff2;

end;

return(g);
finish GRAD ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
k = ncol(yx) + 1;
n = nrow(y);
theta0 = j(1,k,0);
xx = j(n,1,1) || x;
beta = inv(t(xx)*xx)*t(xx) *y;
e = y - xx*beta;
s = sqrt(ssq(e)/(n- k));
theta0[1:(k-1)] = t(beta);
theta0[k] = s;
print 'OLS starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
       con1 = {  .    . , 
                .    . };
*				theta = {1 1 1 1};
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t3 = t2//t1;
con = t3;
*call nlpnrr(rc, theta, 'LL',theta0,optn) grd='GRAD';
*call nlpqn(rc, theta, 'LL',theta0,optn) grd='GRAD';
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
*call nlpqn(rc, theta, 'LL',theta0,optn, , tc);
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
print rc;

yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

yx1 = yx[loc(( limit = .)),];
n1 = nrow(yx1);
y1 = yx1[,1];
x1 = j(n1,1,1) || yx1[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma * sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do i = 1 to n1;
   mu1 = x1[i,]*beta;
   /* Found Boundary of mu */
   mu = mu1;
   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);   
   chk1 = probnorm((a - mu)/sigma);
   if chk1 <= 0 then chk1 = 0.1e8;
   sum1 = sum1 + log(chk1);
end;

do m = 1 to n0;
   mu1 = x0[m,]*beta;
   mu = mu1;
   t1 = 1;
   t2 = probnorm((a - mu)/sigma);
   t3 = 0;
   t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
   tt1 = 0;
   tt2 = (-1) * t4/sigma;
   tt3 = 0;
   tt4 = t4 * (a - mu) / sigma2;

   f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
   f2 = (t1 - t2) + (tt1 - tt2) * mu + (tt4 - tt3) * sigma + tt2 * a;
   bound = t(j(2, 1, 1));
   if f1 <= 0 then do;
/*      bound[1] = mu;*/
/*      bound[2] = 300;*/
      if f2 >= 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 + r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[2] = mu;
	  end;

      else if f2 < 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 <= 0);
		   mu = mu1 - r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[1] = mu;
	  end;
   end;

   else if f1 > 0 then do;
      if f2 >= 0 then do;
        bound[2] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 - r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[1] = mu;
	  end;

      else if f2 < 0 then do;
        bound[1] = mu1;
		do r = 1 to 2000 while (f1 > 0);
		   mu = mu1 + r/100 *abs(mu1);
	       t1 = 1;
           t2 = probnorm((a - mu)/sigma);
           t3 = 0;
           t4 = pdf ('normal', (a - mu)/sigma, 0, 1);
           f1 = (t1 - t2) * mu + (t4 - t3) * sigma + t2 * a - mu1; 
        end;
        bound[2] = mu;
	  end;
   end;
/*   else if f1 > 0 then do;*/
/*      bound[1] = -200;*/
/*      bound[2] = mu;*/
/*   end;*/
*bound = {-200 300};
*   print sigma;
   mu = froot( "Func", bound);  
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - mu)**2/(2*sigma2);
end;

aic = sum1*(-2) + 2 * ncol(theta);
bic = ncol(theta) * log(nrow(yx0)) + sum1*(-2);

sum1 = theta[k-1];
finalsd = sd[k-1];
postprobs=aic||bic||sum1||finalsd||rc;
 cname = {"MARBRAIC" "MARBRBIC" "MARBRE" "MARBRSD" "MARBRST"};

create ff13 from postprobs  [ colname=cname ];
append from postprobs;

quit;
run;
%mend;

%if &analysis ne 1 and &analysis ne 2 %then %do;
  proc iml;
    print,"WARNING: NEED TO SPECIFY CORRECT NUMBER, 1 OR 2","PROGRAM WILL TERMINATE",;
  quit;
%end;

%else %if &dataset= | &y= | &x= %then %do;
    proc iml;
      print,"WARNING: NEED TO SPECIFY DATASET, COVARIATE, OUTCOME AND BOTH BOUNDARIES","PROGRAM WILL TERMINATE",;
	quit;
%end;

%else %if &analysis=1 %then %do;

data tt1;
set &dataset;
if &y = &cutoff1 then vis1 = 1; else vis1 = 0;
if &y = &cutoff2 then vis2 = 1; else vis2 = 0;
run;

proc means data = tt1 noprint;
var vis1;
output out = sim1 sum(vis1) = na sum(vis2) = nb;
run;

data _null_;
set sim1;
call symput('NA', na);
call symput('NB', nb);
run;

%if  &NA gt 0  and &NB gt 0 %then %do;

%tobit1(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1, cutoff2 = &cutoff2);
%tobit2(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1, cutoff2 = &cutoff2);
%tobit3(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1, cutoff2 = &cutoff2);

%end;

%else %if  &NA = 0  and &NB gt 0 %then %do;

%tobit4(dsn=&dataset, y=&y, x=&x, cutoff2 = &cutoff2);
%tobit5(dsn=&dataset, y=&y, x=&x, cutoff2 = &cutoff2);
%tobit6(dsn=&dataset, y=&y, x=&x, cutoff2 = &cutoff2);

%end;

%else %if  &NA gt 0  and &NB = 0 %then %do;

%tobit7(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1);
%tobit8(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1);
%tobit9(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1);

%end;

%else %if  &NA = 0  and &NB = 0 %then %do;

%tobit0(dsn=&dataset, y=&y, x=&x);

%end;

%if  &NA = 0  and &NB = 0 %then %do;

data temp;
merge ff11;
if NOBOUNDst = . then do;
NOBOUNDst = NOBOUNDsd;
NOBOUNDsd = .;
end;
if NOBOUNDst lt 0 then do;
NOBOUNDaic = .;
NOBOUNDe = .;
NOBOUNDsd = .;
end;
drop NOBOUNDst;
NOBOUNDCI =trim(left(put((NOBOUNDe - 1.959964 * NOBOUNDsd), 10.4)))||'-'||trim(left(put((NOBOUNDe + 1.959964 * NOBOUNDsd), 10.4))); 
NOBOUNDP = put((2 - 2* probnorm(abs(NOBOUNDe)/NOBOUNDsd)), 5.3);
format method $50.;
Method = 'Regular Linear Regression with Noboundary Data';
run;

data &out;
retain Method NOBOUNDAIC NOBOUNDbic NOBOUNDE NOBOUNDSD NOBOUNDCI NOBOUNDP;
SET TEMP;
RUN;

proc print data = &out noobs label;
label 
NOBOUNDAIC = "AIC"
NOBOUNDbic = "BIC"
NOBOUNDE = "Exposure Effects Estimate"
NOBOUNDSD = "Standard Error of Exposure Effects Estimate"
NOBOUNDCI = "95% Confidence Interval of Exposure Effects Estimate"
NOBOUNDP = "p-Value";
run;

%end;

%else %do;

data temp;
merge ff12 ff13 ff11;
if apvst = . then do;
apvst = apvsd;
apvsd = .;
end;
if marnrst = . then do;
marnrst = marnrsd;
marnrsd = .;
end;
if marbrst = . then do;
marbrst = marbrsd;
marbrsd = .;
end;
if apvst lt 0 then do;
apvaic = .;
apve = .;
apvsd = .;
end;
if marnrst lt 0 then do;
marnraic = .;
marnre = .;
marnrsd = .;
end;
if marbrst lt 0 then do;
marbraic = .;
marbre = .;
marbrsd = .;
end;
drop apvst marnrst marbrst;
run;

data temp1;
set temp;
keep method APVAIC apvbic APVE APVSD CI P;
rename apvaic = aic;
rename apvbic = bic;
rename apve = est;
rename apvsd = std;
CI =trim(left(put((APVe - 1.959964 * APVsd), 10.4)))||'-'||trim(left(put((APVe + 1.959964 * APVsd), 10.4))); 
P = put((2 - 2* probnorm(abs(APVe)/APVsd)), 5.3);
format method $50.;
Method = 'APV Approach';
run;

data temp2;
set temp;
keep method marnrAIC marnrbic marnrE marnrSD CI P;
rename marnraic = aic;
rename marnrbic = bic;
rename marnre = est;
rename marnrsd = std;
CI =trim(left(put((marnre - 1.959964 * marnrsd), 10.4)))||'-'||trim(left(put((marnre + 1.959964 * marnrsd), 10.4))); 
P = put((2 - 2* probnorm(abs(marnre)/marnrsd)), 5.3);
format method $50.;
Method = 'Direct-Marginalization Approach using Newton-Raphson Method';
run;

data temp3;
set temp;
keep method marbrAIC marbrbic marbrE marbrSD CI P;
rename marbraic = aic;
rename marbrbic = bic;
rename marbre = est;
rename marbrsd = std;
CI =trim(left(put((marbre - 1.959964 * marbrsd), 10.4)))||'-'||trim(left(put((marbre + 1.959964 * marbrsd), 10.4))); 
P = put((2 - 2* probnorm(abs(marbre)/marbrsd)), 5.3);
format method $50.;
Method = "Direct-Marginalization Approach using Brent's Method";
run;

data &out;
retain Method AIC BIC est std CI P;
set temp1 temp2 temp3;
run;

proc print data = &out noobs label;
label AIC = "AIC"
BIC = "BIC"
Est = "Overall Exposure Effects Estimate"
StD = "Standard Error of Overall Exposure Effects Estimate"
CI = "95% Confidence Interval of Overall Exposure Effects Estimate"
P = "p-Value";
run;

%end;

%end;

%else %if &analysis=2 %then %do;

data tt1;
set &dataset;
if &y = &cutoff1 then vis1 = 1; else vis1 = 0;
if &y = &cutoff2 then vis2 = 1; else vis2 = 0;
run;

proc means data = tt1 noprint;
var vis1;
output out = sim1 sum(vis1) = na sum(vis2) = nb;
run;

data _null_;
set sim1;
call symput('NA', na);
call symput('NB', nb);
run;

%if  &NA gt 0  and &NB gt 0 %then %do;

%tobit1(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1, cutoff2 = &cutoff2);
%tobit2(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1, cutoff2 = &cutoff2);

%end;

%else %if  &NA = 0  and &NB gt 0 %then %do;

%tobit4(dsn=&dataset, y=&y, x=&x, cutoff2 = &cutoff2);
%tobit5(dsn=&dataset, y=&y, x=&x, cutoff2 = &cutoff2);

%end;

%else %if  &NA gt 0  and &NB = 0 %then %do;

%tobit7(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1);
%tobit8(dsn=&dataset, y=&y, x=&x, cutoff1 = &cutoff1);

%end;

%else %if  &NA = 0  and &NB = 0 %then %do;

%tobit0(dsn=&dataset, y=&y, x=&x);

%end;

%if  &NA = 0  and &NB = 0 %then %do;

data temp;
merge ff11;
if NOBOUNDst = . then do;
NOBOUNDst = NOBOUNDsd;
NOBOUNDsd = .;
end;
if NOBOUNDst lt 0 then do;
NOBOUNDaic = .;
NOBOUNDe = .;
NOBOUNDsd = .;
end;
drop NOBOUNDst;
NOBOUNDCI =trim(left(put((NOBOUNDe - 1.959964 * NOBOUNDsd), 10.4)))||'-'||trim(left(put((NOBOUNDe + 1.959964 * NOBOUNDsd), 10.4))); 
NOBOUNDP = put((2 - 2* probnorm(abs(NOBOUNDe)/NOBOUNDsd)), 5.3);
format method $50.;
Method = 'Regular Linear Regression with Noboundary Data';
run;

data &out;
retain Method NOBOUNDAIC NOBOUNDbic NOBOUNDE NOBOUNDSD NOBOUNDCI NOBOUNDP;
SET TEMP;
RUN;

proc print data = &out noobs label;
label 
NOBOUNDAIC = "AIC"
NOBOUNDbic = "BIC"
NOBOUNDE = "Exposure Effects Estimate"
NOBOUNDSD = "Standard Error of Exposure Effects Estimate"
NOBOUNDCI = "95% Confidence Interval of Exposure Effects Estimate"
NOBOUNDP = "p-Value";
run;

%end;

%else %do;

data temp;
merge ff12 ff11;
if apvst = . then do;
apvst = apvsd;
apvsd = .;
end;
if marnrst = . then do;
marnrst = marnrsd;
marnrsd = .;
end;
if apvst lt 0 then do;
apvaic = .;
apve = .;
apvsd = .;
end;
if marnrst lt 0 then do;
marnraic = .;
marnre = .;
marnrsd = .;
end;
drop apvst marnrst;
run;

data temp1;
set temp;
keep method APVAIC APVBIC APVE APVSD CI P;
rename apvaic = aic;
rename apvbic = bic;
rename apve = est;
rename apvsd = std;
CI =trim(left(put((APVe - 1.959964 * APVsd), 10.4)))||'-'||trim(left(put((APVe + 1.959964 * APVsd), 10.4))); 
P = put((2 - 2* probnorm(abs(APVe)/APVsd)), 5.3);
format method $50.;
Method = 'APV Approach';
run;

data temp2;
set temp;
keep method marnrAIC marnrBIC marnrE marnrSD CI P;
rename marnraic = aic;
rename marnrbic = bic;
rename marnre = est;
rename marnrsd = std;
CI =trim(left(put((marnre - 1.959964 * marnrsd), 10.4)))||'-'||trim(left(put((marnre + 1.959964 * marnrsd), 10.4))); 
P = put((2 - 2* probnorm(abs(marnre)/marnrsd)), 5.3);
format method $50.;
Method = 'Direct-Marginalization Approach using Newton-Raphson Method';
run;

data &out;
retain Method AIC bic est std CI P;
set temp1 temp2;
run;

proc print data = &out noobs label;
label AIC = "AIC"
BIC = "BIC"
Est = "Overall Exposure Effects Estimate"
StD = "Standard Error of Overall Exposure Effects Estimate"
CI = "95% Confidence Interval of Overall Exposure Effects Estimate"
P = "p-Value";
run;

%end;

%end;

%mend;
/**/
/*LIBNAME TE 'C:\Users\wwang\Desktop\SM_Results\TOBIT Exposure Jan 2014\Jan 2014\Continuous Non Zero 200';*/
/**/
/*data sim2;*/
/*set te.sim2;*/
/**if no ge 1 and no le 20 or no ge 101 and no le 120;*/
/*if vis le 60 then vis = 60;*/
/*drop vis1 vis2;*/
/*run;*/
/**/
/*proc export data = sim2 outfile = 'C:\Users\wwang\Box Sync\Archive Done\Tobit Overall Exposure Effect\data\example.xls' dbms = excel replace;*/
/*run;*/

proc import datafile = 'C:\Users\wei wang\Box Sync\Archive Done\Tobit Overall Exposure Effect\data\example.xls' out = sim2 replace;
run;

%tobitoverall(analysis =1, dataset=sim2, y=vis, x=w ztrt, cutoff1 =60, cutoff2 =100, out =out1);
%tobitoverall(analysis =2, dataset=sim2, y=vis, x=w ztrt, cutoff1 =60, cutoff2 =100, out =out2);
%tobitoverall(analysis =2, dataset=sim2, y=vis, x=w ztrt, cutoff1 =60, out =out3);
%tobitoverall(analysis =2, dataset=sim2, y=vis, x=w ztrt, cutoff2 =100, out =out4);
