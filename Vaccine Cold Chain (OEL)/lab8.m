%% lab8.m  -  Vaccine cold-chain temperature excursion  (Numerical Methods, Open Ended Lab 8)
% Run with:   lab8
% Opens a tabbed window (MATLAB) and prints the same report in the Command Window.
% Only base MATLAB is needed (no toolbox): the natural cubic spline is coded by hand.
%
% ASSUMPTIONS (change them in section 0 if your instructor says otherwise)
%   * Origin / nodes: t = 0,2,...,12 h, equally spaced, h = 2 h.  Origins are chosen per task below.
%   * "Out of specification" means T > 8.0 deg C (strictly above).
%   * The lab gives no allowed excursion time.  MAX_ALLOWED_H = 0 means "any time above 8 deg C
%     means the shipment is not released automatically" (usual rule for 2-8 deg C vaccines).
%   * Degrees used: full table = degree 6 (7 readings); Bessel = 5; backward in the validation = 4.
clear; clc; close all;

%% 0. DATA AND SETTINGS
t = 0:2:12;                                  % time (h)
T = [5.2 6.4 7.9 9.5 10.1 9.3 7.6];          % logged product temperature (deg C)
n = numel(t);
LIMIT         = 8.0;                         % deg C
MAX_ALLOWED_H = 0;                           % assumed allowable time above LIMIT (h)
DEG_STIR = 6;  DEG_BESS = 5;  DEG_BACK = 4;  % degrees used in the validation at t = 6.4 h
W = @(pre,s) wrapText(pre,s,104);            % line wrapper for the report

%% 1. DIFFERENCE TABLE
h  = diff(t);
equalSp = max(abs(h - h(1))) < 1e-12;
hh = h(1);
D  = buildDiff(T);                           % D(i,k+1) = Delta^k y_i
L1 = {};
L1{end+1} = 'DIFFERENCE TABLE (forward differences, 4 d.p.)';
hdr = sprintf('%6s %9s','t (h)','T');
for k = 1:6, hdr = [hdr sprintf(' %9s',sprintf('Delta^%d',k))]; end %#ok<AGROW>
L1{end+1} = hdr;
tabCell = cell(n,n+1);
for i = 1:n
    s = sprintf('%6.1f %9.4f',t(i),D(i,1));
    tabCell{i,1} = sprintf('%.1f',t(i));
    for k = 1:n
        if i <= n-k+1
            tabCell{i,k+1} = sprintf('%.4f',D(i,k));
        else
            tabCell{i,k+1} = '';
        end
    end
    for k = 2:n
        if i <= n-k+1, s = [s sprintf(' %9.4f',D(i,k))]; else, s = [s sprintf(' %9s','')]; end %#ok<AGROW>
    end
    L1{end+1} = s; %#ok<SAGROW>
end
L1{end+1} = '';
if equalSp
    L1 = [L1, W('', sprintf('Spacing: the readings ARE equally spaced, h = %.4f h (largest deviation of any step from h = %.1e). Forward, backward and central (Stirling, Bessel) difference formulas are all allowed.',hh,max(abs(h-hh))))];
else
    L1 = [L1, W('', 'Spacing: NOT equally spaced - only divided differences / Lagrange are allowed.')];
end
L1{end+1} = 'Size of each difference order (min .. max of the column):';
for k = 2:n
    col = D(1:n-k+1,k);
    L1{end+1} = sprintf('   Delta^%d : %9.4f .. %9.4f   (spread %.4f)',k-1,min(col),max(col),max(col)-min(col)); %#ok<SAGROW>
end
L1 = [L1, W('', ['Degree supported: no difference column is constant or close to zero - the sixth difference (' sprintf('%.4f',D(1,7)) ') is as large as the lower orders. ' ...
    'The table therefore gives no evidence that a polynomial of degree below 6 reproduces the readings, so the only polynomial through all seven points is the degree-6 one (n-1 = 6). ' ...
    'Because the high differences do not shrink (each order amplifies the 0.1 deg C logger resolution), the high-order terms are used with caution and every estimate is cross-checked.'])];

%% shared tools: degree-6 polynomial (centred on t = 6 h) and natural cubic spline
tc0 = t(4);                                  % centre of the table
cp  = polyfit((t-tc0)/hh, T, n-1);           % exact interpolating polynomial (scaled variable)
P   = @(x) polyval(cp,(x-tc0)/hh);
dP  = polyder(cp);
S   = nsCoef(t,T);                           % natural cubic spline

%% 2. ESTIMATES AT t = 5.0, 6.4, 11.0 h
est(1) = struct('x',5.0, 'name','Bessel',          'm',3, 'deg',5,        'nodes',1:6);
est(2) = struct('x',6.4, 'name','Stirling',        'm',4, 'deg',6,        'nodes',1:7);
est(3) = struct('x',11.0,'name','Newton backward', 'm',n, 'deg',6,        'nodes',1:7);
L2 = {};
L2{end+1} = 'ESTIMATES';
estVal = zeros(1,3);
for i = 1:3
    e = est(i);  x = e.x;  p = (x - t(e.m))/hh;
    switch e.name
        case 'Bessel',          [v,lab,tv] = bessF(D,e.m,p,e.deg);
        case 'Stirling',        [v,lab,tv] = stirF(D,e.m,p,e.deg);
        otherwise,              [v,lab,tv] = backF(D,p,e.deg);
    end
    estVal(i) = v;
    vLag  = lagr(t(e.nodes),T(e.nodes),x);   % second method: Lagrange on the same nodes
    vFull = P(x);                            % degree-6 polynomial through all readings
    vSpl  = nsEval(S,x);                     % natural cubic spline
    L2{end+1} = ''; %#ok<SAGROW>
    L2{end+1} = sprintf('--- Estimate %d :  T(%.4f h) ---',i,x); %#ok<SAGROW>
    L2{end+1} = sprintf('Formula : %s,  origin t0 = %.1f h (row %d),  p = (x - t0)/h = %.4f,  degree %d,  nodes t = %g .. %g h', ...
                        e.name,t(e.m),e.m,p,e.deg,t(e.nodes(1)),t(e.nodes(end))); %#ok<SAGROW>
    switch i
        case 1
            why = sprintf(['The point lies in the middle of the interior interval 4-6 h (p = %.4f). Bessel is the central formula built for the middle of an interval: ' ...
                'at p = 1/2 every term with the factor (p - 1/2) vanishes, so the odd differences drop out. Stirling needs |p| <= 0.25 and forward/backward formulas suit only the ends of the table.'],p);
        case 2
            why = sprintf(['The point lies %.1f h from the node t = 6 h at the centre of the table, p = %.4f, inside the Stirling range |p| <= 0.25. The node set 0-12 h is symmetric about t0, ' ...
                'so Stirling can use all seven readings (degree 6).'],x-t(e.m),p);
        otherwise
            why = sprintf(['The point is in the last interval, p = %.4f from the final node. There are no readings to the right, so central formulas cannot be centred here; ' ...
                'Newton backward starts from the last row (the best-supported differences) and uses all seven readings. This is the least trustworthy region of the table.'],p);
    end
    L2 = [L2, W('Why     : ',why)]; %#ok<AGROW>
    L2{end+1} = 'Terms   :'; %#ok<SAGROW>
    for q = 1:numel(tv)
        L2{end+1} = sprintf('   term %d  %-34s = %10.4f',q-1,lab{q},tv(q)); %#ok<SAGROW>
    end
    L2{end+1} = sprintf('Result  : T(%.4f h) = %.4f deg C',x,v); %#ok<SAGROW>
    L2{end+1} = sprintf('Check   : Lagrange, same nodes      = %.4f   (difference %.1e)',vLag,abs(v-vLag)); %#ok<SAGROW>
    L2{end+1} = sprintf('          degree-6 polynomial (all) = %.4f   (difference %.4f)',vFull,abs(v-vFull)); %#ok<SAGROW>
    L2{end+1} = sprintf('          natural cubic spline      = %.4f   (difference %.4f)',vSpl,abs(v-vSpl)); %#ok<SAGROW>
    acc = max([abs(v-vFull),abs(v-vSpl)]);
    L2 = [L2, W('Accuracy: ',sprintf(['The two formulas agree to rounding (same polynomial), so the arithmetic is right. The independent models differ by up to %.4f deg C; ' ...
        'that is the realistic uncertainty of the estimate. The logger resolves only 0.1 deg C, so the 4th decimal is shown only because the lab asks for it.'],acc))]; %#ok<AGROW>
    chk(i).vLag = vLag; chk(i).vFull = vFull; chk(i).vSpl = vSpl; %#ok<SAGROW>
end

%% 3. MAXIMUM
rc = roots(dP);  rc = rc(abs(imag(rc)) < 1e-9);
tcrit = sort(tc0 + hh*real(rc));  tcrit = tcrit(tcrit >= t(1) & tcrit <= t(n));
[Tlog,kmax] = max(T);
lo = t(max(kmax-1,1));  hi = t(min(kmax+1,n));
cand = [tcrit(tcrit>=lo & tcrit<=hi); lo; hi];
[Pm,ix] = max(P(cand));  tmax = cand(ix);
gg = linspace(lo,hi,40001);  [gm,gi] = max(P(gg));
[tsm,fsm] = fminbnd(@(x) -nsEval(S,x), lo, hi, optimset('TolX',1e-10));  Tsm = -fsm;
L3 = {};
L3{end+1} = 'MAXIMUM TEMPERATURE';
L3 = [L3, W('',sprintf('Degree used: %d (the interpolating polynomial through all %d readings, written in s = (t - %g)/h and checked against the Stirling value at 6.4 h: difference %.1e).',n-1,n,tc0,abs(P(6.4)-estVal(2))))];
L3{end+1} = sprintf('Critical points of dT/dt = 0 inside [0,12] h : %s', sprintf('%.4f ',tcrit));
L3{end+1} = sprintf('Highest logged reading                    : %.4f deg C at t = %.4f h', Tlog, t(kmax));
L3{end+1} = sprintf('Interpolated maximum                      : T = %.4f deg C at t = %.4f h', Pm, tmax);
L3{end+1} = sprintf('Check on the adjacent nodes [%g h, %g h]     : P(%g) = %.4f, P(%g) = %.4f, dense scan (40001 points) maximum = %.4f at t = %.4f h', lo,hi,lo,P(lo),hi,P(hi),gm,gg(gi));
L3{end+1} = sprintf('   one critical point in the interval: %d ; largest scanned value minus Pm = %.1e  ->  no larger value between the nodes.', sum(tcrit>=lo & tcrit<=hi), gm-Pm);
L3{end+1} = sprintf('Second method (natural cubic spline)      : T = %.4f deg C at t = %.4f h', Tsm, tsm);
L3{end+1} = sprintf('Difference between the two methods       : %.4f deg C in T, %.4f h in time', abs(Pm-Tsm), abs(tmax-tsm));
L3 = [L3, W('Accuracy: ',sprintf(['The peak lies between the readings, so it is an estimate: it exceeds the highest reading by only %.4f deg C and the two methods agree within %.4f deg C. ' ...
    'The peak is on a flat top, so its time is far less certain (%.4f h between methods) than its height.'],Pm-Tlog,abs(Pm-Tsm),abs(tmax-tsm)))];

%% 4. CROSSINGS OF 8.0 deg C AND RELEASE DECISION
cc = cp;  cc(end) = cc(end) - LIMIT;
r  = roots(cc);  r = r(abs(imag(r)) < 1e-9);
tx = sort(tc0 + hh*real(r));  tx = tx(tx >= t(1) & tx <= t(n));
g  = linspace(t(1),t(n),12001);  sg = nsEval(S,g) - LIMIT;
ib = find(sg(1:end-1).*sg(2:end) < 0);
txs = zeros(numel(ib),1);
for q = 1:numel(ib), txs(q) = fzero(@(x) nsEval(S,x)-LIMIT, [g(ib(q)) g(ib(q)+1)]); end
tabove = timeAbove(P,tx,t(1),t(n),LIMIT);
tabovs = timeAbove(@(x) nsEval(S,x),txs,t(1),t(n),LIMIT);
L4 = {};
L4{end+1} = sprintf('CROSSINGS OF %.1f deg C AND RELEASE DECISION',LIMIT);
L4{end+1} = 'Method 1: roots of P(t) - 8 = 0 for the degree-6 polynomial; Method 2: root of the natural cubic spline minus 8.';
for q = 1:numel(tx)
    if polyval(dP,(tx(q)-tc0)/hh) > 0, dirn = 'rising (excursion begins)'; else, dirn = 'falling (excursion ends)'; end
    L4{end+1} = sprintf('   crossing %d : polynomial t = %.4f h  |  spline t = %.4f h  |  difference %.4f h   %s', q,tx(q),txs(q),abs(tx(q)-txs(q)),dirn); %#ok<SAGROW>
end
L4{end+1} = sprintf('Time above %.1f deg C :  polynomial %.4f h  |  spline %.4f h  |  difference %.4f h', LIMIT,tabove,tabovs,abs(tabove-tabovs));
ia = find(T > LIMIT);
L4 = [L4, W('',sprintf('Bracket from the raw readings alone: readings at t = %g..%g h are above the limit and the neighbours (%g h, %g h) are below, so the true time above lies between %.4f h and %.4f h. The interpolated value falls inside that bracket.', ...
    t(ia(1)),t(ia(end)),t(ia(1)-1),t(ia(end)+1),t(ia(end))-t(ia(1)),t(ia(end)+1)-t(ia(1)-1)))];
% inverse-interpolation check (reversed variables) - NOT used for the answer
L4{end+1} = 'Inverse-interpolation check (t as a polynomial in T on the nodes nearest each crossing) - for comparison only:';
for q = 1:numel(tx)
    k = find(t <= tx(q),1,'last');  k = min(k,n-1);
    ids = max(k-1,1):min(k+2,n);
    if ~(all(diff(T(ids))>0) || all(diff(T(ids))<0)), ids = k:k+1; end
    tinv = lagr(T(ids),t(ids),LIMIT);
    L4{end+1} = sprintf('   crossing %d : nodes t = %s h -> reversed polynomial gives %.4f h,  direct interpolant gives %.4f h,  difference %.4f h', ...
        q,mat2str(t(ids)),tinv,tx(q),abs(tinv-tx(q))); %#ok<SAGROW>
end
L4 = [L4, W('',['The reversed polynomial is not the inverse function and is badly conditioned where the data are flat (near the peak), so the crossing times quoted above are the roots of the DIRECT interpolant, with the reversed result used only as a check.'])];
rel = tabove <= MAX_ALLOWED_H;
if rel, dec = 'RELEASE'; else, dec = 'DO NOT RELEASE - hold in quarantine and refer to QA / the manufacturer for a stability assessment'; end
L4{end+1} = '';
L4{end+1} = sprintf('DECISION : %s',dec);
L4 = [L4, W('Basis    : ',sprintf(['The pallet spent %.4f h (about %.1f h) above %.1f deg C, peaking at %.4f deg C, against an assumed allowance of %.1f h. ' ...
    'The decision does not hinge on the interpolation details: the two methods differ by %.4f h in duration and even the raw readings alone show at least %.1f h above the limit. ' ...
    'MAX_ALLOWED_H is an assumption (the lab sets no number) and can be edited at the top of the script.'],tabove,tabove,LIMIT,Pm,MAX_ALLOWED_H,abs(tabove-tabovs),t(ia(end))-t(ia(1))))];

%% 5. VALIDATION AT t = 6.4 h
xv = 6.4;
mS = 4;  pS = (xv - t(mS))/hh;                         % central formulas: origin t = 6 h
[vSt,~,~]  = stirF(D,mS,pS,DEG_STIR);
[vBe,~,~]  = bessF(D,mS,pS,DEG_BESS);
pB = (xv - t(n))/hh;                                   % backward formula: origin t = 12 h
[vBa,~,~]  = backF(D,pB,DEG_BACK);
vSp = nsEval(S,xv);
vals  = [vSt vBe vBa vSp];
names = {'Stirling','Bessel','Backward','Natural cubic spline'};
setup = {sprintf('origin 6 h, p = %.4f, degree %d (t = 0-12 h)',pS,DEG_STIR), ...
         sprintf('origin 6 h, p = %.4f, degree %d (t = 2-12 h)',pS,DEG_BESS), ...
         sprintf('origin 12 h, p = %.4f, degree %d (t = %g-12 h)',pB,DEG_BACK,t(n-DEG_BACK)), ...
         'cubic pieces, zero 2nd derivative at both ends'};
[vmax,imx] = max(vals);  [vmin,imn] = min(vals);
valTab = cell(4,4);
for q = 1:4
    valTab(q,:) = {names{q}, setup{q}, sprintf('%.4f',vals(q)), sprintf('%.4f',vals(q)-vmin)};
end
% how the backward formula behaves as terms are added
L5 = {};
L5{end+1} = sprintf('VALIDATION AT t = %.1f h',xv);
L5{end+1} = sprintf('%-22s %-48s %10s %12s','Method','Settings','T (deg C)','minus lowest');
for q = 1:4, L5{end+1} = sprintf('%-22s %-48s %10.4f %12.4f',names{q},setup{q},vals(q),vals(q)-vmin); end %#ok<SAGROW>
L5{end+1} = '';
L5{end+1} = sprintf('Spread (highest - lowest) = %.4f deg C   (highest: %s, lowest: %s)',vmax-vmin,names{imx},names{imn});
L5{end+1} = 'Backward formula as terms are added (origin 12 h, p = -2.8000):';
for dg = 1:n-1
    vb = backF(D,pB,dg);
    L5{end+1} = sprintf('   degree %d : %.4f   (%.4f from the Stirling value)',dg,vb,abs(vb-vSt)); %#ok<SAGROW>
end
% spread of the four curves over the whole interval
gx = 0:0.05:12;
Cst = arrayfun(@(x) stirF(D,mS,(x-t(mS))/hh,DEG_STIR),gx);
Cbe = arrayfun(@(x) bessF(D,mS,(x-t(mS))/hh,DEG_BESS),gx);
Cba = arrayfun(@(x) backF(D,(x-t(n))/hh,DEG_BACK),gx);
Csp = nsEval(S,gx);
CC  = [Cst; Cbe; Cba; Csp];
spread = max(CC,[],1) - min(CC,[],1);
r1 = gx <= 4;  r2 = gx >= 4;
[s1,i1] = max(spread(r1));  g1 = gx(r1);  [s2,i2] = max(spread(r2));  g2 = gx(r2);
nodeIdx = arrayfun(@(x) find(abs(gx-x)<1e-9),t(3:end));
L5{end+1} = '';
L5{end+1} = 'Spread of the four curves over the interval (max - min):';
L5{end+1} = sprintf('   at the nodes t = 4,6,8,10,12 h : at most %.1e   (all four pass through the readings)',max(spread(nodeIdx)));
L5{end+1} = sprintf('   inside 4-12 h                 : at most %.4f deg C, at t = %.2f h',s2,g2(i2));
L5{end+1} = sprintf('   inside 0-4 h                  : up to  %.4f deg C, at t = %.2f h',s1,g1(i1));
L5{end+1} = sprintf('   at t = 6.4 h                  : %.4f deg C',vmax-vmin);
L5 = [L5, W('',sprintf(['Interpretation: at 6.4 h the four values agree to %.4f deg C because the point is well inside the table and every curve is anchored by readings on both sides of it. ' ...
    'All four pass through the readings at t = 4-12 h, so they agree exactly at those nodes, and between them the spread stays below %.4f deg C (largest at t = %.2f h, in the last interval where there are no readings to the right to steady the curves). ' ...
    'They part company at the left end: the backward polynomial (nodes 4-12 h) and the Bessel polynomial (nodes 2-12 h) have no readings at t = 0 h and are extrapolating there, so inside 0-4 h the spread grows to %.4f deg C at t = %.2f h, while the spline stays close to the data. ' ...
    'Stirling and the backward formula would be the same polynomial at degree 6; here the backward formula is truncated at degree %d, where p = %.1f makes its terms large and it leans on the far end of the table.'],vmax-vmin,s2,g2(i2),s1,g1(i1),DEG_BACK,pB))];

%% REPORT IN THE COMMAND WINDOW
allL = {L1,L2,L3,L4,L5};
for q = 1:numel(allL)
    fprintf('\n%s\n',repmat('=',1,110));
    fprintf('%s\n',allL{q}{:});
end

%% PACK RESULTS AND OPEN THE WINDOW
R.t = t; R.T = T; R.D = D; R.h = hh; R.S = S; R.P = P; R.tabCell = tabCell;
R.L = allL; R.valTab = valTab; R.gx = gx; R.CC = CC; R.spread = spread; R.names = names;
R.v64 = vals; R.notes1 = L1(n+4:end); R.LIMIT = LIMIT; R.tx = tx; R.tmax = tmax; R.Pm = Pm; R.estx = [est.x]; R.estVal = estVal;
useTabs = ~exist('OCTAVE_VERSION','builtin');
if useTabs
    try
        buildWindow(R);
    catch ME
        warning('Window could not be built (%s). Showing plain figure instead.',ME.message);
        f = figure('Name','Cold chain lab','NumberTitle','off');
        drawPlots(subplot(2,1,1),subplot(2,1,2),R);
    end
else
    f = figure('Name','Cold chain lab','NumberTitle','off','Position',[100 100 900 700]); %#ok<NASGU>
    drawPlots(subplot(2,1,1),subplot(2,1,2),R);
end

%% LOCAL FUNCTIONS (must stay at the end of the file)
function D = buildDiff(y)
% D(i,k+1) = k-th forward difference starting at row i (NaN where it does not exist)
n = numel(y);  D = nan(n,n);  D(:,1) = y(:);
for k = 2:n
    D(1:n-k+1,k) = diff(D(1:n-k+2,k-1));
end
end

function [v,lab,tv] = stirF(D,m,p,deg)
% Stirling central formula, origin row m, p = (x - x0)/h, terms up to 'deg'
d = @(k,j) D(m+j,k+1);
lab = {'y0','p*mu(d1)','p^2/2! * d2','p(p^2-1)/3! * mu(d3)','p^2(p^2-1)/4! * d4','p(p^2-1)(p^2-4)/5! * mu(d5)','p^2(p^2-1)(p^2-4)/6! * d6'};
tv = zeros(1,deg+1);
tv(1) = d(0,0);
if deg >= 1, tv(2) = p*(d(1,-1)+d(1,0))/2; end
if deg >= 2, tv(3) = p^2/2*d(2,-1); end
if deg >= 3, tv(4) = p*(p^2-1)/6*(d(3,-2)+d(3,-1))/2; end
if deg >= 4, tv(5) = p^2*(p^2-1)/24*d(4,-2); end
if deg >= 5, tv(6) = p*(p^2-1)*(p^2-4)/120*(d(5,-3)+d(5,-2))/2; end
if deg >= 6, tv(7) = p^2*(p^2-1)*(p^2-4)/720*d(6,-3); end
lab = lab(1:deg+1);  v = sum(tv);
end

function [v,lab,tv] = bessF(D,m,p,deg)
% Bessel central formula for the interval [x0, x0+h], origin row m, terms up to 'deg' (max 5)
d = @(k,j) D(m+j,k+1);
lab = {'(y0+y1)/2','(p-1/2)*d1','p(p-1)/2! * mu(d2)','(p-1/2)p(p-1)/3! * d3','(p+1)p(p-1)(p-2)/4! * mu(d4)','(p-1/2)(p+1)p(p-1)(p-2)/5! * d5'};
tv = zeros(1,deg+1);
tv(1) = (d(0,0)+d(0,1))/2;
if deg >= 1, tv(2) = (p-0.5)*d(1,0); end
if deg >= 2, tv(3) = p*(p-1)/2*(d(2,-1)+d(2,0))/2; end
if deg >= 3, tv(4) = (p-0.5)*p*(p-1)/6*d(3,-1); end
if deg >= 4, tv(5) = (p+1)*p*(p-1)*(p-2)/24*(d(4,-2)+d(4,-1))/2; end
if deg >= 5, tv(6) = (p-0.5)*(p+1)*p*(p-1)*(p-2)/120*d(5,-2); end
lab = lab(1:deg+1);  v = sum(tv);
end

function [v,lab,tv] = backF(D,p,deg)
% Newton backward formula from the last node, p = (x - x_n)/h
n = size(D,1);  tv = zeros(1,deg+1);  lab = cell(1,deg+1);  c = 1;
for k = 0:deg
    tv(k+1) = c*D(n-k,k+1);
    lab{k+1} = sprintf('[p(p+1)..(p+%d)]/%d! * grad^%d y_n',max(k-1,0),k,k);
    c = c*(p+k)/(k+1);
end
lab{1} = 'y_n';
v = sum(tv);
end

function [v,lab,tv] = fwdF(D,p,deg)
% Newton forward formula from the first node, p = (x - x_0)/h
tv = zeros(1,deg+1);  lab = cell(1,deg+1);  c = 1;
for k = 0:deg
    tv(k+1) = c*D(1,k+1);
    lab{k+1} = sprintf('[p(p-1)..(p-%d)]/%d! * d%d y_0',max(k-1,0),k,k);
    c = c*(p-k)/(k+1);
end
lab{1} = 'y_0';
v = sum(tv);
end

function y = lagr(xn,yn,x)
% Lagrange interpolation polynomial through (xn,yn), evaluated at x
y = zeros(size(x));
for i = 1:numel(xn)
    Li = ones(size(x));
    for j = [1:i-1, i+1:numel(xn)]
        Li = Li.*(x - xn(j))/(xn(i) - xn(j));
    end
    y = y + yn(i)*Li;
end
end

function S = nsCoef(x,y)
% Natural cubic spline coefficients (second derivative zero at both ends)
x = x(:);  y = y(:);  n = numel(x);  hs = diff(x);
A = zeros(n);  rhs = zeros(n,1);
A(1,1) = 1;  A(n,n) = 1;
for i = 2:n-1
    A(i,i-1) = hs(i-1);  A(i,i) = 2*(hs(i-1)+hs(i));  A(i,i+1) = hs(i);
    rhs(i) = 3*((y(i+1)-y(i))/hs(i) - (y(i)-y(i-1))/hs(i-1));
end
c = A\rhs;
b = zeros(n-1,1);  d = zeros(n-1,1);
for i = 1:n-1
    b(i) = (y(i+1)-y(i))/hs(i) - hs(i)*(2*c(i)+c(i+1))/3;
    d(i) = (c(i+1)-c(i))/(3*hs(i));
end
S.x = x;  S.a = y(1:n-1);  S.b = b;  S.c = c(1:n-1);  S.d = d;
end

function yq = nsEval(S,xq)
% Evaluate the natural spline S at xq (end pieces are continued outside the table)
yq = zeros(size(xq));
nseg = numel(S.a);
for q = 1:numel(xq)
    i = find(S.x(1:nseg) <= xq(q),1,'last');
    if isempty(i), i = 1; end
    dx = xq(q) - S.x(i);
    yq(q) = S.a(i) + S.b(i)*dx + S.c(i)*dx^2 + S.d(i)*dx^3;
end
end

function ta = timeAbove(F,tx,t0,t1,lim)
% Total length of [t0,t1] on which F > lim, given the sorted crossing times tx
edges = [t0; tx(:); t1];  ta = 0;
for i = 1:numel(edges)-1
    mid = (edges(i)+edges(i+1))/2;
    if F(mid) > lim, ta = ta + edges(i+1) - edges(i); end
end
end

function out = wrapText(pre,s,w)
% Break the sentence s into lines of at most w characters, first line starts with pre
words = strsplit(strtrim(s),' ');
ind = repmat(' ',1,numel(pre));
out = {};  cur = pre;  fresh = true;
for i = 1:numel(words)
    if ~fresh && numel(cur)+1+numel(words{i}) > w
        out{end+1} = cur;  cur = [ind words{i}]; %#ok<AGROW>
    elseif fresh
        cur = [cur words{i}];
    else
        cur = [cur ' ' words{i}];
    end
    fresh = false;
end
out{end+1} = cur;
end

function drawPlots(ax1,ax2,R)
% Top: readings and the four curves.  Bottom: spread of the four curves.
gx = R.gx;  CC = R.CC;
axes(ax1); cla(ax1); hold(ax1,'on');
cols = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56];
sty  = {'-','--','-.',':'};
for q = 1:4
    plot(ax1,gx,CC(q,:),sty{q},'Color',cols(q,:),'LineWidth',1.8);
end
plot(ax1,R.t,R.T,'ko','MarkerFaceColor','k','MarkerSize',7);
plot(ax1,[R.t(1) R.t(end)],[R.LIMIT R.LIMIT],'r-','LineWidth',1.2);
plot(ax1,R.tx,R.LIMIT*ones(size(R.tx)),'rs','MarkerFaceColor','r','MarkerSize',8);
plot(ax1,R.tmax,R.Pm,'m^','MarkerFaceColor','m','MarkerSize',9);
plot(ax1,6.4*ones(1,4),R.v64,'kd','MarkerSize',8);
ylim(ax1,[floor(min(R.T))-1, ceil(max(R.T))+2]);
xlim(ax1,[0 12]);  grid(ax1,'on');  box(ax1,'on');
xlabel(ax1,'t (h)');  ylabel(ax1,'T (deg C)');
title(ax1,'Product temperature: readings and the four interpolants');
legend(ax1,[R.names {'Readings','8.0 deg C limit','Crossings','Maximum','t = 6.4 h'}],'Location','northwest');
hold(ax1,'off');
axes(ax2); cla(ax2);
plot(ax2,gx,R.spread,'k-','LineWidth',1.5);
hold(ax2,'on'); plot(ax2,R.t,zeros(size(R.t)),'ko','MarkerFaceColor','k'); hold(ax2,'off');
xlim(ax2,[0 12]);  grid(ax2,'on');
xlabel(ax2,'t (h)');  ylabel(ax2,'spread (deg C)');
title(ax2,'Highest minus lowest of the four curves');
end

function buildWindow(R)
% Tabbed window with the report, the tables, a quick estimator and the plots
f = figure('Name','Cold chain lab','NumberTitle','off','Position',[80 60 1180 740]);
tg = uitabgroup(f,'Position',[0 0 1 1]);
nm = {'Difference table','Estimates','Maximum','Crossings','Validation','Plots'};
tb = cell(1,6);
for q = 1:6, tb{q} = uitab(tg,'Title',nm{q}); end
fw = get(0,'FixedWidthFontName');
mk = @(par,pos,lines) uicontrol(par,'Style','listbox','Units','normalized','Position',pos,'String',lines,'FontName',fw,'FontSize',10,'BackgroundColor','w','HorizontalAlignment','left');

% tab 1: difference table
cn = {'t (h)','T','Delta 1','Delta 2','Delta 3','Delta 4','Delta 5','Delta 6'};
uitable(tb{1},'Units','normalized','Position',[0.01 0.52 0.98 0.46],'Data',R.tabCell,'ColumnName',cn,'ColumnWidth',{70,90,90,90,90,90,90,90},'FontSize',11);
mk(tb{1},[0.01 0.01 0.98 0.49],R.notes1);

% tab 2: estimates plus quick estimator
mk(tb{2},[0.01 0.14 0.98 0.85],R.L{2});
ui.edit = uicontrol(tb{2},'Style','edit','Units','normalized','Position',[0.12 0.04 0.08 0.05],'String','7.0','FontSize',11);
uicontrol(tb{2},'Style','text','Units','normalized','Position',[0.01 0.04 0.10 0.05],'String','Quick estimate t (h):','HorizontalAlignment','right','FontSize',10);
ui.pop = uicontrol(tb{2},'Style','popupmenu','Units','normalized','Position',[0.22 0.04 0.20 0.05],'FontSize',10, ...
    'String',{'Newton forward (deg 6)','Newton backward (deg 6)','Stirling (deg 6)','Bessel (deg 5)','Lagrange (deg 6)','Natural cubic spline'});
ui.out = uicontrol(tb{2},'Style','text','Units','normalized','Position',[0.01 0.005 0.98 0.035],'String','','HorizontalAlignment','left','FontSize',11,'FontName',fw);
uicontrol(tb{2},'Style','pushbutton','Units','normalized','Position',[0.44 0.04 0.10 0.05],'String','Estimate','FontSize',11, ...
    'Callback',@(s,e) quickEst(ui,R));

% tabs 3 and 4: text
mk(tb{3},[0.01 0.01 0.98 0.98],R.L{3});
mk(tb{4},[0.01 0.01 0.98 0.98],R.L{4});

% tab 5: validation table + text
uitable(tb{5},'Units','normalized','Position',[0.01 0.70 0.98 0.28],'Data',R.valTab,'ColumnName',{'Method','Settings','T(6.4 h) deg C','Minus lowest'}, ...
    'ColumnWidth',{150,420,120,110},'FontSize',11);
mk(tb{5},[0.01 0.01 0.98 0.67],R.L{5});

% tab 6: plots
ax1 = axes('Parent',tb{6},'Position',[0.07 0.36 0.90 0.58]);
ax2 = axes('Parent',tb{6},'Position',[0.07 0.07 0.90 0.20]);
drawPlots(ax1,ax2,R);
end

function quickEst(ui,R)
% Callback of the Estimate button: value at any t with the chosen formula
x = str2double(get(ui.edit,'String'));
if isnan(x), set(ui.out,'String','Enter a number for t.'); return; end
t = R.t;  D = R.D;  h = R.h;  n = numel(t);
sel = get(ui.pop,'Value');  note = '';
switch sel
    case 1, p = (x-t(1))/h;   v = fwdF(D,p,6);
    case 2, p = (x-t(n))/h;   v = backF(D,p,6);
    case 3, p = (x-t(4))/h;   v = stirF(D,4,p,6);
    case 4
        if x < t(4), m = 3; else, m = 4; end
        p = (x-t(m))/h;       v = bessF(D,m,p,5);
    case 5, p = NaN;          v = lagr(t,R.T,x);
    otherwise, p = NaN;       v = nsEval(R.S,x);
end
if x < t(1) || x > t(n), note = '  (EXTRAPOLATION - not reliable)'; end
if isnan(p), ps = ''; else, ps = sprintf('   p = %.4f',p); end
if sel == 4 && (p < 0 || p > 1), note = [note '  (outside the Bessel interval)']; end
set(ui.out,'String',sprintf('T(%.4f h) = %.4f deg C%s%s',x,v,ps,note));
end
