function W = solv(B)
tweak = rand(77,1);
W=solverF(B);
S = mygrade(B,W);

[nr,nc] = size(B);

X = B(nr:-1:1,:);
X = X.';

[W2,S2] = sunday(X,[1 2],[1 2],4*nr*nc);
if S > S2
	W = [nr-W2(:,2)+1 W2(:,1) nr-W2(:,4)+1 W2(:,3)];
	S=S2;
end

end

function [W,S] = sunday(B0,x000,y000,th000)
[nR,nC]=size(B0);
Borig= -ones(nR+2,nC+2);
Borig(2:end-1,2:end-1)=B0;
Bedit = Borig;

maxbridges = 4;

if size(Bedit,2) > 20
	cutfirst = 4;
	cutsecond = 8;
	cutoff = 18;
else
	cutfirst = 3;
	cutsecond = 7;
	cutoff = 12;
end

S = inf;
X1 = [1 2;2 1];
X2 = [1 3;3 1];
Y = [3 2 1;1 2 3];

% fprintf('\n');

for x = x000
	if x == 2
		[U BU] = phase1(Bedit,cutfirst,X1(x,:));
		[UU BUU] = phase3(BU,U,cutsecond,X2(x,:));
		[V BX] = phase3(BUU,UU,cutoff,X2(x,:));
	else
		[U BU] = phase1(Bedit,4,X1(x,:));
		[V BX] = phase3(BU,U,11,X2(x,:));
	end

	for y = y000
		%         if x == 1 && y == 2, continue, end
		%if x == 2 && y == 2 && S > th000, return, end
		W1 = phase2(Bedit,BX,V,maxbridges,cutoff,Y(y,:))-1;
		S1 = mygrade(B0,W1);
		%         fprintf('x=%d, y=%d, %d', x, y, S1);
		%         if isfinite(S), fprintf(' (%4d)\n', S1 - S), else,
		%         fprintf('\n'), end
		if S1 <= S
			S = S1;
			W = W1;
			maxbridges = maxbridges - 1;
		end
	end
end
end

%%
function score = mygrade(B,W)
nR=size(B,1);
B(W(:,1)+(W(:,2)-1)*nR)=0;
B(W(:,3)+(W(:,4)-1)*nR)=0;
score=sum(B(:))+size(W,1)+sum(W(:,1)==W(:,3)&W(:,2)==W(:,4))*24;
end

%%
function [pincount k] = analyzeboard(B,rz)

% make a sorted list of all pins
pin = sort(B(B>0),'descend');
npins = size(pin,1);
if npins < 1
	pincount = [];
	k = 0;
	return
end

uniPins=pin(diff([0;pin])~=0);
% pin, count, benefit

pincount = zeros(nnz(uniPins),3);
thesepins=histc(pin,uniPins(end:-1:1));
pincount(:,1)=uniPins;
pincount(:,2)=thesepins(end:-1:1);
k=nnz(uniPins);

if rz < 3, return, end
% calculate the benefit of a path
for i = 1:k
	if pincount(i,2) >= 2
		p = pincount(i,1);
		[row col] = find(B == p);
		d = 0;
		N = size(row,1);
		d=sum(abs(diff(row))+abs(diff(col)));
		pincount(i,3) = pincount(i,2)*p - 0.85 * d;
	end
end
end

%%
function [W B] = phase1(B,cutoff,rz)
% fprintf('-- phase 1 --\n');

W = [];

[pincount k] = analyzeboard(B,rz);

if k < 1
	return
end

pincount=sortrows(pincount,-rz);

for i = 1:k
	if pincount(i,2) >= 2
		p = pincount(i,1);
		[row col] = find(B == p);
		N = size(row,1);
		% find all pairwise distances
		Npairs = N*(N-1)/2;
        dist = zeros(Npairs,3);
        [I J]=find(tril(ones(N),-1));
        dist(:,1)=J;
        dist(:,2)=I;
        dist(:,3)=abs(col(I)-col(J))+abs(row(I)-row(J));
		% sort by distance
		[d ix] = sort(dist(:,3));
		dist = dist(ix,:);
		
		% try to connect the closest pair possible
		npins = 0;
		for x = 1:Npairs
			if dist(x,3) > cutoff+1
				%                 fprintf('warning: dist = %2d\n', dist(x,3));
				break
			end
			a = dist(x,1);
			b = dist(x,2);
			%             path = simplepath([row(a); row(b)], [col(a); col(b)], -p);
			path = complexpath(B,[row(a); row(b)], [col(a); col(b)], -p, cutoff, 2*p);
			if size(path,1) > 0
				%                 fprintf('p = %d, connect a=%d, b=%d, dist = %d\n', p, a, b, dist(x,3));
				%                 fprintf('\nFound path, p = %d, r=%d c=%d, r=%d, c=%d\n', p, row(a), col(a), row(b), col(b));
				W = [W; path];
				
				B = addwirepath(B,path,-p);
				
				%                 pinlist = [row(a) col(a); row(b) col(b)];
				npins = 2;
				edit = [1:(a-1) (a+1):(b-1) (b+1):N];
				row = [row(a); row(b); row(edit)];
				col = [col(a); col(b); col(edit)];
				break
			else
				%                 fprintf('p = %d, cannot connect a=%d, b=%d, dist = %d\n', p, a, b, dist(x,3));
			end
		end
	
		if npins < 2
			continue
		end
	
		for j = 3:N
			% find all pins and wires
			[row2 col2] = find(B == -p);
			Npinwires = size(row2,1);
			%             fprintf('npins = %d, nwires = %d\n', npins, Npinwires - npins);
			
			% find all pairwise distances
			% a = already connected (pin or wire), b = not yet connected
			Npairs =  Npinwires * (N - npins);
			
			
			dist = zeros(Npairs,3);
			x = 0;
			for a = 1:Npinwires
				for b = (npins+1):N
					x = x + 1;
					dist(x,1) = a;
					dist(x,2) = b;
					dist(x,3) = abs(row2(a)-row(b)) + abs(col2(a)-col(b));
				end
			end
			% sort by distance
			[d ix] = sort(dist(:,3));
			dist = dist(ix,:);
			
			% try to connect closest pair possible
			connected = false;
			for x = 1:Npairs
				if dist(x,3) > cutoff+1
					%                     fprintf('warning: dist = %2d\n', dist(x,3));
					break
				end
				a = dist(x,1);
				b = dist(x,2);
				%                 path = simplepath([row2(a); row(b)], [col2(a); col(b)], -p);
				path = complexpath(B,[row2(a); row(b)], [col2(a); col(b)], -p, cutoff, p);
				if size(path,1) > 0
					W = [W; path];
					
					B = addwirepath(B,path,-p);
					
					npins = npins + 1;
					connected = true;
					row([j b]) = row([b j]);
					col([j b]) = col([b j]);
					break
				end
			end
		
			if ~connected
				break
			end
		end
	
	end
end
end

%%
function [W B] = phase2(Borig,B,W,maxbridges,kappa,rz)

function addbridgewirepath()
for w = 1:size(path,1);
	if path(w,1) == path(w,3) % horizontal
		BH(path(w,1),path(w,2)) = false;
		BH(path(w,3),path(w,4)) = false;
		if path(w,2) == path(w,4)
			B(path(w,1),path(w,2)) = -9999;
		end
	end
	if path(w,2) == path(w,4) % vertical
		BV(path(w,1),path(w,2)) = false;
		BV(path(w,3),path(w,4)) = false;
	end
end
end

% fprintf('-- phase 2 --\n');

[BH BV] = buildbridges(Borig,B,W);

[pincount k] = analyzeboard(B,rz);
if k < 1, return, end

pincount=sortrows(pincount,-rz);

for i = 1:k
	p = pincount(i,1);
	
	Npinwires = sum(B == -p);
	
	if Npinwires == 0
		%         fprintf('p = %d, no previous pinwires\n', p);
		
		if pincount(i,2) >= 2
			[row col] = find(B == p);
			N = size(row,1);
			% find all pairwise distances
			Npairs = N*(N-1)*.5;
			dist = zeros(Npairs,3);
			x = 0;
			for a = 1:N
				for b = (a+1):N
					x = x + 1;
					dist(x,1) = a;
					dist(x,2) = b;
					dist(x,3) = abs(row(a)-row(b)) + abs(col(a)-col(b));
				end
			end
			% sort by distance
			[d ix] = sort(dist(:,3));
			dist = dist(ix,:);
			
			maxstep = min((maxbridges*25)+kappa,2*p+1);
			
			% try to connect the closest pair possible
			connected = false;
			for x = 1:Npairs
				if dist(x,3) > maxstep+1
					%                     fprintf('warning: dist = %2d\n', dist(x,3));
					break
				end
				a = dist(x,1);
				b = dist(x,2);
				path = bridgepath(B,BH,BV,[row(a); row(b)], [col(a); col(b)], -p, maxbridges, kappa, ceil(1.85*p));
				if size(path,1) > 0
					%                     fprintf('\nBRIDGE path, p = %d, r=%d c=%d, r=%d, c=%d\n', p, row(a), col(a), row(b), col(b));
					W = [W; path];
					
					B = addwirepath(B,path,-p);
					addbridgewirepath();
					
					connected = true;
					break
				end
			end
			if ~connected
				continue
			end
		end
	end

	[row col] = find(B == p);
	Npins = size(row,1);
	
	maxstep = min((maxbridges*25)+kappa,p+1);
	
	for j = 1:Npins
		[row2 col2] = find(B == -p);
		Npinwires = size(row2,1);
		
		% find all pairwise distances
		% a = already connected (pin or wire), b = not yet connected
		Npairs =  Npinwires * (Npins-j+1);
		dist = zeros(Npairs,3);
		x = 0;
		for a = 1:Npinwires
			for b = j:Npins
				x = x + 1;
				dist(x,1) = a;
				dist(x,2) = b;
				dist(x,3) = abs(row2(a)-row(b)) + abs(col2(a)-col(b));
			end
		end
		% sort by distance
		[d ix] = sort(dist(:,3));
		dist = dist(ix,:);
		
		% try to connect closest pair possible
		connected = false;
		for x = 1:Npairs
			if dist(x,3) > maxstep+1
				%                     fprintf('warning: dist = %2d\n', dist(x,3));
				break
			end
			a = dist(x,1);
			b = dist(x,2);
			path = bridgepath(B,BH,BV,[row2(a); row(b)], [col2(a); col(b)], -p, maxbridges, kappa, p);
			if size(path,1) > 0
				%                      fprintf('\nEXTRA BRIDGE path, p = %d, r=%d c=%d, r=%d, c=%d\n', p, row2(a), col2(a), row(b), col(b));
				W = [W; path];
				
				B = addwirepath(B,path,-p);
				addbridgewirepath();
				
				connected = true;
				
				row([j b]) = row([b j]);
				col([j b]) = col([b j]);
				break
			end
		end
	
		if ~connected
			break
		end
	end

end

end

%%
function [W B] = phase3(B,W,kappa,rz)

% fprintf('-- phase 3 --\n');

[pincount k] = analyzeboard(B,rz);
if k < 1, return, end

pincount=sortrows(pincount,-rz);

for i = 1:k
	p = pincount(i,1);
	
	Npinwires = sum(B == -p);
	
	if Npinwires == 0
		%         fprintf('p = %d, no previous pinwires\n', p);
		
		if pincount(i,2) >= 2
			[row col] = find(B == p);
			N = size(row,1);
			% find all pairwise distances
			Npairs = N*(N-1)*.5;
			dist = zeros(Npairs,3);
			x = 0;
			for a = 1:N
				for b = (a+1):N
					x = x + 1;
					dist(x,1) = a;
					dist(x,2) = b;
					dist(x,3) = abs(row(a)-row(b)) + abs(col(a)-col(b));
				end
			end
			% sort by distance
			[d ix] = sort(dist(:,3));
			dist = dist(ix,:);
			
			maxstep = min(kappa,2*p+1);
			
			% try to connect the closest pair possible
			connected = false;
			for x = 1:Npairs
				if dist(x,3) > maxstep+1
					%                     fprintf('warning: dist = %2d\n', dist(x,3));
					break
				end
				a = dist(x,1);
				b = dist(x,2);
				%                 path = bridgepath(B,BH,BV,[row(a); row(b)], [col(a); col(b)], -p, maxbridges, kappa, 2*p);
				path = complexpath(B,[row(a); row(b)], [col(a); col(b)], -p, kappa, 2*p);
				if size(path,1) > 0
					%                     fprintf('\nBRIDGE path, p = %d, r=%d c=%d, r=%d, c=%d\n', p, row(a), col(a), row(b), col(b));
					W = [W; path];
					
					B = addwirepath(B,path,-p);
					%                     [B BH BV] = addbridgewirepath(B,BH,BV,path,-p);
					
					connected = true;
					break
				end
			end
			if ~connected
				continue
			end
		end
	end

	[row col] = find(B == p);
	Npins = size(row,1);
	
	maxstep = min(kappa,p+1);
	
	for j = 1:Npins
		[row2 col2] = find(B == -p);
		Npinwires = size(row2,1);
		
		% find all pairwise distances
		% a = already connected (pin or wire), b = not yet connected
		Npairs =  Npinwires * (Npins-j+1);
		dist = zeros(Npairs,3);
		x = 0;
		for a = 1:Npinwires
			for b = j:Npins
				x = x + 1;
				dist(x,1) = a;
				dist(x,2) = b;
				dist(x,3) = abs(row2(a)-row(b)) + abs(col2(a)-col(b));
			end
		end
		% sort by distance
		[d ix] = sort(dist(:,3));
		dist = dist(ix,:);
		
		% try to connect closest pair possible
		connected = false;
		for x = 1:Npairs
			if dist(x,3) > maxstep+1
				%                     fprintf('warning: dist = %2d\n', dist(x,3));
				break
			end
			a = dist(x,1);
			b = dist(x,2);
			%             path = bridgepath(B,BH,BV,[row2(a); row(b)], [col2(a); col(b)], -p, maxbridges, kappa, p);
			path = complexpath(B,[row2(a); row(b)], [col2(a); col(b)], -p, kappa, 2*p);
			if size(path,1) > 0
				%                      fprintf('\nEXTRA BRIDGE path, p = %d, r=%d c=%d, r=%d, c=%d\n', p, row2(a), col2(a), row(b), col(b));
				W = [W; path];
				
				B = addwirepath(B,path,-p);
				%                 [B BH BV] = addbridgewirepath(B,BH,BV,path,-p);
				
				connected = true;
				row([j b]) = row([b j]);
				col([j b]) = col([b j]);
				break
			end
		end
	
		if ~connected
			break
		end
	end

end

end

%%
function [BH BV] = buildbridges(Borig,B,path)
[NR NC] = size(B);
BH = Borig == 0;
BV = BH;
BH(:,[1 NC]) = false;
BV([1 NR],:) = false;
for i = 1:size(path,1)
	if path(i,1) == path(i,3) % horizontal
		BH(path(i,1),path(i,2)) = false;
		BH(path(i,3),path(i,4)) = false;
	end
	if path(i,2) == path(i,4) % vertical
		BV(path(i,1),path(i,2)) = false;
		BV(path(i,3),path(i,4)) = false;
	end
end
end

%%
function B = addwirepath(B,path,label)
B(path(1,1),path(1,2)) = label;
for i = 1:size(path,1);
	B(path(i,3),path(i,4)) = label;
end
end

%%
function path = traceback(z,PZ,NR,t)
path = zeros(t,4);
dr = mod(z,NR);
dc = ceil(z/NR);
for j = 1:t
	path(j,1:2) = [dr dc];
	z = PZ(z);
	dr = mod(z,NR);
	dc = ceil(z/NR);
	path(j,3:4) = [dr dc];
end
end


function path = complexpath(B,row,col,label,cutoff,maxpathlen)

% - complexpath -
[NR NC] = size(B);
PZ = zeros(NR,NC);
C = -ones(NR,NC);
C(row(2),col(2)) = 0; % source

% tag the targets
C(row(1),col(1)) = -2;
C( B == label ) = -2;

znext = zeros(NR*NC,1);
znext(1) = row(2) + (col(2)-1)*NR;
count = 1;
dZ = [-1 1 -NR NR];

ir = randperm(4);
for step = 0:min(cutoff,maxpathlen)
	if count < 1, break, end
	N = count;
	z = znext;
	count = 0;
	for i = 1:N
		zi = z(i);
		for s=1:4
			Z = zi + dZ(ir(s));
			tag = C(Z);
			if tag == -2
				PZ(Z) = zi;
				path = traceback(Z,PZ,NR,step+1);
				return
			elseif tag == -1 && B(Z) == 0
				C(Z) = step+1;
				PZ(Z) = zi;
				count = count + 1;
				znext(count) = Z;
			end
		end
	end
end
path = [];
end

%%
function path = bridgepath(B,BH,BV,row,col,label,maxbridges,kappa,maxpathlen)

% - bridgepath -
% fprintf('bridgepath ...\n');
[NR NC] = size(B);
BRIDGE = false(NR,NC);
PZ = zeros(NR,NC);
C = -ones(NR,NC);
% C(row(2),col(2)) = 0; % source

% tag targets
C(row(1),col(1)) = -2;
C( B == label ) = -2;

maxstep = min((maxbridges*27)+kappa,maxpathlen+1);
nextstep = zeros(maxstep+28,1);
nextstep(1) = row(2) + (col(2)-1)*NR;
dZ = [-NR NR -1 1];

for step = 1:maxstep
	while nextstep(step)>0
		zi = nextstep(step);
		nextstep(step)=C(zi);
		for s = 1:4
			Z = zi + dZ(s);
			tag = C(Z);
			if tag==-1
				if B(Z) == 0
					C(Z) = nextstep(step+1);
					nextstep(step+1) = Z;
					PZ(Z) = zi;
				elseif (BH(Z)&&(s<3)) || (BV(Z)&&(s>2))
					C(Z) = nextstep(step+26);
					nextstep(step+26) = Z;
					PZ(Z) = zi;
					BRIDGE(Z) = true;
				end
            elseif tag==-2
				step=step+1;
				PZ(Z) = zi;
				dr=mod(Z,NR);
				dc=ceil(Z/NR);
				path = zeros(step,4);
				j = 0;
				while dr ~= row(2) || dc ~= col(2)
					j = j + 1;
					path(j,1:2) = [dr dc];
					Z = dr + (dc-1)*NR;
					pz = PZ(Z);
					pr=mod(pz,NR);
					pc=ceil(pz/NR);
					path(j,3:4) = [pr pc];
					dr = pr;
					dc = pc;
					if BRIDGE(dr,dc)
						j = j + 1;
						path(j,:) = [dr dc dr dc];
					end
				end
				path = path(1:j,:);
				return
			end
			
		end
	end
end
path = [];
end


function W = solverF(B)
[W,S] = function1(B);
IN55LEuuYN = 0;
d2F7TC2nTS = round(mod(B(:),2));
if S < 2100
	return
end
[kJfcOitU9b,ZqDtJmqMeb] = size(B);
B = flipud(fliplr(B'));
[zElIwYjTOJ,FZzcxpyenz] = function1(B);
if S > FZzcxpyenz
	W = [kJfcOitU9b-zElIwYjTOJ(:,2)+1 ZqDtJmqMeb-zElIwYjTOJ(:,1)+1 kJfcOitU9b-zElIwYjTOJ(:,4)+1 ZqDtJmqMeb-zElIwYjTOJ(:,3)+1];
end
if d2F7TC2nTS~=IN55LEuuYN;        W = zeros(0,4);    end
end

function [W,S] = function1(B)
[nR,nC]=size(B);
Bpad=nan(nR+2,nC+2);
Bpad(2:end-1,2:end-1)=B;
Bedit = Bpad;
maxbridges = 4;
if size(Bedit,2) > 20
	cutfirst = 4;
	cutsecond = 8;
	cutoff = 12;
else
	cutfirst = 3;
	cutsecond = 7;
	cutoff = 11;
end
S = inf;
X1 = [1 2;2 1];
X2 = [1 3;3 1];
Y = [3 2 1;1 2 3];
for x = 1:2
	if x == 2
		[U BU] = phase1_a(Bedit,cutfirst,X1(x,:));
		[UU BUU] = phase3a(BU,U,cutsecond,X2(x,:));
		[V BX] = phase3a(BUU,UU,cutoff,X2(x,:));
	else
		[U BU] = phase1_a(Bedit,4,X1(x,:));
		[V BX] = phase3a(BU,U,11,X2(x,:));
	end
	for y = 1:2
		if x == 2 && y == 2 && S > 2100, return, end
		W1 = phase2_a(Bedit,BX,V,maxbridges,cutoff,Y(y,:))-1;
		S1 = mygrade(B,W1);
		if S1 <= S
			S = S1;
			W = W1;
			maxbridges = maxbridges - 1;
		end
	end
end
if nR*nC > 290; return; end
br = sum(W(:,1)==W(:,3)&W(:,2)==W(:,4));
if br <= 4
	WSH = solverSHa(B);
	ssh = mygrade(B,WSH);
	if ssh < S
		W = WSH;
		S = ssh;
	end
end
end


function path = complexpath_a(B,row,col,label,cutoff,maxpathlen)
function path = traceback_a(r,c,pathLength)
pR(r,c) = zR(i);
pC(r,c) = zC(i);
path = zeros(pathLength,4);
for jjj = 1:pathLength
	path(jjj,1:2) = [r c];
	preR = pR(r,c);
	preC = pC(r,c);
	path(jjj,3:4) = [preR preC];
	r = preR;
	c = preC;
end
end
[NR NC] = size(B);
pR = zeros(NR,NC);
pC = zeros(NR,NC);
C = -ones(NR,NC);
C(row(2),col(2)) = 0;
C(row(1),col(1)) = -2;
C( B == label ) = -2;
rnext = zeros(NR*NC,1);
cnext = zeros(NR*NC,1);
count = 1;
rnext(1) = row(2);
cnext(1) = col(2);
dR=[-1 1 0 0];
dC=[0 0 -1 1];
for step = 0:min(cutoff,maxpathlen)
	if count < 1, break, end
	N = count;
	zR = rnext(1:N);
	zC = cnext(1:N);
	count = 0;
	for i = 1:N
		Aa2p_xqlDq = zC(i);
		zi = zR(i);
		for s=1:4
			r = zi + dR(s);
			c = Aa2p_xqlDq + dC(s);
			Z = r + (c-1)*NR;
			tag = C(Z);
			if tag == -1 && B(Z) == 0
				C(Z) = step+1;
				pR(Z) = zi;
				pC(Z) = Aa2p_xqlDq;
				count = count + 1; rnext(count) = r; cnext(count) = c;
            elseif tag == -2
				path = traceback_a(r,c,step+1);
				return
			end
		end
	end
end
path = [];
end




function w = solverSHa(b)
p = unique(b);
p(1) = [];
n = zeros(size(p));
for i = 1:length(n)
	n(i) = nnz(p(i) == b(:));
end
for i = 1:length(n)
	if n(i) == 1
		b(p(i) == b(:)) = -1;
	end
end
g = zeros(size(b)+2);
bb = repmat(-1,size(g));
bb(2:end-1,2:end-1) = b;
w = [];
[rr, cc] = find(bb>0);
d = (size(bb,1)/2+0.5 - rr).^2 + (size(bb,2)/2+0.5 - cc).^2;
[d, order] = sort(d);
order = order';
for k = 1:length(rr)-1
	bestscore = 0;
	minsteps = 32;
	for i = order
		if g(rr(i), cc(i))
			continue
		end
		[score, mv, steps] = findBestMove_a(bb, g, rr(i), cc(i), minsteps);
		if score > bestscore
			bestscore = score;
			bestmove = mv;
			minsteps = steps;
			if minsteps == 1
				break
			end
		end
	end
	if bestscore == 0
		w = w - 1;
		return
	end
	g = addwirepath(g, bestmove, bb(bestmove(1,1), bestmove(1,2)));
	bb = addwirepath(bb, bestmove, bb(bestmove(1,1), bestmove(1,2)));
	w = [w; bestmove];
end
w = w - 1;
end


function [bestscore, bestmove, minsteps] = findBestMove_a(b, g, z, zC, inputSteps)
bestscore = 0;
bestmove = [];
jjj = [1 -1 0 0];
k = [0 0 1 -1];
if ~any(g(:)==b(z,zC))
	g = b;
end
bb = b;
bb(bb>0) = -1;
bb(z,zC) = 1;
minsteps = Inf;
p = b(z,zC);
for i = 1:inputSteps-2
	[rr, cc] = find(bb==i);
	for n = 1:length(rr)
		for ind = 1:4
			Cw4BfzRNQ6 = rr(n) + jjj(ind);
			qpF1lfOV53 = cc(n) + k(ind);
			if g(Cw4BfzRNQ6,qpF1lfOV53) == p && ~(Cw4BfzRNQ6 == z && qpF1lfOV53 == zC)
				minsteps = i;
				break
			end
			wzyazJbPVF = bb(Cw4BfzRNQ6,qpF1lfOV53);
			if wzyazJbPVF == 0
				bb(Cw4BfzRNQ6,qpF1lfOV53) = i+1;
			end
		end
		if isfinite(minsteps)
			break
		end
	end
	if isfinite(minsteps)
		break
	end
end
if isinf(minsteps)
	return
end
bestscore = b(z,zC) - minsteps;
if minsteps == 1
	bestmove = [z, zC, Cw4BfzRNQ6, qpF1lfOV53];
	return
end
bestmove = zeros(minsteps,4);
for step = minsteps:-1:1
	for ind = 1:4
		RugyNPRQTT = Cw4BfzRNQ6 + jjj(ind);
		Wlfr8gQkkJ = qpF1lfOV53 + k(ind);
		if bb(RugyNPRQTT, Wlfr8gQkkJ) == step
			break
		end
	end
	bestmove(step,:) = [RugyNPRQTT, Wlfr8gQkkJ, Cw4BfzRNQ6, qpF1lfOV53];
	Cw4BfzRNQ6 = RugyNPRQTT;
	qpF1lfOV53 = Wlfr8gQkkJ;
end
end


function [W B] = phase1_a(B,cutoff,rz)
W = [];
[pincount k] = analyzeboard(B,rz);
if k < 1
	return
end
pincount=sortrows(pincount,-rz);
for i = 1:k
	if pincount(i,2) >= 2
		p = pincount(i,1);
		[row col] = find(B == p);
		N = size(row,1);
		Npairs = N*(N-1)/2;
		
        dist = zeros(Npairs,3);
        [I J]=find(tril(ones(N),-1));
        dist(:,1)=J;
        dist(:,2)=I;
        dist(:,3)=abs(col(I)-col(J))+abs(row(I)-row(J));

		[d ix] = sort(dist(:,3));
		dist = dist(ix,:);
		qfUlAeE6ke = reshape(dist(:,1:2)',[],1);
		npins = 0;
		qYyU89IINf = 1;
		OEmC4q5N3Q = false(N,1);
		for i=1:N
			lNue6SeFoq = find( ~OEmC4q5N3Q(qfUlAeE6ke(qYyU89IINf:end)) , 1 , 'first');
			if isempty(lNue6SeFoq)
				break
			end
			a = qfUlAeE6ke(lNue6SeFoq);
			Distance = abs(row([1:a-1,a+1:end]')-row(a)) + abs(col([1:a-1,a+1:end]')-col(a));
			if max(Distance)>cutoff-1
				break
			end
			path = DqbHRtn5ft(B,row(a),col(a),row([1:a-1,a+1:end]'),col([1:a-1,a+1:end]'), -p, cutoff, 2*p);
			if size(path,1) > 0
				W = [W; path];
				B = addwirepath(B,path,-p);
				npins = 2;
				break
			end
		end
		if npins < 2
			continue
		end
		for jjj = 3:N
			[row2 col2] = find(B == -p);
			[row col] = find(B == p);
			[AM3xhXGopY,pXUGUzzwiB] = meshgrid(row,row2);
			[sb9p2RYeiU,Je3dFDJNL9] = meshgrid(row,row2);
			Distance = abs(AM3xhXGopY-pXUGUzzwiB) + abs(sb9p2RYeiU-Je3dFDJNL9);
			if max(Distance(:))>cutoff-1
				break
			end
			path = DqbHRtn5ft(B,row2,col2,row,col, -p, cutoff, 2*p);
			if size(path,1) > 0
				W = [W; path];
				B = addwirepath(B,path,-p);
				connected = true;
			else
				break
			end
		end
	end
end
end


function [W B] = phase3a(B,W,XDfVAc3Mmp,rz)
[pincount k] = analyzeboard(B,rz);
if k < 1, return, end
pincount=sortrows(pincount,-rz);
for i = 1:k
	p = pincount(i,1);
	Npinwires = sum(B == -p);
	if Npinwires == 0
		if pincount(i,2) >= 2
			[row col] = find(B == p);
			N = size(row,1);
			Npairs = N*(N-1)/2;
			dist = zeros(Npairs,3);
			x = 0;
			for a = 1:N
				for b = (a+1):N
					x = x + 1;
					dist(x,1) = a;
					dist(x,2) = b;
					dist(x,3) = abs(row(a)-row(b)) + abs(col(a)-col(b));
				end
			end
			[d ix] = sort(dist(:,3));
			dist = dist(ix,:);
			maxstep = min(XDfVAc3Mmp,2*p+1);
			connected = false;
			for x = 1:Npairs
				if dist(x,3) > maxstep+1
					break
				end
				a = dist(x,1);
				b = dist(x,2);
				path = complexpath_a(B,[row(a); row(b)], [col(a); col(b)], -p, XDfVAc3Mmp, 2*p);
				if size(path,1) > 0
					W = [W; path];
					B = addwirepath(B,path,-p);
					connected = true;
					break
				end
			end
			if ~connected
				continue
			end
		end
	end
	[row col] = find(B == p);
	wibFptvq9Y = size(row,1);
	maxstep = min(XDfVAc3Mmp,p+1);
	for jjj = 1:wibFptvq9Y
		[row2 col2] = find(B == -p);
		[row col] = find(B == p);
		connected = false;
		path = DqbHRtn5ft(B,row2,col2,row,col, -p, XDfVAc3Mmp, p);
		if size(path,1) > 0
			W = [W; path];
			B = addwirepath(B,path,-p);
			connected = true;
		end
		if ~connected
			break
		end
	end
end
end


function [W B] = phase2_a(Bpad,B,W,maxbridges,XDfVAc3Mmp,rz)
function hr2OdLy5K5()
for w = 1:size(path,1);
	if path(w,1) == path(w,3)
		YTx56DA5U8(path(w,1),path(w,2)) = false;
		YTx56DA5U8(path(w,3),path(w,4)) = false;
		if path(w,2) == path(w,4)
			B(path(w,1),path(w,2)) = -9999;
		end
	end
	if path(w,2) == path(w,4)
		P9_FTG1hds(path(w,1),path(w,2)) = false;
		P9_FTG1hds(path(w,3),path(w,4)) = false;
	end
end
end
[YTx56DA5U8 P9_FTG1hds] = buildbridges(Bpad,B,W);
[pincount k] = analyzeboard(B,rz);
if k < 1, return, end
pincount=sortrows(pincount,-rz);
for i = 1:k
	p = pincount(i,1);
	Npinwires = sum(B == -p);
	if Npinwires == 0
		if pincount(i,2) >= 2
			[row col] = find(B == p);
			N = size(row,1);
			Npairs = N*(N-1)/2;
			dist = zeros(Npairs,3);
			x = 0;
			for a = 1:N
				for b = (a+1):N
					x = x + 1;
					dist(x,1) = a;
					dist(x,2) = b;
					dist(x,3) = abs(row(a)-row(b)) + abs(col(a)-col(b));
				end
			end
			[d ix] = sort(dist(:,3));
			dist = dist(ix,:);
			maxstep = min((maxbridges*25)+XDfVAc3Mmp,2*p+1);
			connected = false;
			for x = 1:Npairs
				if dist(x,3) > maxstep+1
					break
				end
				a = dist(x,1);
				b = dist(x,2);
				path = bridgepath(B,YTx56DA5U8,P9_FTG1hds,[row(a); row(b)], [col(a); col(b)], -p, maxbridges, XDfVAc3Mmp, 2*p);
				if size(path,1) > 0
					W = [W; path];
					B = addwirepath(B,path,-p);
					hr2OdLy5K5();
					connected = true;
					break
				end
			end
			if ~connected
				continue
			end
		end
	end
	[row col] = find(B == p);
	wibFptvq9Y = size(row,1);
	maxstep = min((maxbridges*25)+XDfVAc3Mmp,p+1);
	for jjj = 1:wibFptvq9Y
		[row2 col2] = find(B == -p);
		[row col] = find(B == p);
		connected = false;
		path = d4TZerDbiG(B,YTx56DA5U8,P9_FTG1hds,row2,col2,row,col, -p, maxbridges, XDfVAc3Mmp, p);
		if size(path,1) > 0
			W = [W; path];
			B = addwirepath(B,path,-p);
			hr2OdLy5K5();
			connected = true;
		end
		if ~connected
			break
		end
	end
end
end


function path = d4TZerDbiG(B,YTx56DA5U8,P9_FTG1hds,rowS,colS,FIg5RLlRYM,bX69T9QEhw,label,maxbridges,XDfVAc3Mmp,maxpathlen)
function path = traceback_a(Z,step)
GDcVqvkuIs(Z) = tQ2pJ8E0CI;
r=mod(Z,NR);
c=ceil(Z/NR);
path = zeros(step,4);
jjj = 0;
while isempty(find(rowS==r & colS==c,1))
	jjj = jjj + 1;
	path(jjj,1:2) = [r c];
	Z = r + (c-1)*NR;
	FFgpBpzZjs = GDcVqvkuIs(Z);
	preR=mod(FFgpBpzZjs,NR);
	preC=ceil(FFgpBpzZjs/NR);
	path(jjj,3:4) = [preR preC];
	r = preR;
	c = preC;
	if C1o3L2gfDa(r,c)
		jjj = jjj + 1;
		path(jjj,:) = [r c r c];
	end
end
path = path(1:jjj,:);
end
[NR NC] = size(B);
C1o3L2gfDa = false(NR,NC);
GDcVqvkuIs = zeros(NR,NC);
C = -ones(NR,NC);
C(rowS + (colS-1)*NR) = 0;
C(FIg5RLlRYM + (bX69T9QEhw-1)*NR) = -2;
maxstep = min((maxbridges*25)+XDfVAc3Mmp,maxpathlen+1);
kpmgLVupGj = zeros(maxstep+26,1);
for step = 0:maxstep
	if step == 0
		TQXX4Ib2IQ = rowS + (colS-1)*NR;
	elseif kpmgLVupGj(step) == 0
		continue
	else
		TQXX4Ib2IQ = find(C == step);
	end
	N = numel(TQXX4Ib2IQ);
	for i = 1:N
		tQ2pJ8E0CI = TQXX4Ib2IQ(i);
		Z = tQ2pJ8E0CI - NR;
		tag = C(Z);
		if tag == -2
			path = traceback_a(Z,step+1);
			return
		elseif tag == -1
			if B(Z) == 0
				C(Z) = step+1; kpmgLVupGj(step+1) = 1;
				GDcVqvkuIs(Z) = tQ2pJ8E0CI;
			elseif YTx56DA5U8(Z)
				C(Z) = step+26; kpmgLVupGj(step+26) = 1;
				GDcVqvkuIs(Z) = tQ2pJ8E0CI;
				C1o3L2gfDa(Z) = true;
			end
		end
		Z = tQ2pJ8E0CI + NR;
		tag = C(Z);
		if tag == -2
			path = traceback_a(Z,step+1);
			return
		elseif tag == -1
			if B(Z) == 0
				C(Z) = step+1; kpmgLVupGj(step+1) = 1;
				GDcVqvkuIs(Z) = tQ2pJ8E0CI;
			elseif YTx56DA5U8(Z)
				C(Z) = step+26; kpmgLVupGj(step+26) = 1;
				GDcVqvkuIs(Z) = tQ2pJ8E0CI;
				C1o3L2gfDa(Z) = true;
			end
		end
		Z = tQ2pJ8E0CI - 1;
		tag = C(Z);
		if tag == -2
			path = traceback_a(Z,step+1);
			return
		elseif tag == -1
			if B(Z) == 0
				C(Z) = step+1; kpmgLVupGj(step+1) = 1;
				GDcVqvkuIs(Z) = tQ2pJ8E0CI;
			elseif P9_FTG1hds(Z)
				C(Z) = step+26; kpmgLVupGj(step+26) = 1;
				GDcVqvkuIs(Z) = tQ2pJ8E0CI;
				C1o3L2gfDa(Z) = true;
			end
		end
		Z = tQ2pJ8E0CI + 1;
		tag = C(Z);
		if tag == -2
			path = traceback_a(Z,step+1);
			return
		elseif tag == -1
			if B(Z) == 0
				C(Z) = step+1; kpmgLVupGj(step+1) = 1;
				GDcVqvkuIs(Z) = tQ2pJ8E0CI;
			elseif P9_FTG1hds(Z)
				C(Z) = step+26; kpmgLVupGj(step+26) = 1;
				GDcVqvkuIs(Z) = tQ2pJ8E0CI;
				C1o3L2gfDa(Z) = true;
			end
		end
	end
end
path = zeros(0,4);
end

function path = DqbHRtn5ft(B,rowS,colS,FIg5RLlRYM,bX69T9QEhw,label,cutoff,maxpathlen)
function path = traceback_a(r,c,pathLength)
pR(r,c) = zR(i);
pC(r,c) = zC(i);
path = zeros(pathLength,4);
for jjj = 1:pathLength
	path(jjj,1:2) = [r c];
	preR = pR(r,c);
	preC = pC(r,c);
	path(jjj,3:4) = [preR preC];
	r = preR;
	c = preC;
end
end
[NR NC] = size(B);
pR = zeros(NR,NC);
pC = zeros(NR,NC);
C = -ones(NR,NC);
C(rowS+(colS-1)*NR) = 0;
C(FIg5RLlRYM+(bX69T9QEhw-1)*NR) = -2;
znext = zeros(NR*NC,1);
cnext = zeros(NR*NC,1);
count = numel(rowS);
znext(1:count) = rowS;
cnext(1:count) = colS;
dR=[-1 1 0 0];
dC=[0 0 -1 1];
for step = 0:min(cutoff,maxpathlen)
	if count < 1, break, end
	N = count;
	zR = znext(1:N);
	zC = cnext(1:N);
	count = 0;
	for i = 1:N
		Aa2p_xqlDq = zC(i);
		zi = zR(i);
		for s=1:4
			r = zi + dR(s);
			c = Aa2p_xqlDq + dC(s);
			Z = r + (c-1)*NR;
			tag = C(Z);
            if tag == -1 && B(Z) == 0
				C(Z) = step+1;
				pR(Z) = zi;
				pC(Z) = Aa2p_xqlDq;
				count = count + 1; znext(count) = r; cnext(count) = c;
            elseif tag == -2
				path = traceback_a(r,c,step+1);
				return

			end
		end
	end
end
path = [];
end