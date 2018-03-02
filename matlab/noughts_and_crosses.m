function noughts_and_crosses(xind,yind,varargin)
figure;
if nargin >2
	if nargin  ~= 4
		error('If nargin >2 it should equal 4');
		return;
	end
	xs=varargin{1};
	ys=varargin{2};
else
	xs=[300 500 700];
	ys=[5 10 15];
end


[Xs,Ys]=meshgrid(xs,ys);

% nested for loop
for i=1:length(xind)
	for j=1:length(yind)
		Ys(yind(j),xind(i))=NaN;
	end
end

plot(Xs,Ys,'kx','markersize',20)
axis([0 800 0 20])

[Xs1,Ys1]=meshgrid(xs,ys);

Ys1(find(~isnan(Ys)))=NaN
hold on;
plot(Xs1,Ys1,'ko','markersize',20)
xlabel('Semi-major axis (AU)');
ylabel('Mass (Earth Masses)')
