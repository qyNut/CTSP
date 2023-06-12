function reconImage = iterateFunc(dataloss, fdctPara, iterPara)

%% initial Curvelet Transform
C = fdct(dataloss,fdctPara); 
nscales = length(C);

%% Find the first several largest curvelet coefficients
coeff = 0;
for i = 1:nscales
    for j = 1:length(C{i})
        coeff = max(coeff,max(C{i}{j}(:)));
    end
end

m = mean(dataloss,'all');

L1 = iterPara.fitL1.a*m^iterPara.fitL1.b+iterPara.fitL1.c;
L2 = iterPara.fitL2.a*m^iterPara.fitL2.b+iterPara.fitL2.c;

%% Define regularization parameter
x = C;
u = dataloss;
M=zeros(size(u));
M(dataloss~=0)=1;

counter = 1;
outerloops = iterPara.outerloops;
innerloops = iterPara.innerloops;

CurDecR = L1*(1:-1/outerloops:1/outerloops);
SpatDecR = L2*(1/outerloops:1/outerloops:1);

lambda = coeff(1)*CurDecR(1);
lambda2 = max(dataloss(:))*SpatDecR(1);
mu = iterPara.mu;

%% Start iteration
while counter <= outerloops
    disp(['loop ',num2str(counter)]);
    
    for i = 1:innerloops 
        coeff = 0;
        temp1 = ifdct(x,fdctPara).*M;
        temp2 = dataloss-temp1;
        temp3 = fdct((temp2.*M),fdctPara);
        
        for j = 1:nscales
            for k = 1:length(C{j}) 
                dummy = x{j}{k} + temp3{j}{k};
                x{j}{k} = sign(dummy).*max(0,abs(dummy) - abs(lambda));
                coeff = max(coeff,max(x{j}{k}(:)));
            end
        end
        
        u = ifdct(x,fdctPara);
        temp4 = u-M.*u;
        u = M.*u+sign(temp4).*max(0,abs(temp4) - abs(lambda2));
        
        [ux,uy] = gradient(u);
        tot = sqrt(sum(ux.^2+uy.^2,'all'));
        [uxx,~] = gradient(ux./tot);
        [~,uyy] = gradient(uy./tot);
        TV = uxx+uyy;
        u = u-mu*TV;
            
        x = fdct(u,fdctPara);
        
    end
    coeff = sort(coeff,'descend');
    lambda = coeff(1)*CurDecR(counter);
    lambda2 = max(u(:))*SpatDecR(counter);
    counter = counter + 1;

end

reconImage = abs(ifdct(x,fdctPara));




