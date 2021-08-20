function [p1onsets] = seqofseq(X,lag)
% Adapted from Yunzhe's code

%% Sequence of Sequence


shifted_forward=X((1+lag):end,:);    % move up
shifted_backward=[zeros(lag,nStates);X(1:end-2*lag,:)];      % move down

p1onsets = zeros(nStates,nStates,size(X,1));
p1offsets = zeros(nStates,nStates,size(X,1));

p2onsets=zeros(nStates,nStates,length(X));
p2offsets=zeros(nStates,nStates,length(X));

for i=1:nStates
    for j=1:nStates
        p1onsets(i,j,1:(end-lag))=X(1:(end-lag),i).*shifted_forward(:,j);
        p2onsets(i,j,1:(end-lag))=X(1:(end-lag),i).*shifted_forward(:,j);
        
        p1offsets(i,j,1:(end-lag))=X(1:(end-lag),j).*shifted_backward(:,i); 
        p2offsets(i,j,1:(end-lag))=X(1:(end-lag),j).*shifted_backward(:,i); 
    end
end

%% PolicySeq -> NonPolicySeq
[I,J] = find(T1);
p1transitions=[I,J];

[I,J] = find(T2);
p2transitions = [I,J];

p2idx = find(T2);
YM = reshape(p2onsets, [nStates*nStates,length(p2onsets)]);
X = YM(p2idx,:)';

p1betas = nan(size(p1transitions,1),length(p2idx),maxLag+1);
p2betas = nan(size(p1transitions,1),length(p2idx),maxLag+1);

for path = 1:size(p1transitions,1) % number of sequences of interest
    
    ii = p1transitions(path,1);
    jj = p1transitions(path,2);
    
    Yinterest = squeeze(p1offsets(ii,jj,:));
    
    BP = zeros(1,size(Y,1));  %% shifted states
    BP(1:(end-lag)) = Y(1:(end-lag),jj); %% shifted states
    
    AP = zeros(1,size(Y,1));  %% shifted states
    AP(1:(end-lag)) = shifted_backward(:,ii); %% shifted states
    
    R = [Yinterest,AP',BP'];
    
    %% Core Sequenceness

    warning off
    dm=[toeplitz(R(:,1),[zeros(nbins,1)])];
    dm=dm(:,2:end);
    for kk=2:size(R,2)
        tmp=toeplitz(R(:,kk),[zeros(nbins,1)]);
        tmp=tmp(:,2:end);
        dm=[dm tmp];
    end
    warning on
    
    for ilag=1:bins
        
        idx = (1:bins:size(R,2)*maxLag) + ilag - 1;
        
        for st=1:size(Y,2)
            
            yi=p2transitions(st,1);
            yj=p2transitions(st,2);
            
            YOI1=squeeze(p2onsets(yi,:,:))';
            YOI2=squeeze(p2onsets(:,yj,:))';
            
            YOI1(:,yj)=[];
            YOI2(:,yi)=[];
            YR=squeeze(p2onsets(yi,yj,:));
            
            YOI=[YR,YOI1,YOI2];
            
            temp1 = pinv([dm(:,idx) ones(length(dm(:,idx)),1)])*YOI;
            
            TOI=zeros(size(YOI,2),1);
            TOI(1)=1;
            tempbeta=pinv([TOI,ones(size(YOI,2),1)])*temp1(1,:)';
            p1betas(path,st,ilag+1)= tempbeta(1);
            
            %                         YA=zeros(1,length(X));
            %                         YA(1:(end-NPlag))=X(1:(end-NPlag),yi);
            %                         temp3 = pinv([dm(:,temp_zinds) ones(length(dm(:,temp_zinds)),1)])*YA';
            %
            %                         YB=zeros(1,length(X));
            %                         YB(1:(end-Plag))=Xnplag_onset(:,yj);
            %                         temp4 = pinv([dm(:,temp_zinds) ones(length(dm(:,temp_zinds)),1)])*YB';
            %                         betaSeq1(iseq,iY,ilag+1)= temp3(1)*temp4(1);
        end
        
        tmp = pinv([dm(:,idx) ones(length(dm(:,idx)),1)])*Y;
        p2betas(path,:,ilag+1)=tmp(1,:);
    end
end

beta_p1 = reshape(p1betas,[size(p1betas,1)*size(p1betas,2),size(p1betas,3)]);
betaM_P2P_1 = nanmean(beta_p1,1);

beta_p2 = reshape(p2betas,[size(p2betas,1)*size(p2betas,2),size(p2betas,3)]);
betaM_P2P_2 = nanmean(beta_p2,1);

%% Get path 1 to predict path 2

[I,J] = find(T2);
p1transitions = [I J];
p2idx = find(T1);

[I,J]=find(T1);
p2transitions = [I J];

YM=reshape(p1onsets, [nStates*nStates,length(p1onsets)]);
X = YM(p2idx,:)';

p1betas=nan(size(p1transitions,1),length(p2idx),maxLag+1);
p2betas=nan(size(p1transitions,1),length(p2idx),maxLag+1);

for path=1:size(p1transitions,1) % number of sequence of interest
    ii=p1transitions(path,1);
    jj=p1transitions(path,2);
    
    Yinterest=squeeze(p2offsets(ii,jj,:));
    
    BP=zeros(1,length(X));  %% shifted states
    BP(1:(end-lag))=X(1:(end-lag),jj); %% shifted states
    
    AP=zeros(1,length(X));  %% shifted states
    AP(1:(end-lag))= Xnplag_end(:,ii); %% shifted states
    
    R=[Yinterest,AP',BP'];
    %                 Rmatrix=[Xinterest,BP'];
    %                 Rmatrix=[Xinterest,AP'];
    
    %                 Rmatrix=[Xinterest,Xcontrol];
    %                 Rmatrix=[Xinterest,nanmean([XcontrolA'],2),nanmean([XcontrolB'],2)];
    
    %% Core Sequenceness
    nbins=maxLag+1;
    
    warning off
    dm=[toeplitz(R(:,1),[zeros(nbins,1)])];
    dm=dm(:,2:end);
    
    for kk=2:size(R,2)
        tmp=toeplitz(R(:,kk),[zeros(nbins,1)]);
        tmp=tmp(:,2:end);
        dm=[dm tmp];
    end
    warning on
    
    bins=maxLag;
    
    for ilag=1:bins
        idx = (1:bins:size(R,2)*maxLag) + ilag - 1;
        
        for st=1:size(X,2)
            
            yi=p2transitions(st,1);
            yj=p2transitions(st,2);
            
            YOI1=squeeze(p1onsets(yi,:,:))';
            YOI2=squeeze(p1onsets(:,yj,:))';
            
            YOI1(:,yj)=[];
            YOI2(:,yi)=[];
            YR=squeeze(p1onsets(yi,yj,:));
            
            YOI=[YR,YOI1,YOI2];
            
            temp1 = pinv([dm(:,idx) ones(length(dm(:,idx)),1)])*YOI;
            
            TOI=zeros(size(YOI,2),1);
            TOI(1)=1;
            tempbeta=pinv([TOI,ones(size(YOI,2),1)])*temp1(1,:)';
            p1betas(path,st,ilag+1)= tempbeta(1);
            
            %                         YA=zeros(1,length(X));
            %                         YA(1:(end-Plag))=X(1:(end-Plag),yi);
            %                         temp3 = pinv([dm(:,temp_zinds) ones(length(dm(:,temp_zinds)),1)])*YA';
            %
            %                         YB=zeros(1,length(X));
            %                         YB(1:(end-Plag))=Xplag_onset(:,yj);
            %                         temp4 = pinv([dm(:,temp_zinds) ones(length(dm(:,temp_zinds)),1)])*YB';
            %
            %                         betaSeq1(iseq,iY,ilag+1)= temp3(1)*temp4(1);
        end
        
        
        tmp = pinv([dm(:,idx) ones(length(dm(:,idx)),1)])*X;
        p2betas(path,:,ilag+1)=tmp(1,:);
    end
end

beta_p1=reshape(p1betas,[size(p1betas,1)*size(p1betas,2),size(p1betas,3)]);
betaM_NP2P_1=nanmean(beta_p1,1);
%             betaM_NP2P(2:end)=detrend(betaM_NP2P(2:end));

beta_p2=reshape(p2betas,[size(p2betas,1)*size(p2betas,2),size(p2betas,3)]);
betaM_NP2P_2=nanmean(beta_p2,1);

sp1{iSj}(1,itrial,iShuf,:) = betaM_P2P_1;  % policy->non-policy, excluding zero-lag correlation for now
sn1{iSj}(1,itrial,iShuf,:) = betaM_NP2P_1;  % non-policy->policy, excluding zero-lag correlation for now

sp2{iSj}(1,itrial,iShuf,:) = betaM_P2P_2;  % policy->non-policy, excluding zero-lag correlation for now
sn2{iSj}(1,itrial,iShuf,:) = betaM_NP2P_2;  % non-policy->policy, excluding zero-lag correlation for now

end