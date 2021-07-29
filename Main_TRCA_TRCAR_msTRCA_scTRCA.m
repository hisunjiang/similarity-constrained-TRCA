
clear all;
close all;

%% Parameter setting
Fs=250;         % sample rate
ch_used = [48 54 55 56 57 58 61 62 63]; % Pz, PO5, PO3, POz, PO4, PO6, O1,Oz, O2 (in SSVEP benchmark dataset)
Nchannel=length(ch_used);

num_of_trials = 2;                  % Number of training trials
num_of_harmonics=5;                 % for all cca-based methods
TRCAwithRef=1;                      % For TRCA-R (1: with Ref., 0: without Ref.)
num_of_signal_templates2=2;         % for ms-etrca (1<=num_of_signal_templates<=40)  
dataLength = 0.7;                   % time-window length (s)
sig_len = round(dataLength*Fs); 
num_of_subbands=5;                  % for filter bank analysis
FB_coef0=[1:num_of_subbands].^(-1.25)+0.25; 
enable_bit=[1 1 1 1];               % Select the algorithms: bit 1: TRCA, bit 2: TRCAR, bit 3: msTRCA, bit 4: scTRCA
is_center_std=1;                    % 0: without , 1: with (zero mean, and unity standard deviation)

%% Filter bank
% Chebyshev Type I filter
for k=1:num_of_subbands
    bandpass1(1)=8*k;
    bandpass1(2)=90;
    [b2(k,:),a2(k,:)] = cheby1(4,1,[bandpass1(1)/(Fs/2) bandpass1(2)/(Fs/2)],'bandpass');
end

%% Calculation
seed = RandStream('mt19937ar','Seed','shuffle');
num_of_subj=35;                     % Number of subjects
accuracy = zeros(num_of_subj,4);
for sn=1:num_of_subj
    %% data preprocessing
    % Directory of the SSVEP Dataset
    str_dir='E:\laboratory\codes\DATA\SSVEP data\ssvep_benchmark_dataset'; 
    string=[str_dir, '\S', num2str(sn), '.mat'];
    load(string)
    
    %  pre-stimulus period: 0.5 sec
    %  latency period: 0.14 sec
    eeg=data(ch_used,floor(0.5*Fs+0.14*Fs):floor(0.5*Fs+0.14*Fs)+4*Fs-1,:,:);
    [d1_,d2_,d3_,d4_]=size(eeg);
    d1=d3_;d2=d4_;d3=d1_;d4=d2_;
    no_of_class=d1;
    % d1: num of stimuli
    % d2: num of blocks
    % d3: num of channels       
    % d4: num of sampling points
    for i=1:1:d1
        for j=1:1:d2
            y=reshape(eeg(:,:,i,j),d3,d4);
            SSVEPdata(:,:,j,i)=reshape(y,d3,d4,1,1);
            
            for sub_band=1:num_of_subbands
                
                for ch_no=1:d3
                    if (num_of_subbands==1)
                        y_sb(ch_no,:)=y(ch_no,:);
                    else
                        y_sb(ch_no,:)=filtfilt(b2(sub_band,:),a2(sub_band,:),y(ch_no,:));
                    end
                end
                subband_signal(sub_band).SSVEPdata(:,:,j,i)=reshape(y_sb,d3,d4,1,1);
            end
            
        end
    end
    
    clear eeg
    
    %% stimulation parameters initialization    
    pha_val=[0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 ...
        0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5 0 0.5 1 1.5]*pi;
%     pha_val=zeros(1,40);
    sti_f=[8:0.2:15.8];
    n_sti=length(sti_f);            % number of stimulus frequency
    temp=reshape([1:40],8,5);
    temp=temp';
    target_order=temp(:)';          % re-order the data, it's necessary if you do multi-stimulus learning 
    SSVEPdata=SSVEPdata(:,:,:,target_order);
    for sub_band=1:num_of_subbands
        subband_signal(sub_band).SSVEPdata=subband_signal(sub_band).SSVEPdata(:,:,:,target_order); % To sort the orders of the data as 8.0, 8.2, 8.4, ..., 15.8 Hz
    end
    
    FB_coef=FB_coef0'*ones(1,n_sti);
    n_correct=zeros(1,4);           % Count how many correctly detected trials for each method
       
    %% repeat several times of random validation
    CVtimes = d2;
    seq_0=zeros(CVtimes,num_of_trials);
    for run=1:CVtimes
        if (num_of_trials==1)
            seq1=run;
        elseif (num_of_trials==d2-1)
            seq1=[1:d2];
            seq1(run)=[];
        else
            % Randomly select num_of_trials blocks for training
            isOK=0;
            while (isOK==0)
                seq=randperm(seed,d2);
                seq1=seq(1:num_of_trials);
                seq1=sort(seq1);
                if isempty(find(sum((seq1'*ones(1,d2)-seq_0').^2)==0))
                    isOK=1;
                end
            end
        end
        idx_traindata=seq1;     % index of the training blocks
        idx_testdata=1:d2;   
        idx_testdata(seq1)=[];  % index of the testing blocks
        
        for i=1:no_of_class
            for k=1:num_of_subbands
                if length(idx_traindata)>1
                    subband_signal(k).signal_template(i,:,:)=mean(subband_signal(k).SSVEPdata(:,:,idx_traindata,i),3);
                else
                    subband_signal(k).signal_template(i,:,:)=subband_signal(k).SSVEPdata(:,:,idx_traindata,i);
                end
            end
        end
        
        %% test on each testing block
        for run_test=1:length(idx_testdata) % number of testing block
            test_signal=zeros(d3,sig_len);
            fprintf('No.crossvalidation: %d, No.test block: %d \n',run, run_test);

            for i=1:no_of_class     % each testing trial
                for sub_band=1:num_of_subbands
                    test_signal=subband_signal(sub_band).SSVEPdata(:,1:sig_len,idx_testdata(run_test),i);
                    if (is_center_std==1)
                        test_signal=test_signal-mean(test_signal,2)*ones(1,length(test_signal));
                        test_signal=test_signal./(std(test_signal')'*ones(1,length(test_signal)));
                    end
                    for j=1:no_of_class     % each stimulus
                        template=reshape(subband_signal(sub_band).signal_template(j,:,[1:sig_len]),d3,sig_len);
                        if (is_center_std==1)
                            template=template-mean(template,2)*ones(1,length(template));
                            template=template./(std(template')'*ones(1,length(template)));
                        end

                        % Generate the sine-cosine reference signal
                        ref1=ref_signal_nh(sti_f(j),Fs,pha_val(j),sig_len,num_of_harmonics);
                        
                        %===============eTRCA==================
                        if (enable_bit(1)==1)
                            if (num_of_trials==1)
                                TRCAR(sub_band,j)=0;
                            else
                                if ((i==1) && (j==1))
                                    W_TRCA(sub_band).val=[];
                                    for jj=1:no_of_class
                                        X1=[];
                                        TRCA_X1=zeros(d3,sig_len);
                                        for tr=1:num_of_trials
                                            X0=reshape(subband_signal(sub_band).SSVEPdata(:,1:sig_len,idx_traindata(tr),jj),d3,sig_len);
                                            if (is_center_std==1)
                                                X0=X0-mean(X0,2)*ones(1,length(X0));
                                                X0=X0./(std(X0')'*ones(1,length(X0)));
                                            end
                                            TRCA_X1=TRCA_X1+X0;
                                            X1=[X1;X0'];
                                        end
                                        S=TRCA_X1*TRCA_X1'-X1'*X1;
                                        Q=X1'*X1;
                                        [eig_v1,eig_d1]=eig(Q\S);
                                        [eig_val,sort_idx]=sort(diag(eig_d1),'descend');
                                        eig_vec=eig_v1(:,sort_idx);
                                        W_TRCA(sub_band).val=[W_TRCA(sub_band).val; eig_vec(:,1)'];
                                    end
                                end
                                cr1=corrcoef(W_TRCA(sub_band).val*test_signal,W_TRCA(sub_band).val*template);
                                TRCAR(sub_band,j)=cr1(1,2);
                            end
                        else
                            TRCAR(sub_band,j)=0;
                        end
                        
                        %===============eTRCA-R==================
                        if (enable_bit(2)==1)
                            if (num_of_trials==1)
                                TRCARR(j)=0;
                            else
                                if ((i==1) && (j==1))
                                    W_TRCAR(sub_band).val=[];
                                    for jj=1:no_of_class
                                        TRCA_X=[];
                                        if (TRCAwithRef==1)
                                            Ref=ref_signal_nh(sti_f(jj),Fs,pha_val(jj),sig_len,num_of_harmonics);
                                            [Q_ref1,R_ref1]=qr(Ref',0);
                                            ref_matrix=Q_ref1*Q_ref1';
                                        else
                                            ref_matrix=eye(sig_len);
                                        end
                                        LL=repmat(ref_matrix,num_of_trials);
                                        if (num_of_trials==5)
                                            LL=LL-blkdiag(ref_matrix,ref_matrix,ref_matrix,ref_matrix,ref_matrix);
                                        elseif (num_of_trials==4)
                                            LL=LL-blkdiag(ref_matrix,ref_matrix,ref_matrix,ref_matrix);
                                        elseif (num_of_trials==3)
                                            LL=LL-blkdiag(ref_matrix,ref_matrix,ref_matrix);
                                        elseif (num_of_trials==2)
                                            LL=LL-blkdiag(ref_matrix,ref_matrix);
                                        else
                                        end

                                        for tr=1:num_of_trials
                                            X0=reshape(subband_signal(sub_band).SSVEPdata(:,1:sig_len,idx_traindata(tr),jj),d3,sig_len);
                                            if (is_center_std==1)
                                                X0=X0-mean(X0,2)*ones(1,length(X0));
                                                X0=X0./(std(X0')'*ones(1,length(X0)));
                                            end
                                            TRCA_X=[TRCA_X;X0'];
                                        end
                                        S=TRCA_X'*LL*TRCA_X;
                                        Q=TRCA_X'*TRCA_X;
                                        [eig_v1,eig_d1]=eig(Q\S);
                                        [eig_val,sort_idx]=sort(diag(eig_d1),'descend');
                                        eig_vec=eig_v1(:,sort_idx);
                                        W_TRCAR(sub_band).val = [W_TRCAR(sub_band).val; eig_vec(:,1)'];
                                    end
                                end
                                cr1=corrcoef(W_TRCAR(sub_band).val*test_signal,W_TRCAR(sub_band).val*template);
                                TRCARR(sub_band,j)=cr1(1,2);
                            end
                        else
                            TRCARR(sub_band,j)=0;
                        end
                        
                        %===============emsTRCA==================
                        if (enable_bit(3)==1) 
                            if (num_of_trials==1)
                                msTRCAR(sub_band,j)=0;
                            else
                                if ((i==1) && (j==1))
                                    W_msTRCA(sub_band).val=[];
                                    for jj=1:no_of_class
                                        d0=floor(num_of_signal_templates2/2);
                                        if jj<=d0
                                            template_st=1;
                                            template_ed=num_of_signal_templates2;
                                        elseif ((jj>d0) && jj<(d1-d0+1))
                                            template_st=jj-d0;
                                            template_ed=jj+(num_of_signal_templates2-d0-1);
                                        else
                                            template_st=(d1-num_of_signal_templates2+1);
                                            template_ed=d1;
                                        end
                                        template_seq=[template_st:template_ed];
                                        msTRCA_X1=[];
                                        msTRCA_X2=[];

                                        for n_temp=1:num_of_signal_templates2
                                            jj=template_seq(n_temp);
                                            TRCA_X1=zeros(d3,sig_len);
                                            template2=zeros(d3,sig_len);
                                            TRCA_X2=[];
                                            for tr=1:num_of_trials
                                                X0=reshape(subband_signal(sub_band).SSVEPdata(:,1:sig_len,idx_traindata(tr),jj),d3,sig_len);
                                                if (is_center_std==1)
                                                    X0=X0-mean(X0,2)*ones(1,length(X0));
                                                    X0=X0./(std(X0')'*ones(1,length(X0)));
                                                end
                                                TRCA_X2=[TRCA_X2;X0'];
                                                TRCA_X1=TRCA_X1+X0;
                                            end
                                            msTRCA_X1=[msTRCA_X1 TRCA_X1];
                                            msTRCA_X2=[msTRCA_X2 TRCA_X2'];
                                        end
                                        S=msTRCA_X1*msTRCA_X1'-msTRCA_X2*msTRCA_X2';
                                        Q=msTRCA_X2*msTRCA_X2';
                                        [eig_v1,eig_d1]=eig(Q\S);
                                        [eig_val,sort_idx]=sort(diag(eig_d1),'descend');
                                        eig_vec=eig_v1(:,sort_idx);
                                        W_msTRCA(sub_band).val=[W_msTRCA(sub_band).val; eig_vec(:,1)'];
                                    end
                                end
                                cr1=corrcoef(W_msTRCA(sub_band).val*test_signal, W_msTRCA(sub_band).val*template);
                                msTRCAR(sub_band,j)=cr1(1,2);
                            end
                        else
                            msTRCAR(sub_band,j)=0;
                        end

                        %===============escTRCA==================
                        if (enable_bit(4)==1)
                            if (num_of_trials==1)
                                scTRCAR(sub_band,j)=0;
                            else
                                if ((i==1) && (j==1))
                                    for jj=1:no_of_class
                                        % current reference
                                        ref1_now=ref_signal_nh(sti_f(jj),Fs,pha_val(jj),sig_len,num_of_harmonics);
                                        % scTRCA
                                        scTRCA_X1=[];
                                        scTRCA_X2=[];
                                        for tr=1:num_of_trials
                                            X0=reshape(subband_signal(sub_band).SSVEPdata(:,1:sig_len,idx_traindata(tr),jj),d3,sig_len);
                                            if (is_center_std==1)
                                                X0=X0-mean(X0,2)*ones(1,length(X0));
                                                X0=X0./(std(X0')'*ones(1,length(X0)));
                                            end
                                            scTRCA_X1=[scTRCA_X1,X0];
                                            scTRCA_X2(:,:,tr)=X0;
                                        end
                                        Xa{1}=scTRCA_X1;
                                        Xa{2}=ref1_now;
                                        Xb{1}=scTRCA_X2;
                                        Xb{2}=ref1_now;
                                        % wn
                                        wn = scTRCA(Xa, Xb);
                                        W_scTRCA(sub_band).w1(:, jj) = wn(1:Nchannel, :);
                                        W_scTRCA(sub_band).w2(:, jj) = wn(Nchannel+1:end, :);
                                    end
                                end
                                cr1=corrcoef((W_scTRCA(sub_band).w1' * test_signal),(W_scTRCA(sub_band).w2' * ref1));
                                cr2=corrcoef((W_scTRCA(sub_band).w1' * test_signal),(W_scTRCA(sub_band).w1' * template));
                                scTRCAR(sub_band,j)=sign(cr1(1,2))*cr1(1,2)^2+sign(cr2(1,2))*cr2(1,2)^2;
                            end
                        else
                            scTRCAR(sub_band,j)=0;
                        end
                    end

                end

                TRCAR1=sum((TRCAR).*FB_coef,1);
                TRCARR1=sum((TRCARR).*FB_coef,1);
                msTRCAR1=sum((msTRCAR).*FB_coef,1);
                scTRCAR1=sum((scTRCAR).*FB_coef,1);

                [~,idx]=max(TRCAR1);
                if idx==i
                    n_correct(1)=n_correct(1)+1;
                end
                [~,idx]=max(TRCARR1);
                if idx==i
                    n_correct(2)=n_correct(2)+1;
                end
                [~,idx]=max(msTRCAR1);
                if idx==i
                    n_correct(3)=n_correct(3)+1;
                end
                [~,idx]=max(scTRCAR1);
                if idx==i
                    n_correct(4)=n_correct(4)+1;
                end
            end
        end
        seq_0(run,:)=seq1;
    end

    %% Save results
    accuracy(sn,:)=100*n_correct/n_sti/CVtimes/length(idx_testdata);
    fprintf('Subject:%d, Data length:%.2fs, accuracy:%.4f, %.4f, %.4f, %.4f,\n', ...
        sn, dataLength, accuracy(sn,1), accuracy(sn,2), accuracy(sn,3), accuracy(sn,4));
    % column 1: TRCA
    % column 2: TRCA-R
    % column 3: msTRCA
    % column 4: scTRCA
end


