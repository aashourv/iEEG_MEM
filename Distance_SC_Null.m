%  ========================================================================
%  This script generates structural connectivity nulls for a subset of 
%  anatomical regions of interests (ROIs). In this script the empirical 
%  subset of ROIs correspond to ROIs closest to the sample patinet's 
%  intra-cranial EEG (iEEG) electrodes. This algorithm creates structural 
%  nulls by iteratively resampling ROIs aimed at keeping the pairwise 
%  distance between the null ROI pairs closest to the empirical distance 
%  profiles. "NULL_ROIS" and "Null_Dist" variables contain the null ROI 
%  labels (i.e., numbers) and the distance between null ROI pairs, 
%  respectively.
% 
%  
%  Arian Ashourvan, University of Pennsylvania, july 2020 
%  ========================================================================
%
%% Uncomment if running on cluster (e.g., using qsub command)
rr = runNo;
%% Comment if running on cluster (e.g., using qsub command)
% rr=1;
%% Add Nifti toolbox to path
addpath(genpath('/data/jux/aashourvan/InterictalVR/Cluster_Codes/Comm_Bio_Codes/NIfTI_20140122'));
%% Select a patient for analysis
Pat_Name='HUP094';
outpath=['/data/jux/aashourvan/InterictalVR/Results'];
mkdir(outpath);
%% Load Patient's iEEG electrode (i.e, channel) labels
load('/data/jux/aashourvan/InterictalVR/Cluster_Codes/Comm_Bio_Codes/All_Pat_Channels.mat')
channels= All_Pat_channels{1};
%% Remove bad channels
channels=channels(3:end);
channels=channels([1:10 ,13:end]);
%% Load anatomical parcelations (AAL 600) and find the average MNI
%  coordinate of every region of interest (ROI)
NII=load_nii('/data/jux/aashourvan/InterictalVR/Cluster_Codes/Comm_Bio_Codes/AAL3_1mm_flipx.nii.gz')'; % AAL 600 1mm
Num_Anat_ROIs=600;
for ROI =1 :Num_Anat_ROIs
    [r,c,v] = ind2sub(size(NII.img),find(NII.img == ROI));
    CORDIN(ROI,1)=mean(r);
    CORDIN(ROI,2)=mean(c);
    CORDIN(ROI,3)=mean(v);
end
%% Load MNI coordinates of all channels
filename = ['/data/jux/aashourvan/InterictalVR/Cluster_Codes/Comm_Bio_Codes/centroids_mni_' Pat_Name '.csv'];
delimiter = ',';
formatSpec = '%f%f%f%f%s%[^\n\r]';
fileID = fopen(filename,'r');
%dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
fclose(fileID);
centroidsmni = table(dataArray{1:end-1}, 'VariableNames', {'VarName1','X','Y','Z','Label'});
clearvars filename delimiter formatSpec fileID dataArray ans;

%% Find the Euclidean distance between all possible pairs of ROIs and channels
for elect = 1:length(channels)
    opaa=find(ismember(centroidsmni.Label,channels{elect}));
    for ROI=1 :Num_Anat_ROIs
        DISTT(elect,ROI)=pdist([CORDIN(ROI,:)',[centroidsmni{opaa,2:4}]']');
    end
end
%% Find the closest ROI to each electrode
imagesc(DISTT)
hold on
for elect = 1:length(channels)
    
    plot(find(DISTT(elect,:)==min(DISTT(elect,:))), elect, 'r*')
    
    ELect_ROI(elect)=find(DISTT(elect,:)==min(DISTT(elect,:)));
    ELect_ROI_Dist(elect)= min(DISTT(elect,:));
    
end
ELect_ROI_orig=ELect_ROI;

% For ROIs with more than one closest electrode only pick the closet one
% to  ROIs' centroid
[ aa, bb,cc]= unique(ELect_ROI);
clear Inddexx
for jj = 1: length(aa)
    ttteempp=find(ELect_ROI==aa(jj));
    dissttt=ELect_ROI_Dist(ELect_ROI==aa(jj));
    Inddexxjj=ttteempp(find(dissttt== min(dissttt)));
    Inddexx(jj)=Inddexxjj(1);
    clear ttteempp dissttt Inddexxjj
end
%% Find the Euclidean distance between all electrode pairs
for elect1 = 1:length(channels)
    opaa=find(ismember(centroidsmni.Label,channels{elect1}));
    for elect2=1:length(channels)
        opaa2=find(ismember(centroidsmni.Label,channels{elect2}));
        DISTT_pairs(elect1,elect2)=pdist([[centroidsmni{opaa2,2:4}]',[centroidsmni{opaa,2:4}]']');
        
    end
end
Emp_Dist=DISTT_pairs(Inddexx,Inddexx).*~eye(length(Inddexx));
Emp_Dist=Emp_Dist(find(triu(ones(size(Emp_Dist)),1)));
Emp_Dist2=DISTT_pairs(Inddexx,Inddexx).*~eye(length(Inddexx));
Emp_Dist22= triu(Emp_Dist2,1);
%%  List the electrodes that should to be prioritized by the algorithm
Special_elect_all=[];
%%  Initialize
iterations=10000;  % Number of Nulls ROI sets 
ROI_thresh=0;      % Use higher numbers if you dont want any of the N (default=0) selected 
                   % Null ROIs after the first pair of null ROIs to contain
                   % any of the empirical ROIs.
repeat_ratio=0.1;  % Accpeted ratio of the same null ROI selected for an electrode across iterations 
Pair_Dist_all_rois=zeros(Num_Anat_ROIs,Num_Anat_ROIs);
FULL_pairs=zeros(3,1);
qqq=1;
%% Find the distance between all ROI pairs
for roi1 = 1: Num_Anat_ROIs
    for roi2 = 1: Num_Anat_ROIs
        if roi1< roi2
            Pair_Dist_all_rois(roi1,roi2)=pdist([CORDIN(roi1,:)',CORDIN(roi2,:)']');
        end
    end  
end
Pair_Dist_all_rois=Pair_Dist_all_rois+Pair_Dist_all_rois';
%% Check for saved nulls
if exist([ outpath '/' Pat_Name '_Geometric_SC_Null_1.mat'], 'file') == 2
    load([ outpath '/' Pat_Name '_Geometric_SC_Null_1'],'NULL_ROIS','Null_Dist');
    NULL_ROIS=NULL_ROIS(:,~isnan(sum(NULL_ROIS,1)));
    Null_Dist=Null_Dist(:,:,~isnan(sum(NULL_ROIS,1)));
    NULL_ROIS2=nan(length(Inddexx),iterations);
    iter= size(NULL_ROIS,2);
    NULL_ROIS2(:,1:iter)=NULL_ROIS;
    clear NULL_ROIS;
    NULL_ROIS=NULL_ROIS2;
    clear NULL_ROIS2;
else
    iter=0;
    NULL_ROIS=nan(length(Inddexx),iterations);
end
%% Uncomment if running on cluster to change the seed for randi.m function
pause(rr)
rng('shuffle');
%%
for Iter = 1:iterations
   %%
    iter=iter+1;
    flag=0;
    while flag==0
        % Select a random pair of ROIs
        ROIs=randi(Num_Anat_ROIs,1,2);
        Pair_Dist=pdist([CORDIN(ROIs(1),:)',CORDIN(ROIs(2),:)']');
        % Accpet if the selectef ROI pair's distance is bigger than minimum
        % distance between all electrode pairs
        if ~isempty(find((Emp_Dist-Pair_Dist)<0))
            % Accpet if the distance between the selected ROI pair is between
            % -15 and 10 mm from the maximun distance between electode pairs
            if Pair_Dist > (max(Emp_Dist)-15) && Pair_Dist <= (max(Emp_Dist)+10)
                % Find the electrode pair that the distance between them is
                % closest to the distance between the selected ROI pair
                [nn mm]=find(abs(Emp_Dist22-Pair_Dist)==min(abs(Emp_Dist-Pair_Dist)));
                if length(nn)>1
                    lol=randi(length(nn),1,1);
                    nn=[nn(lol);mm(lol)];
                else
                    nn=[nn(1);mm(1)];
                end
                % Assign the electrodes the ROI pair's numbers in NULL_ROIS
                % variable
                if rand(1)>0.5
                    NULL_ROIS(nn,iter)=ROIs;
                else
                    NULL_ROIS(nn,iter)=fliplr(ROIs);
                end
                selcted_elect_index=nn;
                Null_Dist(nn(1),nn(2),iter)=Pair_Dist;
                Null_Dist(nn(2),nn(1),iter)=Pair_Dist;
                flag=1;
                % Accept if any of the selected ROIs are the same as the empirical
                % ROIs
                if ~isempty(find(ismember(NULL_ROIS(nn,iter),ELect_ROI(Inddexx))))
                    flag=0;
                end
                % Accept of non of the null ROIs are repeating more than
                % the predefined repeat_ratio
                if iter >1
                    DIF1= bsxfun(@minus,NULL_ROIS(nn,1:iter-1),NULL_ROIS(nn,iter));
                    DIF1(DIF1~=0)=1;
                    if ~isempty(find((sum(DIF1==0,2)/iter)>repeat_ratio))
                        flag=0;
                    end
                end
                if flag ==0
                    NULL_ROIS(nn,iter)=nan;
                end
            end
        end
    end
    % Check if the selcted null ROIs are found in previous iterations
    DIF1=bsxfun(@minus,NULL_ROIS(nn,1:iter-1),NULL_ROIS(nn,iter));
    DIF1(DIF1~=0)=1;
    if ~isempty(find(sum(DIF1,1)==0))
        flag_exist=1;
    else
        flag_exist=0;
    end    
    %%  Find ROI/s with closest distance profiles to the prioritized electrode/s
    Special_elect= Special_elect_all(~ismember(Special_elect_all,selcted_elect_index'));
    if ~isempty(Special_elect)
        for kkk= 1: length(Special_elect)
            clear Pair_Dist_roi
            tic
            Pair_Dist_roi=Pair_Dist_all_rois(:,NULL_ROIS(selcted_elect_index,iter));
            Euclidean_dis=nan(1,Num_Anat_ROIs);
            clear Null_temp_dist
            for roi = 1:Num_Anat_ROIs
                if ~ismember(roi,NULL_ROIS(:,iter))
                    TTRRM=Null_Dist(selcted_elect_index,selcted_elect_index,iter);
                    Null_temp_dist = [Pair_Dist_roi(roi,:) , TTRRM(logical(triu(ones(size(TTRRM)),1)))'  ];
                    Emp_Dist3=Emp_Dist2( [selcted_elect_index ;Special_elect(kkk) ],[selcted_elect_index ;Special_elect(kkk) ]);
                    Euclidean_dis(roi)=pdist([Null_temp_dist;Emp_Dist3(logical(triu(ones(size(Emp_Dist3)),1)))'],'Euclidean');
                end
            end
            [ooii ppii]=find(Euclidean_dis==min(Euclidean_dis(:)));
            clear  Euclidean_dis
            NULL_ROIS(Special_elect(kkk),iter)=ppii;
            Null_Dist(Special_elect(kkk),selcted_elect_index,iter)=Pair_Dist_roi(ppii,:);
            Null_Dist(selcted_elect_index,Special_elect(kkk),iter)=Pair_Dist_roi(ppii,:)';
            selcted_elect_index= [selcted_elect_index ;Special_elect(kkk) ];
            clear ooii ppii  Pair_Dist_roi
        end
    end
    %% From the remaining ROIs and electrodes find one ROI and one electrode 
    %  with closest distance profiles to the selected null ROIs and thier 
    %  corresponding electrodes, respectively. Repeat untill all electrodes
    %  are assigned a null ROI number.
    for ROI_Elect = 1:length(Inddexx)-(2+length(Special_elect))
        Iter
        iter
        ROI_Elect
        clear Pair_Dist_roi
        tic
        Pair_Dist_roi=Pair_Dist_all_rois(:,NULL_ROIS(selcted_elect_index,iter));
        OTHER_elect=find(~ismember(1:length(Inddexx),selcted_elect_index));
        Euclidean_dis=nan(length(OTHER_elect),Num_Anat_ROIs);
        clear Null_temp_dist
        for other_elect =1:length(OTHER_elect)
            for roi = 1:Num_Anat_ROIs
                if ~ismember(roi,NULL_ROIS(:,iter))
                    Null_temp_dist = [Pair_Dist_roi(roi,:) ];
                    Emp_Dist3=Emp_Dist2( OTHER_elect(other_elect),selcted_elect_index);
                    Euclidean_dis(other_elect,roi)=pdist([Null_temp_dist;Emp_Dist3],'Euclidean');
                end
            end
        end
        Euclidean_dis2=Euclidean_dis;
        [ooii ppii]=find(Euclidean_dis==min(Euclidean_dis(:)));
        if flag_exist ==0
        else % If the selected null ROIs exist in previous iterations
             % select the next best match
            flagg22=0 ;
            while  flagg22 ==0
                DIF2= bsxfun(@minus,FULL_pairs(1:2,:),nn);
                DIF2(DIF2~=0)=1;
                if ~isempty(find(sum(DIF2,1)==0)) && ROI_Elect == FULL_pairs(3,find(sum(DIF2,1)==0))
                    KK=2;
                else
                    KK=1;
                end
                for plp =1: KK
                    Euclidean_dis2(ooii,ppii)=nan;
                    [ooii ppii]=find(Euclidean_dis2==min(Euclidean_dis2(:)));
                end
                flagg22=1;
                if ~isempty(find(ismember(ppii,ELect_ROI(Inddexx)))) &&  ROI_Elect< ROI_thresh
                    flagg22=0 ;
                end
            end
        end
        if iter >1
            NULL_ROIS_temp=NULL_ROIS(:,iter);
            NULL_ROIS_temp(OTHER_elect(ooii))=ppii;
            Index1=find(~isnan(NULL_ROIS_temp));
            DIF1=bsxfun(@minus,NULL_ROIS(Index1,1:iter-1),NULL_ROIS_temp(Index1));
            DIF1(DIF1~=0)=1;
            %  If the selected null ROIs are found in previous iteration with
            %  rates higher than repeat_ratio take the next best match.
            %  (repreat till these conditions are met)
            while ~isempty(find((sum(DIF1==0,2)/iter)>repeat_ratio))
                Euclidean_dis2(ooii,ppii)=nan;
                [ooii ppii]=find(Euclidean_dis2==min(Euclidean_dis2(:)));
                NULL_ROIS_temp=NULL_ROIS(:,iter);
                NULL_ROIS_temp(OTHER_elect(ooii))=ppii;
                Index1=find(~isnan(NULL_ROIS_temp));
                DIF1=bsxfun(@minus,NULL_ROIS(Index1,1:iter-1),NULL_ROIS_temp(Index1));
                DIF1(DIF1~=0)=1;
            end
            NULL_ROIS_temp=NULL_ROIS(:,iter);
            NULL_ROIS_temp(OTHER_elect(ooii))=ppii;
            Index1=find(~isnan(NULL_ROIS_temp));
            DIF1=bsxfun(@minus,NULL_ROIS(Index1,1:iter-1),NULL_ROIS_temp(Index1));
            DIF1(DIF1~=0)=1;
            % check if the selected ROIs exist in previous iterations
            if  ~isempty(find(sum(DIF1,1)==0))
                flag_exist =1;
            else
                flag_exist =0;
            end
        end
        clear  Euclidean_dis Euclidean_dis2 NULL_ROIS_temp
        NULL_ROIS(OTHER_elect(ooii),iter)=ppii;
        Null_Dist(OTHER_elect(ooii),selcted_elect_index,iter)=Pair_Dist_roi(ppii,:);
        Null_Dist(selcted_elect_index,OTHER_elect(ooii),iter)=Pair_Dist_roi(ppii,:)';
        selcted_elect_index =[selcted_elect_index;OTHER_elect(ooii)];
        clear ooii ppii  Pair_Dist_roi
        toc
    end
    %% Check if the selcted null ROIS are found in previous iterations
    DIF1= bsxfun(@minus,NULL_ROIS(:,1:iter-1),NULL_ROIS(:,iter));
    DIF1(DIF1~=0)=1;
    if ~isempty(find(sum(DIF1,1)==0)) % flag the ROI pair if true
        DIF2= bsxfun(@minus,FULL_pairs(1:2,:),nn);
        DIF2(DIF2~=0)=1;
        if ~isempty(find(sum(DIF2,1)==0))
            FULL_pairs(3,find(sum(DIF2,1)==0))=FULL_pairs(3,find(sum(DIF2,1)==0))+1;
        else
            FULL_pairs(1:2,qqq)=nn;
            FULL_pairs(3,qqq)=1;
            qqq=qqq+1;
        end
    else
    end
    clear selcted_elect_index
    %% Save Null_ROIs periodically
    if (rem(Iter,10)-(rr))==0
        save([ outpath '/' Pat_Name '_Geometric_SC_Null_' num2str(rr)],'NULL_ROIS','Null_Dist');
        % Find all the unique ROI nulls among all saved nulls.
        tic
        clear  Null_Dist NULL_ROIS
        DIR=dir([ outpath '/' Pat_Name '_Geometric_SC_Null_*.mat']);
        for kjkj = 1: length(DIR)
            kjkj
            load([outpath '/' DIR(kjkj).name]);
            if kjkj ==1
                NULL_ROIS_old=NULL_ROIS(:,~isnan(sum(NULL_ROIS,1)));
                Null_Dist_old=Null_Dist(:,:,~isnan(sum(NULL_ROIS,1)));
            else
                NULL_ROIS=NULL_ROIS(:,~isnan(sum(NULL_ROIS,1)));
                Null_Dist=Null_Dist(:,:,~isnan(sum(NULL_ROIS,1)));
                NULL_ROIS=[ NULL_ROIS_old(:,1:end),NULL_ROIS ];
                Null_Dist=cat(3, Null_Dist_old(:,:,1:end),Null_Dist );
                clear  Null_Dist_old NULL_ROIS_old
               [C,ia,ic]=unique(NULL_ROIS(:,:)','rows');
                NULL_ROIS_old=NULL_ROIS(:,ia);
                Null_Dist_old=Null_Dist(:,:,ia);
                clear  Null_Dist NULL_ROIS
            end
        end
        clear  Null_Dist NULL_ROIS
        [C,ia,ic]=unique(NULL_ROIS_old(:,:)','rows');
        NULL_ROIS=NULL_ROIS_old(:,ia);
        Null_Dist=Null_Dist_old(:,:,ia);
        clear  Null_Dist_old NULL_ROIS_old
        NULL_ROIS2=nan(length(Inddexx),iterations);
        iter= size(NULL_ROIS,2);
        NULL_ROIS2(:,1:iter)=NULL_ROIS;
        clear NULL_ROIS;
        NULL_ROIS=NULL_ROIS2;
        clear NULL_ROIS2;
        if iter >= iterations
            save([ outpath '/' Pat_Name '_Final_Geometric_SC_Null_' num2str(rr)],'NULL_ROIS','Null_Dist');
            break
        end
        toc        
    end
   %%    
end








