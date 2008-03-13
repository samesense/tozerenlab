function [RX_MASK INDS]=BreakoutPatients(PAT_STRUCT,PatientCutoff)
%   BreakoutPatients
%       Returns the different types of therapies within the database and
%       the indexes to each patient which is on that therapy regimine.
%
%
%   [PATIENT_MAT INDS_CELL]=BreakoutPatients(PAT_STRUCT)
%
%   PAT_STRUCT          A patient structure as created by PatStructHelper.
%   
%
%
%   RX_MASK             An N X DRUGS matrix where 0's indicate 'Does not
%                       matter' and 1's indicate that a Drug MUST be
%                       present for inclusion.
%
%   INDS                A cell of indexes for each patient included in the
%                       corresponding row in PATIENT_MAT.
%
%
%
%   See also: PatientStructHelper.
%
%

[RX_data RX_inds]=PatientStructHelper(PAT_STRUCT,{'RX_vals','explodenumeric'});

RX_inds=RX_inds(RX_data(:,1)==0);
RX_data=RX_data(RX_data(:,1)==0,2:end);

NumCols=size(RX_data,2);

RX_MASK=zeros(500,NumCols);
INDS=cell(500,1);
counter=1;

for i=1:NumCols
    pat_mask=RX_data(:,i)==1;
    if nnz(pat_mask)>=PatientCutoff
        INDS{counter}=RX_inds(pat_mask);
        RX_MASK(counter,i)=1;
        counter=counter+1;
    end
    for j=0:i-1
        for k=0:j-1
            pat_mask=mean(RX_data(:,nonzeros([i j k])'),2)==1;
            if nnz(pat_mask)>=PatientCutoff
                INDS{counter}=RX_inds(pat_mask);
                RX_MASK(counter,nonzeros([i j k])')=1;
                counter=counter+1;
            end
        end
    end
end

RX_MASK=RX_MASK(1:counter,:);
INDS=INDS(1:counter);
