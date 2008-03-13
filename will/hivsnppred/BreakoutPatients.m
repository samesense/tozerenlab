function [RX_MASK INDS]=BreakoutPatients(PAT_STRUCT,LAST_THERAPY,PatientCutoff)
%   BreakoutPatients
%       Returns the different types of therapies within the database and
%       the indexes to each patient which is on that therapy regimine.
%       Functions as an Iterator.  The function will return the next
%       therapy regimine which has more than a thresh-holded number of
%       patients.
%
%
%   [PATIENT_MAT INDS_CELL]=BreakoutPatients(PAT_STRUCT)
%
%   PAT_STRUCT          A patient structure as created by PatStructHelper.
%   
%   LAST_THERAPY        The Vector of the last therapy processed.
%
%
%
%   RX_MASK             An N X DRUGS matrix where NaNs indicate 'Does not
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


start_I=1;
start_J=0;
start_K=0;

if ~isempty(LAST_THERAPY)
    start_I=find(~isnan(LAST_THERAPY),1,'last');
    if any(~isnan(LAST_THERAPY(1:start_I-1)))
        start_J=find(~isnan(LAST_THERAPY(1:start_I-1)),1,'last');
        if any(~isnan(LAST_THERAPY(1:start_J-1)))
            start_K=find(~isnan(LAST_THERAPY(1:start_J-1)),1,'last');
        end
    end
end


[RX_data RX_inds]=PatientStructHelper(PAT_STRUCT,{'RX_vals','explodenumeric'});

RX_inds=RX_inds(RX_data(:,1)==0);
RX_data=RX_data(RX_data(:,1)==0,2:end);

NumCols=size(RX_data,2);


RX_MASK=NaN(1,NumCols);

for i=start_I:NumCols
    pat_mask=RX_data(:,i)==1;
    if nnz(pat_mask)>=PatientCutoff
        INDS=RX_inds(pat_mask);
        RX_MASK(i)=1;
        if ~all(isequalwithequalnans(RX_MASK,LAST_THERAPY))
            return
        else
            RX_MASK=NaN(1,NumCols);
        end
        
    end
    for j=start_J:i-1
        for k=start_K:j-1
            pat_mask=mean(RX_data(:,nonzeros([i j k])'),2)==1;
            if nnz(pat_mask)>=PatientCutoff
                INDS=RX_inds(pat_mask);
                RX_MASK(nonzeros([i j k])')=1;
                if ~all(isequalwithequalnans(RX_MASK,LAST_THERAPY))
                    return
                else
                    RX_MASK=NaN(1,NumCols);

                end
            end
        end
        start_K=0;

    end
    start_J=0;
end


    

    
    
    
















