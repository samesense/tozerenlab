function RunMultiAnal(PAT_STRUCT,EXCEL_FILENAME)
%   RunMultiAnal
%       Runs the PredictPatients function on all therapy regimines and
%       records the output into an Excel file.  If the Excel file already
%       exists then it continues from the last iteration recorded.
%
%
%   RunMultiAnal(PAT_STRUCT,EXCEL_FILENAME)   
%
%   PAT_STRUCT      A patient structure as created by PatientStructHelper
%
%   EXCEL_FILENAME  An excel file to save the data to.
%
%
%


RX_names=PAT_STRUCT{1}.RX_names';


try
    
    [Numeric_2 Text_2 Raw_2]=xlsread(EXCEL_FILENAME,2);
    
    LastRegimine=Numeric_2(end,:);
    LastIND=size(RAW_2,1);
    
catch
    
    xlswrite(EXCEL_FILENAME,RX_names,1);
    xlswrite(EXCEL_FILENAME,RX_names,2);

    LastRegimine=[];
    LastIND=1;

end



while(1)

    [NewRegimine RegimineInds]=BreakoutPatients(PAT_STRUCT,LastRegimine,200);

    NewPreds=PredictPatients(PAT_STRUCT(RegimineInds),500);

    Zero_regimine=NewRegimine;
    Zero_regimine(isnan(Zero_regimine))=0;
    
    xlswrite(EXCEL_FILENAME,Zero_regimine,2,['A' num2str(LastIND+1)]);
    xlswrite(EXCEL_FILENAME,NewPreds,1,['A' num2str(LastIND+1)]);

    LastIND=LastIND+1;
    LastRegimine=NewRegimine;
    
    fprintf('Finished Running Index %d',LastIND);


end


