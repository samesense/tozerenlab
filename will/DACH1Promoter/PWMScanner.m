function RESULT=PWMScanner(INPUT_PWM,INPUT_SEQ,varargin)
%   PWMScanner
%       Scans through input sequence with the given Position Weight Matrix.
%        The resulting vector is the score of the PWM when it is placed at
%        RESULT(i:i+length(PWM)), trailing zeros are left at the end of the
%        sequence where i+length(PWM)>length(INPUT_SEQ)
%
%
%   RESULT=PWMScanner(INPUT_PWM,INPUT_SEQ)
%
%   INPUT_PWM       A matrix of the frequency of each BASE found in the
%                   motif of interest.
%                       Must be an NUM_BASES X MOTIF_LENGTH matrix of
%                       values between 0 and 1.  All columns must sum to 1.
%
%   INPUT_SEQ       The sequences for scanning.
%                       A CHAR array.
%                       A cell of strings.
%                       A structure with a Sequence feild.
%
%
%   RESULT         A vector of the PWM score at each position normalized by
%                  the maximum possible value of INPUT_PWM.  If INPUT_SEQ
%                  holds multiple sequences then RESULT is a cell where
%                  each cell contains the RESULT of the corresponding
%                  sequence.
%
%
