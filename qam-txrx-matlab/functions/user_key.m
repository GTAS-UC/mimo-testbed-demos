function key=user_key(userno)

stg=7;

taps_table={
    {[2 5],[1 2 4 5],[2 3 4 5]}; %5 stages XCORR not tested
    {[1 6],[1 2 5 6],[2 3 5 6]}; %6 stages XCORR not tested
    {[1 7], [3 7], [2 4 6 7], [4 5 6 7]} %7 stages GOOD XCORR
    %[1 2 3 7]; Martin 
    };

taps = taps_table{stg-4}{userno};
inidata = ones(1,stg);
key = mseq(stg, taps, inidata);