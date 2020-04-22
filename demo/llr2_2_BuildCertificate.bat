move /Y lresults.txt lresults.txt.bak
llr2 -d -oFFT_Increment=1 -oDeletePoints=0 -oProofName=proof -oProofCount=16 -oRndSeed=12345678901234567890 -pBuildCert -q"1201*2^1573648+1"
move /Y lresults.txt cert.txt
move /Y lresults.txt.bak lresults.txt
