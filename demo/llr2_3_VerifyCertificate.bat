move /Y lresults.txt lresults.txt.bak
llr2 -d -oProofName=proof -pVerifyCert -q"1201*2^1573648+1"
move /Y lresults.txt cert_dc.txt
move /Y lresults.txt.bak lresults.txt
fc cert.txt cert_dc.txt
pause
