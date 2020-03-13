
      Atnashev Proof demo


  llr2_0_SimpleTest.bat is just a regular test. It uses the following temporary files:

    zXXX - checkpoint; the latest point in test.
    rXXX - recoverypoint; if Gerbicz check fails, the test rolls back to this point.
    lresults.txt for the result.

  llr2_1_SavePoints.bat is what the first tester runs. It produces the residue, as well
as the following files:

    rXXX.0 - contains a^k and is the only thing double-checked by the authority.
    rXXX.YY - 20-40 recoverypoints saved at regular intervals
    rXXX.s - the last recoverypoint. Can be used to quickly produce the residue. But
needs to be proven correct.
    lresults.txt for the result.

  llr2_2_BuildCertificate.bat builds the certificate and its residue, computes the
residue of the test.

    cXXX - "the challenge".
    cert.txt - the expected response.

  llr2_3_VerifyCertificate.bat is THE double check. Calculates the response.

    cert_dc.txt - the double-checked response. Should match the expected response
unless cheating is involved.
