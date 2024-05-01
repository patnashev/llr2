DEPRECATED: The source code is no longer maintained and no further releases are planned. Replaced by [PRST](https://github.com/patnashev/prst) utility.

LLR2 is a primality testing program for numbers of several specific forms. It uses Gwnum library by George Woltman, original LLR code by Jean Penné, hardware error check by Robert Gerbicz and test verification scheme by Pavel Atnashev with verifiable delay function by Krzysztof Pietrzak.

Main differences from the original LLR:
- Gerbicz check in Proth test.
- Gerbicz check in Fermat PRP test for b=2.
- Gerbicz check in Fermat PRP test for b!=2 using Sliding Window exponentiation method (limits performance penalty to 15-20%).
- Native APRCL test (no need to link additional libraries or have satellite programs).
- [Proth/PRP test verification scheme](https://www.mersenneforum.org/showthread.php?t=25323) with certification and certificate validation (fast double-check).
- [Pietrzak VDF](https://eprint.iacr.org/2018/627.pdf) for more efficient test verification (saves on bandwidth).
- Updated checkpoint file format with test fingerprint.
- Checkpoint write through file cache (on Windows).
- Minor changes and fixes.
