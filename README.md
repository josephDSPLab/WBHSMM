# WBHSMM
## Wide Band Harmonic Sinusoidal Modeling in Matlab

A matlab implementation for wide band hamonics harmonic sinusoidal modeling.
The original paper is here
http://mtg.upf.edu/node/1020

Current state:  
Analysis phase and sysnthesis phase was constructed. However, it still work not well(probably caused by different pulse onset detection tech).

### TODO
Testing the **maximally flat phase alignment(MFPA)** as pulse onset detection tech.

### Data/Directory info

`/pitch detection`: **YAAPT** algorithm for pitch detection

`OnsetDetection.m`: Use maximum/minimum local difference criterion. It is necessary to modify the algorithm first.
    
`WBHSM_DEMO.m`:Main script

`WBHSM_ana.m`:Analysis part with periodization method

`WBHSM_syn`:Synthesis part 