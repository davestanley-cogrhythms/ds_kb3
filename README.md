ds_kb3
=====

New version of LFP/unit analysis code, starting after submission of SFN 2016 poster & Cerebral Cortex paper getting accepted.

### Getting started (Mac / Linux)
=====

### Clone main repo
  
    git clone git@github.com:davestanley-cogrhythms/ds_kb3.git
    
OR
    
    git clone --recursive https://github.com/davestanley-cogrhythms/ds_kb3.git
  
### Clone necessary repos into src folder
  
    mkdir ~/src
    cd ~/src
    git clone git@github.com:cogrhythms/chronux
    git clone git@github.com:cogrhythms/spike_field_assoc_dev
    git clone git@github.com:davestanley-cogrhythms/lib_dav
    git clone git@github.com:davestanley/SigProc-Plott.git
    cd SigProc-Plott
    git checkout dev
    git pull
  
Setup instructions are the same for Windows, but using Windows syntax

