#!/bin/bash

conda activate biobombe

python hetnetEnrich_markers.py --collection GpC2CPBIOCARTA;
wait
python hetnetEnrich_markers.py --collection GpC2CPG;
wait
python hetnetEnrich_markers.py --collection GpC2CPKEGG;
wait
python hetnetEnrich_markers.py --collection GpC2CPREACTOME;
wait
python hetnetEnrich_markers.py --collection GpC3TFT;
wait
python hetnetEnrich_markers.py --collection GpC5BP;
wait
python hetnetEnrich_markers.py --collection GpC5CC;
wait
python hetnetEnrich_markers.py --collection GpC5MF;
wait
python hetnetEnrich_markers.py --collection GpC7;
wait
python hetnetEnrich_markers.py --collection GpH;
wait
python hetnetEnrich_markers.py --collection GpMETABASER;
wait
python hetnetEnrich_markers.py --collection GpWIKIPATH;
wait
python hetnetEnrich_markers.py --collection GpXCELL
