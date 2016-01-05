# POSTME
POSTME : POSTprocessing MEtagenomics assembly with hybrid data.

To use the tool, use the following command(asuming illuminaContigs.fasta and pacbioReads.fasta are in workingDirectory). The output file is improved.fasta

    python postme.py \
        -o workingDirectory/ \
        -a mummerPath/ \
        -c illuminaContigs.fasta \
        -r pacbioReads.fasta

Unit tests are continuously tested on Travis CI and the current build status is  
![alt text](https://travis-ci.org/kakitone/postme.svg?branch=master "Current build status")

## Example
1. Clone POSTME

        git clone https://github.com/kakitone/postme.git

2. Install dependency, you can do it through virtualenv
    
        pip install virtualenv
        source venv/bin/activate
        pip install numpy scipy biopython matplotlib
    
3. To make sure you have installed all the dependency,

        python postme/myunittest.py

4. You can then run POSTME then.
