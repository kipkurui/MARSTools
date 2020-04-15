# MARSTools

MARSTools contains the commandline version of the tools hosted in the Motif Assessment and Ranking Suite
([MARS](http://www.bioinf.ict.ru.ac.za)). There are three major set of tools:

1. Scoring and Classification-Based MARS:
    * SCORE (*Assess_by_score.py*)
    * *run_gimme.py* (A wrapper from GimmeMotifs gimme roc)
2. Consistency-Based MARS:
    * *run_tomtom.py* -- From the MEME-Suite
    * *run_fisim.py* -- From FISim Tools
3. Enrichment-Based MARS:
    * *run_centrimo.py*

## Setup
To set up [MARSTools](https://github.com/kipkurui/MARSTools), follow the following instructions. Clone the MARSTools repository into a directory you want to setup and create a conda environment using the requirements.yml file.

    ```
    git clone https://github.com/kipkurui/MARSTools

    #create a MARS conda environment
    conda env create -f environment.yml
    #Activate the environment
    source activate MARS
    ```

## TODO:
- Add example data scalled down for testing the tools
- Provide basic scripts for testing the tool
- Fix the errors that existed in the previous tool, as highlighted by the reviewers
- Package the tool for easy install using conda