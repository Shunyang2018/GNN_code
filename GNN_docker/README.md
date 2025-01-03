# Retention Time Prediction on a CSV Dataset


## Docker run 

A docker image (gnn:1.0) is built on the wormwood server.
To predict RT you just need to change the name and -tp -sc options. 

More detailed info about options are [listed](#a) 

```bash
docker run -v $PWD:/data shunyang2018/gnn -f /data/<name>.csv -sc smiles -tp <c18> -ip data/
```
The output is prediction.csv at the same folder. And then, you can use RT-filter jupyter notebook to process the data. 

[Batch mode](#b) is another option in the next section. 

<a name="b"></a>
## Batch run on the wormwood



We add 
- **Batch mode**: `-batch True` [default=False]
    - whether use batch run to add RT fileter column to the original input csv file, 
    - the output _rt.csv file will be at the same path of input csv.
Run one file in batch mode:      
```bash
docker run -v $PWD:/data shunyang2018/gnn -f /data/<name>.csv -sc smiles -tp <c18> -ip data/ -batch=True
```
Run every files under the same folder: 

run.sh
```bash
#!/bin/sh
mkdir C18
mv *C18*.csv C18/

cd C18/
# loop run and output rt.csv file
for f in *.csv; do echo "running ${f}"; docker run -v $PWD:/data shunyang2018/gnn -f "/data/${f}" -sc query_smiles -tp c18 -ip data/ -batch=True;done
```

    
## Build Docker

```bash 
docker build -f Dockerfile2 -t gnn:1.0 .
```

<a name="a"></a>
## Prediction options

To use the pre-trained model trained above for prediction on new molecules

```bash
python regression_inference.py -f X
```

where `X` specifies the path to a `.csv`/`.txt` file of SMILES strings

Other optional arguments include:
- **SMILES Column**: `-sc column` [default=None]
    - Specifies the column of SMILES strings in the input `.csv` file. Can be omitted if the input file is a
      `.txt` file or the `.csv` file only has one column of SMILES strings
- **Column Type**: `-tp pfp/c18` [default=c18]
    - Path to the training results saved, which will be used for loading the trained model and related configurations
- **Inference Result Path**: `-ip path` [default=regression_inference_results]
    - Path to the inference results, which will be used to save:
        - `prediction.csv`: A file of predicted properties associated with SMILES strings
- **Task**: `-t task1,task2,task3,...`
    - Task names for saving model predictions in the CSV file to output, which should be the same as the ones
    used for training. If not specified, we will simply use task1, task2, ...
- **Num of Processes for Data Loading**: `-nw number` [default=1]
    - Specifies the number of processes to use for data loading. A larger number might yield a faster speed.

    
