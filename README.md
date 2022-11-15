# General

We here deposit R and Python scripts to analyze a pile of monitoring datasets 
that monthly recorded incidences of crop pests and diseases for decades
by using two statistical and seven machine learning methods.

### Statistical methods (R)
    1. Bayes model (Stan)
    2. multiple linear regression
### Machine Learning (Python)
    3. random forest
    4. decision tree
    5. kNN
    6. SVM
    7. neural network
    8. elastic-net regression
    9. ridge regression

# Dataset

Datasets we used are depoosited at the site below.
URL: XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
You can freely download all the files.

After you download the zipped dataset, put all the decompressed csv files into
`data/formatted_data/random` and `data/formatted_data/shuffle` directories.
The script `format_data.classic.py` was used for randomizing and shuffling dataset
by myself, but Kishi's dataset (send210725) has already randomized and shffuled,
therefore we do not need run this script any more.

# Files (for machine learning methods)

```bash
cd ${PROJECT_PATH}
cd data

cd formatted_data
mkdir random
mkdir shuffle

# move all csv files in the dataset into random and shuffle directories

mkdir -p random4dnn
for fpath in `ls random`
do
    mkdir -p "random4dnn/${fpath%_random.csv}"
    python ${PROJECT_PATH}/data/format_data.dnn.py "random/${fpath}" random4dnn/${fpath%_random.csv}
done
    
mkdir -p shuffle4dnn
for fpath in `ls shuffle`
do
    mkdir -p "shuffle4dnn/${fpath%_shuffle.csv}"
    python ${PROJECT_PATH}/data/format_data.dnn.py "shuffle/${fpath}" shuffle4dnn/${fpath%_shuffle.csv}
done
```

# Modeling

## Bayes model

State-space Bayes models are estimated by a MCMC algorithm, Stan on R.
Modelling and standardization processes for a single dataset is performed by "bayes_calc.R",
"process_fies_bayes.R" sequentially calls and processes each dataset with "bayes_calc.R".

## Multiple linear regression

Multiple linear regression models for each dataset are estimated by a single R file, "LM_model.R"

## Classic machine learning models

```bash
cd ${PROJECT_PATH}/models
mkdir -p test_results

for alg in svm rf dc knn lasso elasticnet
do
    for ft in category decimal
    do
        printf "\e[31m# ALGORITHM: ${alg}   FEATURE_TYPE: ${ft}\e[m\n"

        python model_classic.py  \
        --algorithm ${alg} \
        --dataset ../data/send210315/cucumber/cucumber__udonkobyo__hompohasseimenseki.csv \
        --feature-type ${ft} \
        --output test_results/test_${alg}_${ft}.tsv \
        --test-run
    done
done


qsub train_model_classic.sh
```

## Deep neural network models

```bash
cd ${PROJECT_PATH}/models
dpath=${PROJECT_PATH}/data/formatted_data
mkdir -p test_results

# check cv training function for category type
python bake.py --mode cv --feature-type category --model L2 \
               --weight test_results/test_L2_category.pth   \
               --cv-dataset ${dpath}/cucumber/randomize4dnn/cucumber__udonkobyo__hompohatsubyoyoritsu \
               --epochs 5 --batch-size 1024

python bake.py --mode cv --feature-type decimal --model L2 \
               --weight test_results/test_L2_decimal.pth   \
               --cv-dataset ${dpath}/cucumber/randomize4dnn/cucumber__udonkobyo__hompohatsubyoyoritsu \
               --epochs 5 --batch-size 1024


qsub train_model_dnn.sh
```


## Validation

Summarise validation results of DNN and classic models.

```
cd ${PROJECT_PATH}/models

python summarise_valid.py ${PROJECT_PATH}/data/cv_results


cd ${PROJECT_PATH}/data
R
> source('summarise_valid.R')
```


