# General

We here deposit R and Python scripts to analyze a pile of crop pest and disease monitoring data
by using two statistical and seven machine learning methods.

### Statistical methods:
    Bayes model and multiple linear regression model
### Machine Learning:
    random forest, decision tree, kNN, SVM, neural network, elastic-net regression and ridge regression

# Dataset

Monitoring datasets we used are depoosited at the site below.
URL: XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
You can freely download all the files.

After you download the zipped dataset, put all the decompressed csv files into
`data/formatted_data/random` and `data/formatted_data/shuffle` directories.
The script `format_data.classic.py` was used for randomizing and shuffling dataset
by myself, but Kishi's dataset (send210725) has already randomized and shffuled,
therefore we do not need run this script any more.

# Files

```bash
cd ${PROJECT_PATH}
cd data

cd formatted_data
mkdir random
mkdir shuffle

# move all csv files in Kishi's dataset into random and shuffle directories

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


