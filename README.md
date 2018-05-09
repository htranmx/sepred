# Side Effect Prediction

This code is part of an exploratory project in medicinal chemistry for predicting drug side effects.

We explore the first algorithm (canonical correlation analysis) presented in the paper:

Reference: Atias, N, Sharan, R (2011). An algorithmic framework for predicting side effects of drugs. J. Comput. Biol., 18, 3:207-18.

Files:

* README.md                        : this file
* sepred.py                        : driver, main program
* sepred_utils.py                  : implementation of algorithm, heavy lifting is here
* sepred_reader_effects.py         : parse and load input data files
* sepred_reader_binary_effects.py  : parse and load input data files, converting effect values into binary values.
                                   need to merge into sepred_reader_effects.py or remove altogether
* sepred_tester.py                 : testing accuracy of the model

Main (driver) program is sepred.py. Flags for input:
* -a : file of descriptors for each compound in dataset
* -e : file of GI50 values for compounds
* -k : dimension of projected subspace (see reference).

Example:

First, create the model:

```bash
sepred.py -a nci60_set1.csv -e GI50_RAW.txt -k 5
```

This generates a model and saves it into the file "ccmodel.pkl"

Use model to make predictions:

```bash
sepred_tester.py -m ccmodel.pkl -a nci60_set2.csv
```

Output contains (actual, predicted) pairs.
