The scripts and files in here are used to generate data files that are then used by the test suite in Biff. To prepare the test data files, do:

```
make
make run
```

Then you need to move them in to the data path. From the top level Biff directory, do:

```
cp test-data-helper/*.dat.gz biff/data/
cp test-data-helper/*.coeff biff/data/
```
