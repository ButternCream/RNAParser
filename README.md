# RNA Parser

**Parses files which contain an RNA sequence and outputs the found peptides to a file**

**For help**
```
python3 app.py -h
```

## Examples 
**Using files in examples folder**

To simply convert an RNA sequence in a file
```
python3 app.py -f examples/HomoSapien.RNA
```

Specify the output directory with --output or -o, default is 'Peptides'
```
python3 app.py --output=Converted -f examples/HomoSapien.RNA
```

Specify the file extension of the output with --extension or -e
```
python3 app.py -e=out -f examples/HomoSapien.RNA
```
will create HomoSapien.out in Peptides folder

Delete a specified directory with -c or --clean
```
python3 app.py -c Peptides
```

Parse multiple files
```
python3 app.py -f examples/HomoSapien.RNA examples/Drosophila_Melanogaster.RNA examples/Sedum_Takesimense.RNA
```

Clean multiple folders
```
python3 app.py -c Peptides Folder2 Folder3
```
