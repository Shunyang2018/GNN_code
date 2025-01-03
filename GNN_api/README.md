# Retention Time API server on wormwood


## API run example

http://<server_ip>:<port>/gnn?smiles=<**SMILES**>&column=<**C18/PFP**>


```bash
http://<server_ip>:<port>/gnn?smiles=CC(=O)OC1=CC=CC=C1C(=O)O&column=C18

output: {"canonical_smiles": ["CC(=O)Oc1ccccc1C(=O)O"], "inchikey": ["BSYNRYMUTXBXSQ-UHFFFAOYSA-N"], "RT": 3.609}

```

    
## Build docker

```
docker build -f Dockerfile2 -t gnn:1.0 .
docker run  -d -p 9002:9002 gnn_api:1.0
```
