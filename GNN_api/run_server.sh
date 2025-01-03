#!/bin/bash


#python ./api.py
gunicorn --worker-class=gthread -b 0.0.0.0:9002 --timeout 20 api:server --access-logfile ./access.log --max-requests 250

# http://0.0.0.0:9002/gnn?smiles=CC(=O)OC1=CC=CC=C1C(=O)O&column=PFP
