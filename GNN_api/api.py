import json
import os
import pandas as pd
import torch

from dgllife.data import UnlabeledSMILES
from dgllife.utils import MolToBigraph
from torch.utils.data import DataLoader
from tqdm import tqdm

from utils import mkdir_p, collate_molgraphs_unlabeled, load_model, predict, init_featurizer, combine

from rdkit import Chem

from flask import Flask, request
from collections import defaultdict

def smi2rt(smiles_string, column_type):
    args = defaultdict(list)
    args['smiles'] = [smiles_string]
    args['column_type'] = column_type.upper()
    with open(column_type + '/configure.json', 'r') as f:
        args.update(json.load(f))
    if torch.cuda.is_available():
        args['device'] = torch.device('cuda:0')
    else:
        args['device'] = torch.device('cpu')
    args = init_featurizer(args)
    mol_to_g = MolToBigraph(add_self_loop=True,
                            node_featurizer=args['node_featurizer'],
                            edge_featurizer=args['edge_featurizer'])
    dataset = UnlabeledSMILES(args['smiles'], mol_to_graph=mol_to_g)
    dataloader = DataLoader(dataset, batch_size=args['batch_size'],
                            collate_fn=collate_molgraphs_unlabeled, num_workers=1)
    model = load_model(args).to(args['device'])
    checkpoint = torch.load(args['column_type'] + '/model.pth', map_location='cpu')
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()

    smiles_list = []
    predictions = []

    with torch.no_grad():
        for batch_id, batch_data in enumerate(dataloader):
            batch_smiles, bg = batch_data
            smiles_list.extend(batch_smiles)
            batch_pred = predict(args, model, bg)
            predictions.append(batch_pred.detach().cpu())

    predictions = torch.cat(predictions, dim=0)
    inchikey_list = []
    for smi in args['smiles']:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        inchikey_list.append(Chem.MolToInchiKey(mol))

    output_data = {'canonical_smiles': smiles_list,'inchikey':inchikey_list}

    output_data['RT'] = round(predictions.item(),4)
    return output_data


server = Flask(__name__)


@server.route("/gnn")
def classify():
    smiles_string = request.values.get("smiles")
    column_type = request.values.get("column")
    respond_dict = smi2rt(smiles_string, column_type)

    return json.dumps(respond_dict)



if __name__ == "__main__":
	# for debugging locally

	# for production
	server.run(host='0.0.0.0', port=9002)
