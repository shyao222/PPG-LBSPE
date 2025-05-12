from rdkit import Chem
from copy import deepcopy
import numpy as np
from private.hypergraph import Hypergraph, hg_to_mol
import pickle
import pandas as pd

def random_produce_2K(grammar):
    def sample(l, prob=None):
        if prob is None:
            prob = [1/len(l)] * len(l)
        idx =  np.random.choice(range(len(l)), 1, p=prob)[0]
        return l[idx], idx

    def prob_schedule(_iter, selected_idx):
        prob_list = []
        # prob = exp(a * t * x), x = {0, 1}
        a = 0.5
        for rule_i, rule in enumerate(grammar.prod_rule_list):
            x = rule.is_ending
            if rule.is_start_rule:
                prob_list.append(0)
            else:
                prob_list.append(np.exp(a * _iter * x))
        prob_list = np.array(prob_list)[selected_idx]
        prob_list = prob_list / np.sum(prob_list)
        return prob_list

    hypergraph = Hypergraph()
    starting_rules = [(rule_i, rule) for rule_i, rule in enumerate(grammar.prod_rule_list) if rule.is_start_rule]
    iter = 0
    k_count = 0  # 初始化钾原子计数器

    while(True):
        if iter == 0:
            _, idx = sample(starting_rules)
            selected_rule_idx, selected_rule = starting_rules[idx]
            hg_cand, _, avail = selected_rule.graph_rule_applied_to(hypergraph)
            hypergraph = deepcopy(hg_cand)
        else:
            candidate_rule = []
            candidate_rule_idx = []
            candidate_hg = []
            for rule_i, rule in enumerate(grammar.prod_rule_list):
                hg_prev = deepcopy(hypergraph)
                hg_cand, _, avail = rule.graph_rule_applied_to(hypergraph)
                if(avail):
                    candidate_rule.append(rule)
                    candidate_rule_idx.append(rule_i)
                    candidate_hg.append(hg_cand)
            if (all([rl.is_start_rule for rl in candidate_rule]) and iter > 0) or iter > 30:
                break
            prob_list = prob_schedule(iter, candidate_rule_idx)
            hypergraph, idx = sample(candidate_hg, prob_list)
            selected_rule = candidate_rule_idx[idx]

        # 更新钾原子计数器
        k_count = sum(1 for node in hypergraph.nodes if node == 'K')

        # 如果钾原子数量达到2个，则停止生成
        if k_count == 2:
            break
        iter += 1

    # 生成分子并验证钾原子数量
    try:
        mol = hg_to_mol(hypergraph)
        smiles = Chem.MolToSmiles(mol)
        if smiles.count('[K]') == 2:
            print(smiles)
            return mol, iter
        else:
            return None, iter
    except:
        return None, iter
    

def generate(num):
    
    # Limit the requested number of molecules to a maximum of 10;
    # simply comment out or remove this line to disable the cap.
    num = min(num, 10)

    generated = []
    with open('utilities/models/epoch_grammar_44_1.5142501536779158.pkl', 'rb') as fr:
        grammar = pickle.load(fr)

    max_attempts = num * 10  # Prevent infinite loop; adjust multiplier as needed
    attempts = 0

    while len(generated) < num and attempts < max_attempts:
        mol, _ = random_produce_2K(grammar)
        if mol is not None:
            generated.append(mol)
        attempts += 1

    generated_smiles = [Chem.MolToSmiles(mol) for mol in generated]
    df_generated_smiles = pd.DataFrame(generated_smiles, columns=['SMILES'])
    return df_generated_smiles