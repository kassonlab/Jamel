import ankh
import torch
import transformers
from AccessiontoAlignment import create_dictionary_from_alignment

def ankh_embed_aln(aln_file, new_embed_file):
    model, tokenizer = ankh.load_large_model()
    model.eval()
    protein_sequences = create_dictionary_from_alignment(aln_file)
    outputs = tokenizer.batch_encode_plus([list(seq) for seq in protein_sequences.values()],
                                          add_special_tokens=False,
                                          padding=False,
                                          is_split_into_words=True,
                                          return_tensors="pt")
    with torch.no_grad():
        embeddings: transformers.modeling_outputs.BaseModelOutputWithPastAndCrossAttentions = model(
            input_ids=outputs['input_ids'], attention_mask=outputs['attention_mask'])
        embed_dict = {id: tensor for id, tensor in zip(protein_sequences.keys(), embeddings.last_hidden_state)}
        torch.save(embed_dict, new_embed_file)


def load_embed_file(embed_file):
    return torch.load(embed_file)


if __name__ == '__main__':
    ankh_embed_aln(r"/sfs/weka/scratch/jws6pq/Notebook/ESM/Schema_valdation/labeled_schema_aln","/sfs/weka/scratch/jws6pq/Notebook/ESM/Schema_valdation/large_ankh.pkl")
    # schema_data = pd.read_csv(r"C:\Users\jamel\Downloads\schema_data.csv")
    # schema_data['parent_label'] = schema_data['chimera_block_ID'].apply(
    #     lambda parent_id: tuple(sorted(set(str(parent_id)))))
    # schema_data['ankh_score'] = schema_data['parent_label']
    # parent_scores = {'0': 'c0000000000', '1': 'c1111111111', '2': 'c2222222222'}
    # embed_dict=load_embed_file('large_model_amkh.pkl')
    #
    # for row,data in schema_data.iterrows():
    #     if len(data['parent_label'])==3:
    #         parent1,parent2,_=data['parent_label']
    #         parent2_tensor=embed_dict[parent_scores[parent2]]
    #         parent1_tensor = embed_dict[parent_scores[parent1]]
    #         chimera_tensor=embed_dict[data['chimera_block_ID']]
    #         schema_data.loc[row,'ankh_score']=(ESM.tensor_distance(parent1_tensor,chimera_tensor) +ESM.tensor_distance(parent2_tensor,chimera_tensor)).item()
    # schema_data.to_csv('large_ankh.csv')
