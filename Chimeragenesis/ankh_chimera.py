import ankh
import torch
import transformers
import pandas as pd
from AccessiontoAlignment import create_dictionary_from_alignment

model, tokenizer = ankh.load_base_model()
model.eval()
schema_data = pd.read_csv(r"C:\Users\jamel\Downloads\schema_data.csv")
schema_data['parent_label'] = schema_data['chimera_block_ID'].apply(
    lambda parent_id: tuple(sorted(set(str(parent_id)))))
schema_data['ankh_score']=schema_data['parent_label']
protein_sequences = create_dictionary_from_alignment(r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\esm\schema_msa.aln')
parent_scores = {'0': 'c0000000000', '1': 'c1111111111', '2': 'c2222222222'}

outputs = tokenizer.batch_encode_plus(list(protein_sequences.values()),
                                      add_special_tokens=False,
                                      padding=True,
                                      is_split_into_words=False,
                                      return_tensors="pt")

with torch.no_grad():
    embeddings:transformers.modeling_outputs.BaseModelOutputWithPastAndCrossAttentions = model(input_ids=outputs['input_ids'], attention_mask=outputs['attention_mask'])
    embed_dict={id:tensor for id,tensor in zip(protein_sequences.keys(),embeddings.last_hidden_state)}

    for row,data in schema_data.iterrows():
        if len(data['parent_label'])==3:
            parent1,parent2,_=data['parent_label']
            chimera_tensor=embed_dict[data['chimera_block_ID']]
            schema_data.loc[row,'ankh_score']=(abs(torch.sum(chimera_tensor - embed_dict[parent_scores[parent1]])) + abs(torch.sum(chimera_tensor - embed_dict[parent_scores[parent2]]))).item()

schema_data.to_csv('test.csv')