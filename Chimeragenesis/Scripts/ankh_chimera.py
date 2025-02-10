import ankh
import torch
import transformers
from AccessiontoAlignment import create_dictionary_from_alignment

def ankh_embed_aln(aln_file, new_embed_file):
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    model, tokenizer = ankh.load_large_model()
    model.eval()
    protein_sequences = create_dictionary_from_alignment(aln_file)
    outputs = tokenizer.batch_encode_plus([list(seq) for seq in protein_sequences.values()],
                                          add_special_tokens=False,
                                          padding=False,
                                          is_split_into_words=True,
                                          return_tensors="pt").to(device)
    with torch.no_grad():
        embeddings: transformers.modeling_outputs.BaseModelOutputWithPastAndCrossAttentions = model(
            input_ids=outputs['input_ids'], attention_mask=outputs['attention_mask'])
        embed_dict = {id: tensor for id, tensor in zip(protein_sequences.keys(), embeddings.last_hidden_state)}
    torch.save(embed_dict, new_embed_file)


if __name__ == '__main__':
    ankh_embed_aln(r"/sfs/weka/scratch/jws6pq/Notebook/ESM/Schema_valdation/labeled_schema_aln","/sfs/weka/scratch/jws6pq/Notebook/ESM/Schema_valdation/large_ankh.pkl")

