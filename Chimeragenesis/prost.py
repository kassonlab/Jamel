from transformers import T5Tokenizer, T5EncoderModel
import torch

from Chimeragenesis.AccessiontoAlignment import create_dictionary_from_alignment


def prost_embed_aln(aln_file, new_embed_file):
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    tokenizer = T5Tokenizer.from_pretrained('Rostlab/ProstT5', do_lower_case=False)
    model = T5EncoderModel.from_pretrained("Rostlab/ProstT5").to(device)

    # only GPUs support half-precision currently; if you want to run on CPU use full-precision (not recommended, much slower)
    model.full() if device=='cpu' else model.half()

    protein_sequences = create_dictionary_from_alignment(aln_file)
    protein_sequences = {label:f"<AA2fold> {' '.join(list(seq.replace('-','X')))}" for label,seq in protein_sequences.items()}
    ids = tokenizer.batch_encode_plus(protein_sequences.values(),
                                      add_special_tokens=True,
                                      padding="longest",
                                      return_tensors='pt').to(device)

    with torch.no_grad():
        embedding_rpr = model(
                  ids.input_ids,
                  attention_mask=ids.attention_mask
                  )
    embed_dict = {id: tensor for id, tensor in zip(protein_sequences.keys(), torch.split(embedding_rpr.last_hidden_state,1,dim=0))}
    torch.save(embed_dict, new_embed_file)
# tensors=torch.load('save_embed').last_hidden_state
prost_embed_aln(r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\esm\schema_msa.aln','prost.pkl')
