#! /usr/bin/python

import sys
import random
import torch
import esm

def perSeqPerplexity(model,input_tensor,mask_tok_idx,layer):
    repeats = input_tensor.repeat(input_tensor.size(-1)-2,1)
    mask = torch.ones(input_tensor.size(-1) -1).diag(1)[:-2]
    masked_input = repeats.masked_fill(mask == 1, mask_tok_idx)
    with torch.no_grad():
        results = model(masked_input, repr_layers=[layer], return_contacts=True)

    logits = results['logits'][:,1:-1,4:24]
    #softmax = torch.nn.Softmax(dim=2)
    #probs = torch.diagonal(softmax(logits))
    loss = torch.nn.CrossEntropyLoss(ignore_index=1,reduction='mean')
    masked_logits = torch.transpose(torch.diagonal(logits),0,1)
    target = input_tensor[1:-1]-4
    return float(loss(masked_logits,target))


def calPerplexity(data):
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()

    batch_converter = alphabet.get_batch_converter()
    input_labels, input_strs, input_tokens = batch_converter(data)
    ppl = {}
    for i,input_tensor in enumerate(input_tokens):
        print('Calculate perplexity for %s' % input_labels[i])
        score = perSeqPerplexity(model,input_tensor,alphabet.tok_to_idx['<mask>'],33)
        ppl.update({input_labels[i]:score})

    print(ppl)
    

def main():
    data = [("protein2", "KALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEE")]
    for i in range(10):
        copied_seq = list(data[0][1])
        random.shuffle(copied_seq)
        data.append(('protein2_shuffle%d' % (i+1), ''.join(copied_seq)))
    print(data)

    calPerplexity(data)


if __name__=='__main__':
    main()
