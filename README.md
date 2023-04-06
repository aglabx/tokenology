# tokenology
DNA and tokens


# Usage

For GENA_LM tokens:

```python
from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained('aglabx/dna_tokens', force_download=True, use_fast=True)

print(tokenizer.vocab_size)

tokenizer.tokenize(dna_data.upper())
```

For 16S tokens:

```python
from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained('aglabx/16S_1024_bpe_tokens', force_download=True, use_fast=True)

print(tokenizer.vocab_size)

tokenizer.tokenize(dna_data.upper())
```
