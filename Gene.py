from typing import List

class Gene:
    def __init__(self: object,
                 seq: str,
                 gene_type: str, 
                 epitopes: List[str]) -> None:
        self.seq = seq
        self.gene_type = gene_type
        self.epitopes = epitopes

    def __str__(self: object) -> str:
        return f"Gene(seq={self.seq}, gene_type={self.gene_type}, epitopes={self.epitopes})"

    def __json__(self: object) -> dict:
        return {
            'seq': self.seq,
            'gene_type': self.gene_type,
            'epitopes': self.epitopes
        }