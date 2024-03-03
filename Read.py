from typing import List

class Read:
    def __init__(self: object,
                 seq: str, 
                 v_gene: str,
                 d_gene: str, 
                 j_gene: str, 
                 epitopes: List[str]) -> None:
        self.seq = seq
        self.v_gene = v_gene
        self.d_gene = d_gene
        self.j_gene = j_gene
        self.epitopes = epitopes

    def __str__(self: object) -> str:
        return f"Read(seq={self.seq}, v_gene={self.v_gene}, d_gene={self.d_gene}, j_gene={self.j_gene}, epitopes={self.epitopes})"
    
    def __json__(self: object) -> dict:
        return {
            'seq': self.seq,
            'v_gene': self.v_gene,
            'd_gene': self.d_gene,
            'j_gene': self.j_gene,
            'epitopes': self.epitopes
        }
