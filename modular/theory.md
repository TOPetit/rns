# Operations in Plantard's representation

## Introduction
Plantard's representation is similar to what we can find dealing with Montgomery modular multiplication.
We represent any number $A$ as $AC$ where $C = -2^{2^n}$. Let $\widetilde{A} = A(-2^{2^n})$.

## Modular multiplication
Plantard's Algorithm naturally computes $AB(-2^{2^n}) mod P$, so we keep the same representation.
$$\framebox[1.1\width]{$\widetilde{A}\widetilde{B}\ mod\ P\ =\ \widetilde{AB}$}$$

## Multiplication
Without correcting, we don't have the correct reprentation for the result. We can either shift the result correctly, or use modular multiplication with $(P\ >>\ \widetilde{A}\widetilde{B})$.
$$\framebox[1.1\width]{$\widetilde{A}\widetilde{B}\ =\ AB(-2^{4n})\ =\ \widetilde{AB}(-2^{2n})\ =\ \widetilde{\widetilde{AB}}$}$$

## Addition
Addition in trivial.
$$\framebox[1.1\width]{$\widetilde{A}+\widetilde{B}\ =\ A(-2^{2n})+B(-2^{2n})\ =\ (A+B)(-2^{2n})\ =\ \widetilde{A+B}$}$$

## Soustraction
Same as addition.