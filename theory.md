# Pseudo-Mersenne modulo:

\begin{math}
Let\ a = \sum_{i=0}^{127}{a_i2^i}\ and\ m = 2^{63}-k \newline
\exists!\ (q,r)\ where\ 0 \leq r < m\ and\ a = qm+r \newline
a = q(2^{63}-k) + r \newline
a = 2^{63}q - kq + r \newline
r = kq - 2^{63}q - a\ with\ q = \sum_{i=63}^{127}{a_i2^{i-63}} \newline
r = kq - \sum_{i=63}^{127}{a_i2^i} - a \newline
r = k\sum_{i=63}^{127}{a_i2^{i-63}} - \sum_{i=63}^{127}{a_i2^i} - \sum_{i=0}^{127}{a_i2^i} \newline
r = k\sum_{i=63}^{127}{a_i2^{i-63}} - ( - \sum_{i=0}^{62}{a_i2^i}) \newline
r = k\underbrace{\sum_{i=63}^{127}{a_i2^{i-63}}}_{upper\ part} + \underbrace{\sum_{i=0}^{62}{a_i2^i}}_{lower\ part} \newline
\end{math}