"""Perl value semantics used by the ports of the original Perl scripts (makebox, makespheres).

The scripts split lines on whitespace and use the tokens as numbers or strings. These
helpers reproduce how Perl treats such tokens, so the ports behave identically even on
malformed input (e.g. PDB columns that run together).
"""
import re

_LEADING_NUMBER = re.compile(r"\s*[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def field(tokens, i):
    """`$tokens[i]`: a missing token is undef, which Perl treats as ""."""
    return tokens[i] if i < len(tokens) else ""


def num(token):
    """Numeric value of a string as Perl computes it: its leading number, else 0."""
    match = _LEADING_NUMBER.match(token)
    return float(match.group()) if match else 0.0


def read_tokens(path):
    """Whitespace-split tokens of every line of a file (Perl: `push @items, [split]`)."""
    with open(path) as f:
        return [line.split() for line in f]
