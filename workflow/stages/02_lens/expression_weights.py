"""Explicit event-expression proxies; no calibrated event TPM is inferred."""
import math

def bounded_fraction(fusion_fragments, unfused_fragments):
    """No pseudocount. None means no informative denominator, not zero."""
    if fusion_fragments is None or unfused_fragments is None:return None
    f,u=float(fusion_fragments),float(unfused_fragments)
    if not all(math.isfinite(x) and x>=0 for x in (f,u)):
        raise ValueError('Fragment counts must be finite and nonnegative')
    return f/(f+u) if f+u else None

def expression_proxy(supplied_gene_expression, event_fraction=1):
    """Keep missing distinct from zero; accept bounded fractions only."""
    if supplied_gene_expression is None or event_fraction is None:return None
    e,f=float(supplied_gene_expression),float(event_fraction)
    if not (math.isfinite(e) and e>=0 and math.isfinite(f) and 0<=f<=1):
        raise ValueError('Expected nonnegative expression and a fraction in [0,1]')
    return e*f

def fusion_proxies(left_expression,right_expression,fusion_fragments,left_unfused,right_unfused):
    """Retain both estimates, using the native fusion-support count definition."""
    lf=bounded_fraction(fusion_fragments,left_unfused)
    rf=bounded_fraction(fusion_fragments,right_unfused)
    return {'left_fraction':lf,'right_fraction':rf,
            'left_expression_proxy':expression_proxy(left_expression,lf),
            'right_expression_proxy':expression_proxy(right_expression,rf)}

def event_expression_proxies(source,left_expression,right_expression,left_fraction,right_fraction):
    """Splice ranking uses its assigned gene; fractions apply only to fusions."""
    if source=='SPLICE':
        value=expression_proxy(left_expression)
        return value,value
    if source=='CTA/SELF':
        return expression_proxy(left_expression),expression_proxy(right_expression)
    if source=='FUSION':
        return expression_proxy(left_expression,left_fraction),expression_proxy(right_expression,right_fraction)
    raise ValueError('Unknown antigen source: '+source)
