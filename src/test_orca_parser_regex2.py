import re

_POST_COORDS_KEYWORDS = frozenset({'%frag', '%eprnmr', '%nmr', '%rel', '%epr'})

def _split_post_blocks_regex(text: str) -> tuple[str, str]:
    if not text.strip():
        return text, ""
    
    lines = text.split('\n')
    clean_lines = [line.split('#')[0] for line in lines]
    clean_text = '\n'.join(clean_lines)
    
    pre = []
    post = []
    
    # Split before every %, !, or *
    parts = re.split(r'(?=[%!*])', clean_text)
    
    for part in parts:
        stripped = part.strip()
        if not stripped:
            continue
            
        # Get the first token
        first_token = stripped.split()[0].lower()
        
        if first_token in _POST_COORDS_KEYWORDS:
            post.append(stripped)
        else:
            pre.append(stripped)
            
    return "\n".join(pre), "\n".join(post)

text = "! optts %maxcore 4000 %eprnmr nuclei H = all {shift} end end \n! Freq \n%pal nprocs_group 4 end %scf MaxIter 200 end"
pre, post = _split_post_blocks_regex(text)
print("PRE:\n" + pre)
print("POST:\n" + post)
