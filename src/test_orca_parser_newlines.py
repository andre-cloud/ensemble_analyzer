_POST_COORDS_KEYWORDS = frozenset({'%frag', '%eprnmr', '%nmr', '%rel', '%epr'})

def _split_post_blocks_percent(text: str) -> tuple[str, str]:
    if not text.strip():
        return text, ""
    
    lines = text.split('\n')
    clean_lines = [line.split('#')[0] for line in lines]
    clean_text = '\n'.join(clean_lines)
    
    pre = []
    post = []
    
    parts = clean_text.split('%')
    
    # We strip parts to avoid accumulating spaces, and ensure clean newlines
    if parts[0].strip():
        # Handle the possibility of multiple '!' lines glued together without newlines?
        # A simple approach is just replace '!' with '\n!' if it's not the start, 
        # but let's just leave part 0 as is, just stripped
        # Actually, replace '! ' with '\n! ' if we want to force newlines?
        p0 = parts[0].strip()
        # To ensure ! starts on a new line if it's glued to something? But it's part 0, so it's the beginning!
        pre.append(p0)
    
    for part in parts[1:]:
        stripped = part.strip()
        if not stripped: continue
        
        first_word = stripped.split()[0].lower()
        kw = f"%{first_word}"
        
        if kw in _POST_COORDS_KEYWORDS:
            post.append(f"%{stripped}")
        else:
            pre.append(f"%{stripped}")
            
    return "\n".join(pre), "\n".join(post)

text = "! optts %maxcore 4000 %eprnmr nuclei H = all {shift} end end \n! Freq \n%pal nprocs_group 4 end %scf MaxIter 200 end"
pre, post = _split_post_blocks_percent(text)
print("PRE:")
print(pre)
print("POST:")
print(post)
