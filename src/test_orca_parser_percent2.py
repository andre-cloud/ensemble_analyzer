_POST_COORDS_KEYWORDS = frozenset({'%frag', '%eprnmr', '%nmr', '%rel', '%epr'})

def _split_post_blocks_percent(text: str) -> tuple[str, str]:
    if not text.strip():
        return text, ""
    
    # Remove comments
    lines = text.split('\n')
    no_comment_lines = [line.split('#')[0] for line in lines]
    clean_text = '\n'.join(no_comment_lines)
    
    pre = []
    post = []
    
    parts = clean_text.split('%')
    pre.append(parts[0])
    
    for part in parts[1:]:
        first_word = part.split()[0].lower() if part.split() else ""
        kw = f"%{first_word}"
        
        if kw in _POST_COORDS_KEYWORDS:
            post.append(f"%{part}")
        else:
            pre.append(f"%{part}")
            
    return "".join(pre), "".join(post)

text1 = "%scf MaxIter 200 end # ignoring %frag\n%maxcore 4000 %eprnmr nuclei H = all {shift} end end"

print("--- TEXT 1 ---")
pre, post = _split_post_blocks_percent(text1)
print("PRE:\n" + pre)
print("POST:\n" + post)
