_POST_COORDS_KEYWORDS = frozenset({'%frag', '%eprnmr', '%nmr', '%rel', '%epr'})

def _split_post_blocks_percent(text: str) -> tuple[str, str]:
    if not text.strip():
        return text, ""
    
    # We should handle comments carefully, but for now let's assume no % in comments
    # or we can remove comments first? No, we might want to keep comments.
    
    pre = []
    post = []
    
    # Simple lines that don't start with % at all (like ! NMR) usually come before any % block
    # But wait, splitting by % breaks "! NMR \n %scf" into "! NMR \n " and "scf"
    parts = text.split('%')
    
    # The first part is everything before the first % (e.g. ! NMR, * xyz, etc.)
    pre.append(parts[0])
    
    for part in parts[1:]:
        # Find the first word
        first_word = part.split()[0].lower() if part.split() else ""
        kw = f"%{first_word}"
        
        if kw in _POST_COORDS_KEYWORDS:
            post.append(f"%{part}")
        else:
            pre.append(f"%{part}")
            
    return "".join(pre), "".join(post)

text1 = "%maxcore 4000 %eprnmr nuclei H = all {shift} end end %pal nprocs_group 4 end %scf MaxIter 200 end"
text2 = "! NMR\n%eprnmr\n nuclei H = all {shift} end\nend\n%geom end"

for i, t in enumerate([text1, text2]):
    print(f"--- TEXT {i+1} ---")
    pre, post = _split_post_blocks_percent(t)
    print("PRE:\n" + pre)
    print("POST:\n" + post)
    print()
