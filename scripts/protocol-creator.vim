" protocol-creator.vim -- Vim tools for ensemble_analyzer protocol JSON
"
" Install: source this from vimrc or symlink into ~/.vim/plugin/
"
"   :ProtocolNew             Start a new protocol from scratch
"   :ProtocolStep [N]        Insert generic step N (default: auto-number)
"   :ProtocolStepSP [N]      Insert single-point step
"   :ProtocolStepOpt [N]     Insert optimisation step
"   :ProtocolStepFreq [N]    Insert frequency step
"   :ProtocolStepOptFreq [N] Insert optimise+freq step
"   :ProtocolHelp            Show this message
"
"   <C-x><C-u>  Omni-completion (field names, functionals, basis sets)
"   K           Show field documentation (Normal mode)
"   <C-k>       Show field documentation (Insert mode)

if exists('g:protocol_creator_loaded')
  finish
endif
let g:protocol_creator_loaded = 1

" ─── Config ────────────────────────────────────────────────────────
if !exists('g:protocol_creator_project_root')
  let g:protocol_creator_project_root = ''
endif

" ─── Field docs ────────────────────────────────────────────────────
let s:field_docs = {
      \ 'functional':            'XC functional, e.g. B3LYP, r2SCAN-3c',
      \ 'basis':                 'Basis set, e.g. def2-SVP, def2-QZVPP',
      \ 'calculator':            'Program: orca (default), gaussian, tblite, aimnet',
      \ 'solvent':               'Dict: {"solvent": "water", "smd": false}',
      \ 'mult':                  'Multiplicity (default 1)',
      \ 'charge':                'Charge (default 0)',
      \ 'opt':                   'true = geometry optimisation',
      \ 'freq':                  'true = frequency calculation',
      \ 'freq_fact':             'Frequency scaling factor (default 1.0)',
      \ 'constrains':            'Geom constraints: [[a,b]], [[a,b,c]], [[a,b,c,d]], [[a,b,c,d,e]]',
      \ 'read_orbitals':         'Read MOs from step N',
      \ 'add_input':             'Raw string injected into calculator input',
      \ 'graph':                 'true = generate graph pickles',
      \ 'no_prune':              'true = disable pruning for this step',
      \ 'cluster':               'false/true or int(N clusters)',
      \ 'thrG':                  'Energy threshold (kcal/mol)',
      \ 'thrB':                  'Frequency threshold (cm^{-1})',
      \ 'thrGMAX':               'Max energy threshold (kcal/mol)',
      \ 'monitor_internals':     'Track internal coords [[i,j], [i,j,k], ...]',
      \ 'comment':               'Step description',
      \ 'read_population':       'Read Boltzmann pop from step N',
      \ 'skip_opt_fail':         'true = continue if opt fails',
      \ 'block_on_retention_rate': 'true = error if retention < 20%',
      \ }

" ─── Find project root ─────────────────────────────────────────────
function! s:FindProjectRoot(...)
  if g:protocol_creator_project_root !=# '' && isdirectory(g:protocol_creator_project_root)
    return g:protocol_creator_project_root
  endif
  let dir = fnamemodify(a:0 ? a:1 : expand('%:p:h'), ':p')
  let depth = 0
  while depth < 20
    if filereadable(dir . 'setup.py') || filereadable(dir . 'pyproject.toml') || isdirectory(dir . 'src/ensemble_analyzer')
      return dir
    endif
    let next = fnamemodify(dir, ':h')
    if next == dir | break | endif
    let dir = next
    let depth += 1
  endwhile
  return ''
endfunction

" ─── Load functionals & basis sets ─────────────────────────────────
let s:functionals = []
let s:basis_sets = []

function! s:LoadProjectData()
  let root = s:FindProjectRoot()
  if root ==# '' | return | endif
  let func_file = root . 'src/ensemble_analyzer/parameters_file/functionals'
  if filereadable(func_file)
    for line in readfile(func_file)
      let s = substitute(substitute(line, '#.*', '', ''), '^\s*\|\s*$', '', 'g')
      if s !=# '' && s !~ '^#'
        call add(s:functionals, s)
      endif
    endfor
  endif
  let basis_file = root . 'src/ensemble_analyzer/parameters_file/basis_sets'
  if filereadable(basis_file)
    for line in readfile(basis_file)
      let s = substitute(substitute(line, '#.*', '', ''), '^\s*\|\s*$', '', 'g')
      if s !=# '' && s !~ '^#'
        call add(s:basis_sets, substitute(s, '\s*|.*', '', ''))
      endif
    endfor
  endif
endfunction
call s:LoadProjectData()

if empty(s:functionals)
  let s:functionals = [
        \ 'B3LYP', 'B3LYP/G', 'PBE0', 'BP86', 'BLYP', 'PBE', 'TPSS',
        \ 'rSCAN', 'SCAN', 'M06', 'M062X', 'M06L', 'M06L', 'MN15',
        \ 'TPSSH', 'SCAN0', 'r2SCAN0', 'PW6B95',
        \ 'wB97X', 'wB97X-D3BJ', 'wB97X-D4', 'wB97X-D4rev',
        \ 'CAM-B3LYP', 'B2PLYP', 'DSD-PBEP86',
        \ 'B97-D3', 'B97-D4', 'B97M-D4', 'wB97M-D4',
        \ 'HF-3c', 'B97-3c', 'r2SCAN-3c', 'PBEh-3c', 'wB97X-3c',
        \ 'r2SCAN', 'r2SCANh',
        \ ]
endif
if empty(s:basis_sets)
  let s:basis_sets = [
        \ 'def2-SVP', 'def2-SVPD', 'def2-TZVP', 'def2-TZVPP',
        \ 'def2-TZVP(-f)', 'def2-mTZVP', 'def2-mTZVPP',
        \ 'def2-QZVP', 'def2-QZVPP',
        \ 'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ',
        \ 'aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ',
        \ 'SVP', 'TZVP', 'TZVPP', 'QZVP',
        \ '6-31G*', '6-31+G**', '6-311G*', '6-311+G*',
        \ 'pcseg-0', 'pcseg-1', 'pcseg-2',
        \ 'DKH-def2-TZVP', 'DKH-def2-TZVPP', 'DKH-def2-QZVPP',
        \ 'ZORA-def2-TZVP', 'ZORA-def2-TZVPP',
        \ 'x2c-TZVPall', 'x2c-QZVPall',
        \ 'ANO-pVDZ', 'ANO-pVTZ',
        \ 'vDZP', 'MINIX', 'STO-3G',
        \ ]
endif

" ─── Omni-complete ─────────────────────────────────────────────────
function! s:CompleteProtocol(findstart, base)
  if a:findstart
    let line = getline('.')
    let start = col('.') - 1
    while start > 0 && line[start - 1] =~ '[a-zA-Z_":]'
      let start -= 1
    endwhile
    return start
  endif
  let candidates = []
  let base = a:base
  " Detect key being completed by looking backward for "key":
  let before = strpart(getline('.'), 0, col('.') - 1)
  let context_key = substitute(matchstr(before, '"\zs\w\+\ze"\s*:\s*"\{0,1}$'), '"', '', 'g')
  if context_key ==# 'functional'
    for f in s:functionals
      if f =~? '^' . base
        call add(candidates, {'word': f, 'menu': '(functional)', 'dup': 1})
      endif
    endfor
    for f in s:functionals
      if f =~? '^' . base
        for disp in [' D3', ' D3BJ', ' D4']
          call add(candidates, {'word': f . disp, 'menu': '(func+disp)', 'dup': 1})
        endfor
        break
      endif
    endfor
  elseif context_key ==# 'basis'
    for b in s:basis_sets
      if b =~? '^' . base
        call add(candidates, {'word': b, 'menu': '(basis)', 'dup': 1})
      endif
    endfor
  elseif context_key ==# 'calculator'
    for c in ['orca', 'gaussian', 'tblite', 'aimnet']
      if c =~? '^' . base
        call add(candidates, {'word': c, 'menu': '(calc)', 'dup': 1})
      endif
    endfor
  elseif context_key ==# 'mult' || context_key ==# 'cluster' || context_key ==# 'read_orbitals' || context_key ==# 'read_population'
    for n in ['1', '2', '3', '4', '5', '6', '7', '8', '9']
      if n =~ '^' . base
        call add(candidates, {'word': n, 'menu': '(int)', 'dup': 1})
      endif
    endfor
  elseif context_key ==# 'charge'
    for c in ['0', '1', '-1', '2', '-2']
      if c =~ '^' . base
        call add(candidates, {'word': c, 'menu': '(charge)', 'dup': 1})
      endif
    endfor
  elseif context_key ==# 'opt' || context_key ==# 'freq' || context_key ==# 'graph' || context_key ==# 'no_prune' || context_key ==# 'skip_opt_fail'
    for v in ['true', 'false']
      if v =~ '^' . base
        call add(candidates, {'word': v, 'menu': '(bool)', 'dup': 1})
      endif
    endfor
  elseif context_key ==# 'solvent'
    for s in ['water', 'chloroform', 'acetonitrile', 'methanol', 'ethanol', 'THF', 'toluene', 'benzene', 'DMSO', 'DMF', 'acetone']
      if s =~? '^' . base
        call add(candidates, {'word': '"' . s . '"', 'menu': '(solvent)', 'dup': 1})
      endif
    endfor
  elseif context_key ==# 'smd'
    for v in ['true', 'false']
      if v =~ '^' . base
        call add(candidates, {'word': v, 'menu': '(bool)', 'dup': 1})
      endif
    endfor
  else
    for k in keys(s:field_docs)
      if k =~ '^' . base
        call add(candidates, {'word': '"' . k . '": ', 'menu': s:field_docs[k], 'dup': 1})
      endif
    endfor
  endif
  return candidates
endfunction

" ─── Syntax highlighting ───────────────────────────────────────────
function! s:SetupProtocolSyntax()
  setlocal omnifunc=s:CompleteProtocol
  if !exists('g:syntax_on') | syntax enable | endif
  syn keyword eaProtocolField contained
        \ functional basis calculator solvent
        \ mult charge opt freq freq_fact
        \ constrains read_orbitals add_input
        \ graph no_prune cluster
        \ thrG thrB thrGMAX
        \ monitor_internals comment
        \ read_population skip_opt_fail block_on_retention_rate
  hi def link eaProtocolField Special
  syn match eaStepNumber /"\d\+"\s*:/ containedin=ALLBUT,jsonString
  hi def link eaStepNumber Label
endfunction

augroup protocol_creator_auto
  au!
  au BufRead,BufNewFile *protocol*.json call s:SetupProtocolSyntax()
  au BufRead,BufNewFile *.json
        \ if expand('<afile>:p') =~ 'ensemble_analyzer'
        \ \| call s:SetupProtocolSyntax() \| endif
augroup END

" ─── Step insert ───────────────────────────────────────────────────
function! s:NextStepNumber()
  let max = -1
  for line in getline(1, '$')
    let m = matchlist(line, '"\(\d\+\)"\s*:')
    if len(m) > 1
      let n = str2nr(m[1])
      if n > max | let max = n | endif
    endif
  endfor
  return max + 1
endfunction

function! s:InsertStep(step_num, fields)
  let cur = getline('.')
  let top_indent = '    '
  let ins_line = line('.')
  let do_comma = 0

  " Cursor on top-level closing } → insert before it (add at end)
  if cur =~ '^\s*}$'
    let ins_line = line('.') - 1
    let top_indent = matchstr(cur, '^\s*')
    if getline(ins_line) !~ ',$' && getline(ins_line) !~ '^\s*{'
      let do_comma = 1
    endif
  " Cursor on step closing }, → insert after it
  elseif cur =~ '^\s*},$'
    let ins_line = line('.')
    let top_indent = matchstr(cur, '^\s*')
  " Cursor on a step header "N": { → find its closing and insert after
  elseif cur =~ '^\s*"\d\+"\s*:\s*{'
    let cnt = 1
    let l = line('.') + 1
    while cnt > 0 && l <= line('$')
      let ll = getline(l)
      let cnt += count(ll, '{') - count(ll, '}')
      let l += 1
    endwhile
    let ins_line = l - 1
    let top_indent = matchstr(getline(ins_line), '^\s*')
    if getline(ins_line) =~ '^\s*}$'
      call setline(ins_line, substitute(getline(ins_line), '}$', '},', ''))
    endif
  else
    let top_indent = matchstr(cur, '^\s*')
  endif
  let top_indent = top_indent ==# '' ? '    ' : top_indent
  let inner = top_indent . '    '

  " Add comma to previous line if needed
  if do_comma
    call setline(ins_line, getline(ins_line) . ',')
  endif

  " Determine if this will be the last step (before closing })
  let next_nonblank = nextnonblank(ins_line + 2)
  let is_last = (next_nonblank == 0 || getline(next_nonblank) =~ '^\s*}$')
  let step_lines = [top_indent . '"' . a:step_num . '": {']
  for [k, v] in a:fields
    call add(step_lines, inner . '"' . k . '": ' . v . ',')
  endfor
  if len(step_lines) > 1
    let step_lines[-1] = substitute(step_lines[-1], ',$', '', '')
  endif
  if is_last
    call add(step_lines, top_indent . '}')
  else
    call add(step_lines, top_indent . '},')
  endif
  call append(ins_line, step_lines)
endfunction

function! s:Step(...)
  let n = a:0 ? a:1 : s:NextStepNumber()
  call s:InsertStep(n, [
        \ ['functional', '"B3LYP"'],
        \ ['basis', '"def2-SVP"'],
        \ ['opt', 'false'],
        \ ['freq', 'false'],
        \ ['mult', '1'],
        \ ['charge', '0'],
        \ ])
endfunction

function! s:StepSP(...)
  let n = a:0 ? a:1 : s:NextStepNumber()
  call s:InsertStep(n, [
        \ ['functional', '"B3LYP"'],
        \ ['basis', '"def2-TZVP"'],
        \ ['mult', '1'],
        \ ['charge', '0'],
        \ ])
endfunction

function! s:StepOpt(...)
  let n = a:0 ? a:1 : s:NextStepNumber()
  call s:InsertStep(n, [
        \ ['functional', '"r2SCAN-3c"'],
        \ ['opt', 'true'],
        \ ['mult', '1'],
        \ ['charge', '0'],
        \ ])
endfunction

function! s:StepFreq(...)
  let n = a:0 ? a:1 : s:NextStepNumber()
  call s:InsertStep(n, [
        \ ['functional', '"r2SCAN-3c"'],
        \ ['freq', 'true'],
        \ ['mult', '1'],
        \ ['charge', '0'],
        \ ])
endfunction

function! s:StepOptFreq(...)
  let n = a:0 ? a:1 : s:NextStepNumber()
  call s:InsertStep(n, [
        \ ['functional', '"r2SCAN-3c"'],
        \ ['opt', 'true'],
        \ ['freq', 'true'],
        \ ['mult', '1'],
        \ ['charge', '0'],
        \ ])
endfunction

function! s:NewProtocol()
  set ft=json
  call setline(1, [
        \ '{',
        \ '    "0": {',
        \ '        "functional": "b97-3c",',
        \ '        "basis": "def2-mTZVP"',
        \ '    },',
        \ '    "1": {',
        \ '        "functional": "r2scan-3c",',
        \ '        "basis": "def2-mTZVPP",',
        \ '        "opt": true,',
        \ '        "freq": true',
        \ '    },',
        \ '    "2": {',
        \ '        "functional": "wb97x-d4rev",',
        \ '        "basis": "def2-qzvpp"',
        \ '    }',
        \ '}',
        \ ])
  call s:SetupProtocolSyntax()
endfunction

" ─── Commands ──────────────────────────────────────────────────────
command! -nargs=? ProtocolNew        call s:NewProtocol()
command! -nargs=? ProtocolStep       call s:Step(<f-args>)
command! -nargs=? ProtocolStepSP     call s:StepSP(<f-args>)
command! -nargs=? ProtocolStepOpt    call s:StepOpt(<f-args>)
command! -nargs=? ProtocolStepFreq   call s:StepFreq(<f-args>)
command! -nargs=? ProtocolStepOptFreq call s:StepOptFreq(<f-args>)
command! -nargs=0 ProtocolHelp       call s:ShowHelp()

" ─── K mapping ─────────────────────────────────────────────────────
function! s:ShowFieldDoc()
  let word = substitute(expand('<cword>'), '^"\|"$', '', 'g')
  if has_key(s:field_docs, word)
    echohl Title | echo word . ': ' . s:field_docs[word] | echohl None
  elseif word !=# ''
    echo 'Unknown field: ' . word
  endif
endfunction

function! s:ShowHelp()
  echohl Title
  echo 'protocol-creator.vim  --  :ProtocolHelp  for help'
  echohl None
  echo '  :ProtocolNew         create new protocol'
  echo '  :ProtocolStep [N]    insert step (auto-number)'
  echo '  :ProtocolStepSP [N]  insert single-point step'
  echo '  :ProtocolStepOpt [N] insert optimisation step'
  echo '  :ProtocolStepFreq [N] insert frequency step'
  echo '  :ProtocolStepOptFreq [N] insert opt+freq step'
  echo ''
  echo '  <C-x><C-u>  omni-completion for field names, functionals, basis sets'
  echo '  K           show field documentation (Normal mode)'
  echo '  <C-k>       show field documentation (Insert mode)'
endfunction

augroup protocol_creator_maps
  au!
  au FileType json nnoremap <buffer> <silent> K :call <SID>ShowFieldDoc()<CR>
  au FileType json inoremap <buffer> <silent> <C-k> <C-o>:call <SID>ShowFieldDoc()<CR>
augroup END

" ─── Info ──────────────────────────────────────────────────────────
echohl Comment
echomsg 'protocol-creator.vim loaded. :ProtocolHelp for usage'
echohl None